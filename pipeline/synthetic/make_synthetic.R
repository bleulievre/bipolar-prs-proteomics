#!/usr/bin/env Rscript
# Simulated data for testing the pipeline end to end, written in two layouts from the same simulated people:
#   synthetic/glasgow/  the layout of the dataset used for the dissertation (UKB.csv as written by UKB.R, with
#                       proteins in columns 34+, sex_plink coded 1/2, withdrawn participants as negative ids,
#                       and a separate anyMentalIllness file), with a data map derived from config/data_map_glasgow.yaml
#   synthetic/rap/      the UK Biobank RAP template layout of config/data_map.yaml (separate participant, PRS and
#                       protein files, withdrawal list, kinship file, assay manifest, code lists)
# Nothing here is derived from UK Biobank data: every value is drawn from the distributions below.
# Planted truth (truth.csv): 8 assays depend on the BD-PRS (7 inverse, 1 positive), one also through eGFR;
# the others are null. Missingness varies by assay so that some exceed 10%.
#
#   Rscript synthetic/make_synthetic.R [n_participants] [n_assays] [output_dir]

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args) >= 1) as.integer(args[1]) else 8000L
k <- if (length(args) >= 2) as.integer(args[2]) else 60L
here <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
repo <- dirname(here)
out <- if (length(args) >= 3) args[3] else here
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(20261020)

eid <- 1000000L + seq_len(n)
draw <- function(codes, p) sample(codes, n, replace = TRUE, prob = p)
sex01 <- draw(0:1, c(0.54, 0.46))                     # 0 female, 1 male (field 31 coding)
age <- round(runif(n, 40, 70))
pcs <- matrix(rnorm(n * 8), n, 8)
batch <- ifelse(runif(n) < 0.1, -sample(1:11, n, TRUE), sample(1:95, n, TRUE))
townsend <- rnorm(n, -1.4, 3)
prs_raw <- rnorm(n, 0.18, 0.14)
prs_std <- (prs_raw - mean(prs_raw)) / sd(prs_raw)
psych <- rbinom(n, 1, plogis(-2.4 + 0.25 * prs_std))
smoking <- draw(c(0, 1, 2, -3), c(0.55, 0.34, 0.10, 0.01))
alcohol <- draw(c(0, 1, 2, -3), c(0.04, 0.04, 0.91, 0.01))
bmi <- rnorm(n, 27.4, 4.8)
hip <- rnorm(n, 103, 8)
waist <- hip * rnorm(n, 0.87, 0.08)
sbp <- rnorm(n, 138, 18); dbp <- rnorm(n, 82, 10)
bpmed <- rbinom(n, 1, plogis(-1.5 + 0.03 * (age - 55)))
sbp_adj <- sbp + 15 * bpmed; dbp_adj <- dbp + 10 * bpmed
hyper <- rbinom(n, 1, 0.29); lipid <- rbinom(n, 1, 0.19)
stroke <- rbinom(n, 1, 0.02); ihd <- rbinom(n, 1, 0.07)
stroke[runif(n) < 0.05] <- NA; ihd[runif(n) < 0.05] <- NA
t2d <- rbinom(n, 1, plogis(-3 + 0.12 * (bmi - 27)))
egfr <- rnorm(n, 90, 14) - 0.6 * (age - 55) + 1.5 * prs_std   # a small PRS -> filtration path (cf. Sewell 2021)
ancestry <- draw(c(1001, 1002, 1003, 3001, 4001, -3), c(0.90, 0.03, 0.03, 0.02, 0.015, 0.005))
ad_use <- rbinom(n, 1, ifelse(psych == 1, 0.45, 0.05))
ap_use <- rbinom(n, 1, ifelse(psych == 1, 0.10, 0.003))
prs_missing <- runif(n) < 0.09                                   # ~9% missing, as in the dissertation
townsend[runif(n) < 0.02] <- NA                                  # missing covariate (core in the RAP layout)

assays <- sprintf("assay%03d", seq_len(k))
planted <- data.table(assay = assays[1:8], beta_per_sd = c(-0.06, -0.05, -0.05, -0.04, -0.04, -0.035, -0.03, 0.05),
                      via_egfr = c(TRUE, rep(FALSE, 7)))
prot <- matrix(NA_real_, n, k, dimnames = list(NULL, assays))
for (j in seq_len(k)) {
  a <- assays[j]
  mu <- 0.01 * (age - 55) + 0.1 * sex01 + 0.02 * pcs[, 1] + 0.01 * (bmi - 27)
  b <- if (a %in% planted$assay) planted[assay == a, beta_per_sd] else 0
  y <- mu + b * prs_std + rnorm(n, 0, 0.5)
  if (a %in% planted[via_egfr == TRUE, assay]) y <- y - 0.02 * (egfr - 90)
  miss <- if (a %in% planted$assay) 0.02 else c(0.01, 0.03, 0.05, 0.08, 0.10, 0.15, 0.25, 0.45)[(j %% 8) + 1]
  y[runif(n) < miss] <- NA
  prot[, j] <- round(y, 4)
}

# ---------------- Glasgow layout (UKB.csv as written by write.csv in UKB.R) ----------------
g <- file.path(out, "glasgow")
dir.create(file.path(g, "data"), recursive = TRUE, showWarnings = FALSE)
withdrawn_rows <- sample.int(n, 5)
n_eid <- eid
n_eid[withdrawn_rows] <- -(seq_along(withdrawn_rows))            # UK Biobank marks withdrawals with negative ids
lead <- data.table(V1 = seq_len(n), n_eid = n_eid,               # V1 = the row-name column that write.csv adds
                   prs_bp = ifelse(prs_missing, NA, prs_raw), prs_bp_std = ifelse(prs_missing, NA, prs_std),
                   age_0_0 = age, sex_plink = ifelse(sex01 == 1, 1, 2))
for (i in 1:8) lead[, (paste0("PGC", i)) := pcs[, i]]
lead[, `:=`(chip = as.integer(batch > 0), townsend = townsend, n_21001_0_0 = bmi, whr_0 = waist / hip,
            alcoholstatus = alcohol, smokestatus = smoking, sbp_0 = sbp, dbp_0 = dbp, stroke_0 = stroke, ht_0 = hyper,
            ish_0 = ihd, lipidlowering_0 = lipid, htnmed_0 = bpmed, adj_sbp_0 = sbp_adj, adj_dbp_0 = dbp_adj,
            adj_probposs_t2dm_0 = t2d, any_antidepBL = ad_use, any_antipsychBL = ap_use, n_31_0_0 = sex01)]
stopifnot(ncol(lead) == 33)                                       # proteins start at column 34, as in UKB.R
ukb <- cbind(lead, as.data.table(prot), data.table(prs_scz_std = 0.3 * prs_std + rnorm(n, 0, 0.95)))
fwrite(ukb, file.path(g, "data", "UKB.csv"))
mental <- data.table(n_eid = n_eid, anyMentalIllness = ifelse(runif(n) < 0.05, NA, psych))
mental <- mental[runif(n) > 0.03]                                 # some participants absent from the supplement
fwrite(mental, file.path(g, "data", "supplement.csv"))
gmap <- yaml::read_yaml(file.path(repo, "config", "data_map_glasgow.yaml"))
gmap$proteins$columns <- list(from = 34L, to = 33L + k, first = assays[1], last = assays[k])   # no check_names
gmap$variables$egfr <- NULL
writeLines(c("# Generated by synthetic/make_synthetic.R from config/data_map_glasgow.yaml (assay columns 34 to 33 + k).",
             yaml::as.yaml(gmap)), file.path(g, "data_map.yaml"))
fwrite(planted, file.path(g, "truth.csv"))

# ---------------- RAP template layout ----------------
r <- file.path(out, "rap")
dir.create(file.path(r, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(r, "config", "codes"), recursive = TRUE, showWarnings = FALSE)
multi <- function(prefix, width, fill) {
  m <- matrix(NA_real_, n, width, dimnames = list(NULL, paste0(prefix, 0:(width - 1))))
  for (i in seq_len(n)) {
    v <- fill(i)
    if (length(v)) m[i, seq_along(v)] <- v
  }
  m
}
vasc <- multi("p6150_i0_a", 4, function(i) {
  if (runif(1) < 0.01) return(-3)
  v <- c(if (isTRUE(ihd[i] == 1)) 2, if (isTRUE(stroke[i] == 1)) 3, if (hyper[i] == 1) 4)
  if (length(v)) v else -7
})
AD <- c(9100001, 9100002, 9100003); AP <- c(9200001, 9200002); PSY <- c(9300001, 9300002, 9300003)  # not real codes
meds <- multi("p20003_i0_a", 5, function(i)
  c(if (ad_use[i]) sample(AD, 1), if (ap_use[i]) sample(AP, 1), if (runif(1) < 0.5) sample(9900001:9900050, 1)))
ill <- multi("p20002_i0_a", 4, function(i) c(if (psych[i]) sample(PSY, 1), if (runif(1) < 0.3) sample(9800001:9800040, 1)))
pc_cols <- pcs; colnames(pc_cols) <- paste0("p22009_a", 1:8)
participants <- data.table(eid, p31 = sex01, p21022 = age, pc_cols, p22000 = batch, p189 = townsend,
                           p20116_i0 = smoking, p20117_i0 = alcohol, p21001_i0 = bmi, p48_i0 = waist, p49_i0 = hip,
                           sbp_mean_i0 = sbp, dbp_mean_i0 = dbp, sbp_adj_i0 = sbp_adj, dbp_adj_i0 = dbp_adj,
                           hypertension_i0 = hyper, lipid_lowering_i0 = lipid, bp_medication_i0 = bpmed,
                           t2d_eastwood_probable_possible = t2d, vasc, meds, ill, p21000_i0 = ancestry, egfr_ckdepi = egfr)
fwrite(participants, file.path(r, "data", "participants.tsv"), sep = "\t")
fwrite(data.table(eid, prs_bd = ifelse(prs_missing, NA, prs_std), prs_scz = 0.3 * prs_std + rnorm(n, 0, 0.95)),
       file.path(r, "data", "prs.tsv"), sep = "\t")
fwrite(cbind(data.table(eid), as.data.table(prot)), file.path(r, "data", "olink_instance0.tsv"), sep = "\t")
pairs <- matrix(sample(eid, 300), ncol = 2)
fwrite(data.table(ID1 = pairs[, 1], ID2 = pairs[, 2], HetHet = 0.1, IBS0 = 0.001, Kinship = runif(150, 0.05, 0.27)),
       file.path(r, "data", "ukb_rel.dat"), sep = " ")
writeLines(as.character(eid[withdrawn_rows]), file.path(r, "data", "withdrawals.csv"))
writeLines(c("# synthetic codes", AD), file.path(r, "config", "codes", "antidepressants.txt"))
writeLines(c("# synthetic codes", AP), file.path(r, "config", "codes", "antipsychotics.txt"))
writeLines(c("# synthetic codes", PSY), file.path(r, "config", "codes", "psychiatric_history.txt"))
fwrite(data.table(column = assays, assay_id = assays, target = toupper(assays), uniprot = NA, panel = "Synthetic"),
       file.path(r, "config", "assay_manifest.csv"))
file.copy(file.path(repo, "config", "data_map.yaml"), file.path(r, "data_map.yaml"), overwrite = TRUE)
fwrite(planted, file.path(r, "truth.csv"))
message(sprintf("synthetic data: %d participants, %d assays (8 planted) -> %s/{glasgow,rap}", n, k, out))
