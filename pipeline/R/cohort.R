# Building the analysis samples (exclusions in the order set in the plan) and the protein QC.
#   dataset = participants left after the plan's dataset steps (withdrawals, ancestry) and relatedness pruning;
#             the original (dissertation) procedure is reproduced on this sample
#   cohort  = dataset minus participants with a missing BD-PRS or missing core covariates (primary analysis)

# Seeded randomized greedy removal: pairs are put in a canonical order (so the file's row order does not matter),
# shuffled with a fixed seed and RNG, and one member of each still-intact pair is removed at random. No related
# pair remains; this does not guarantee the largest possible sample.
drop_related <- function(ids, pairs, seed) {
  if (is.null(pairs) || !nrow(pairs)) return(character())
  a <- pmin(pairs$ID1, pairs$ID2); b <- pmax(pairs$ID1, pairs$ID2)
  p <- unique(data.table(ID1 = a, ID2 = b))[ID1 != ID2 & ID1 %in% ids & ID2 %in% ids]
  if (!nrow(p)) return(character())
  setorder(p, ID1, ID2)
  old_kind <- RNGkind()
  old_seed <- if (exists(".Random.seed", globalenv(), inherits = FALSE)) get(".Random.seed", globalenv()) else NULL
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (is.null(old_seed)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old_seed, globalenv())
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  p <- p[sample.int(nrow(p))]
  coin <- runif(nrow(p)) < 0.5
  removed <- character()
  for (i in seq_len(nrow(p))) {
    if (p$ID1[i] %in% removed || p$ID2[i] %in% removed) next
    removed <- c(removed, if (coin[i]) p$ID1[i] else p$ID2[i])
  }
  removed
}

# Returns list(data = cohort, dataset = dataset sample, flow, prs_scaling).
build_cohort <- function(inp, plan, map) {
  d <- copy(inp$data)
  core <- unlist(plan$covariates$core)
  exposure <- plan$exposure$variable
  dataset_steps <- c("withdrawn", "not_european", "related")
  flow <- list(data.table(step = "start", label = "Participants with Olink protein data",
                          n_excluded = NA_integer_, n_remaining = nrow(d)))
  dataset <- NULL
  removed_related <- character()
  for (s in plan$cohort$exclusions) {
    if (is.null(dataset) && !s$flag %in% dataset_steps) dataset <- copy(d)
    drop <- switch(s$flag,
      withdrawn = d$eid %in% inp$withdrawn,
      related = {
        removed_related <- drop_related(d$eid, inp$kinship, map$relatedness$seed)
        d$eid %in% removed_related
      },
      missing_prs = is.na(d[[exposure]]),
      missing_core = !complete.cases(d[, ..core]),
      {
        if (!s$flag %in% names(d)) stop("exclusion flag '", s$flag, "' is not a mapped variable", call. = FALSE)
        is.na(d[[s$flag]]) | d[[s$flag]] == 1   # NA counts as "exclude" (see the plan)
      }
    )
    flow[[length(flow) + 1]] <- data.table(step = s$flag, label = s$label, n_excluded = sum(drop),
                                           n_remaining = sum(!drop))
    d <- d[!drop]
  }
  if (is.null(dataset)) dataset <- copy(d)
  dataset <- dataset[!eid %in% removed_related]   # relatedness pruning applies to the dataset sample too
  if (!nrow(d)) stop("no participants left after the exclusions", call. = FALSE)

  scaling <- list(mean = mean(d[[exposure]]), sd = sd(d[[exposure]]), n = nrow(d), scale = plan$exposure$scale)
  if (plan$exposure$scale == "as_supplied" && (abs(scaling$sd - 1) > 0.05 || abs(scaling$mean) > 0.05))
    warning(sprintf("the supplied BD-PRS has mean %.3f and SD %.3f in the cohort: coefficients are per unit of the supplied score, not per SD",
                    scaling$mean, scaling$sd), call. = FALSE)
  d[, prs_bd_supplied := get(exposure)]   # kept for Table 1
  if (plan$exposure$scale == "standardise_in_cohort") {
    d[, (exposure) := (get(exposure) - scaling$mean) / scaling$sd]
    dataset[, (exposure) := (get(exposure) - scaling$mean) / scaling$sd]   # same scale in both samples
  } else if (plan$exposure$scale != "as_supplied") {
    stop("exposure$scale must be standardise_in_cohort or as_supplied", call. = FALSE)
  }
  flow[[length(flow) + 1]] <- data.table(step = "cohort", label = "Eligible analysis cohort",
                                         n_excluded = NA_integer_, n_remaining = nrow(d))
  list(data = d, dataset = dataset, flow = rbindlist(flow), prs_scaling = scaling)
}

# Missing fraction of each assay; assays above the threshold (strictly greater) are excluded, so an assay with
# exactly the threshold fraction is kept (as `<= 10` in the original code).
protein_qc <- function(d, assays, max_missing) {
  qc <- data.table(assay = assays,
                   n_missing = vapply(assays, function(a) sum(is.na(d[[a]])), integer(1)),
                   n_total = nrow(d))
  qc[, missing_fraction := n_missing / n_total]
  qc[, analysed := missing_fraction <= max_missing + 1e-12]
  qc[]
}
