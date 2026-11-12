# Writing the analysis package: aggregate results only, with small counts suppressed.
# These checks are a safety net, not a substitute for UK Biobank's own rules on exporting results.

TABLE1_VARS <- c(prs_bd_supplied = "BD-PRS (as supplied)", age = "Age at recruitment, years", sex = "Sex",
                 townsend = "Townsend deprivation index", bmi = "Body mass index, kg/m2", whr = "Waist-to-hip ratio",
                 smoking_current = "Current smoking", alcohol_current = "Current alcohol drinking",
                 sbp = "Systolic blood pressure, mmHg", dbp = "Diastolic blood pressure, mmHg",
                 hypertension = "Hypertension", stroke = "Stroke", ihd = "Ischaemic heart disease",
                 lipid_lowering = "Lipid-lowering medication", bp_medication = "Antihypertensive medication",
                 t2d = "Type 2 diabetes (probable or possible)",
                 antidepressant = "Antidepressant use", antipsychotic = "Antipsychotic use",
                 psych_history = "Self-reported mental illness (restricted-sample definition)")

# A count and its percentage are shown only if both the count and its complement reach min_cell.
count_cell <- function(k, denom, min_cell) {
  if ((k > 0 && k < min_cell) || (denom - k > 0 && denom - k < min_cell))
    return(sprintf("suppressed (a cell below %d)", min_cell))
  sprintf("%s (%.1f%%)", format(k, big.mark = ","), 100 * k / denom)
}
small <- function(n, min_cell) {
  if (is.na(n)) return("")
  if (n > 0 && n < min_cell) sprintf("<%d", min_cell) else format(n, big.mark = ",")
}

table1 <- function(d, min_cell) {
  rows <- list()
  for (v in intersect(names(TABLE1_VARS), names(d))) {
    x <- d[[v]]
    miss <- small(sum(is.na(x)), min_cell)
    denom <- sum(!is.na(x))
    if (is.factor(x)) {
      for (lv in levels(x))
        rows[[length(rows) + 1]] <- data.table(variable = TABLE1_VARS[[v]], level = lv,
          summary = count_cell(sum(x == lv, na.rm = TRUE), denom, min_cell), n_missing = miss)
    } else if (all(x %in% c(0, 1, NA))) {
      rows[[length(rows) + 1]] <- data.table(variable = TABLE1_VARS[[v]], level = "yes",
        summary = count_cell(sum(x == 1, na.rm = TRUE), denom, min_cell), n_missing = miss)
    } else {
      rows[[length(rows) + 1]] <- data.table(variable = TABLE1_VARS[[v]], level = "mean (SD)",
        summary = if (denom < min_cell) sprintf("suppressed (n < %d)", min_cell) else
          sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)), n_missing = miss)
    }
  }
  rbindlist(rows)
}

# Missingness of each model covariate in the eligible cohort (Table S1), with small counts suppressed.
covariate_missingness <- function(d, plan, min_cell) {
  vars <- unique(c(unlist(plan$covariates$core), unlist(plan$covariates$extended), plan$secondary$restricted$exclude_flag))
  rbindlist(lapply(intersect(vars, names(d)), function(v) {
    k <- sum(is.na(d[[v]]))
    data.table(variable = v, n_missing = small(k, min_cell),
               pct_missing = if (k > 0 && k < min_cell) NA_real_ else round(100 * k / nrow(d), 1))
  }))
}

# Every exported table must be an aggregate: no id column, no column of participant ids, and no more rows than
# assays x analyses.
assert_aggregate <- function(tables, n_assays, ids = character()) {
  for (nm in names(tables)) {
    t <- tables[[nm]]
    if (any(tolower(names(t)) %in% c("eid", "n_eid", "id", "participant", "f.eid")))
      stop("refusing to export '", nm, "': it has a participant id column", call. = FALSE)
    if (nrow(t) > 20 * max(n_assays, 1))
      stop("refusing to export '", nm, "': more rows than an aggregate table should have", call. = FALSE)
    if (length(ids) && any(vapply(t, function(col) any(as.character(col) %in% ids), logical(1))))
      stop("refusing to export '", nm, "': it contains participant ids", call. = FALSE)
  }
  invisible(TRUE)
}

result_keys <- function(cohort, qc, prim, sec, legacy, spear, plan) {
  d <- cohort$data
  sel <- prim[bonferroni == TRUE]
  ext_models <- unlist(plan$secondary$models)
  common <- if (nrow(sec)) sec[sample == "common"] else sec
  ext_common <- if (nrow(common)) common[model %in% ext_models] else common
  avail <- if (nrow(sec)) sec[sample == "available"] else sec
  restr <- if (nrow(avail)) avail[model == "restricted"] else avail
  same_all <- if (nrow(ext_common)) ext_common[, .(ok = all(same_sign)), by = assay][, sum(ok)] else NA
  female <- if ("sex" %in% names(d)) sum(d$sex == "female", na.rm = TRUE) else NA
  rs <- plan$secondary$restricted
  legacy_final <- if (is.null(legacy)) character() else
    legacy$table[stage == tail(legacy$counts$stage, 1) & passed == TRUE, assay]
  list(
    cohort_n = nrow(d), dataset_n = nrow(cohort$dataset),
    age_mean = mean(d$age), age_sd = sd(d$age),
    female_n = female, female_pct = 100 * female / nrow(d),
    prs_scaling = cohort$prs_scaling,
    n_assays_in_data = nrow(qc), n_analysed = sum(qc$analysed), n_excluded_missing = sum(!qc$analysed),
    basic_n_min = min(prim$n), basic_n_max = max(prim$n),
    bonferroni_threshold = prim$threshold[1],
    n_bonf = nrow(sel), n_fdr = sum(prim$fdr), n_inverse = sum(sel$estimate < 0),
    beta_min = if (nrow(sel)) min(sel$estimate) else NA, beta_max = if (nrow(sel)) max(sel$estimate) else NA,
    max_partial_r2 = if (nrow(sel)) 100 * max(sel$partial_r2) else NA,
    max_partial_r2_all = 100 * max(prim$partial_r2),
    # change with added covariates, on the common sample (same rows in every model)
    median_change = if (nrow(ext_common)) median(ext_common$abs_pct_change, na.rm = TRUE) else NA,
    change_range = if (nrow(ext_common)) range(ext_common$abs_pct_change, na.rm = TRUE) else NA,
    n_same_sign_all = same_all,
    common_n_min = if (nrow(common)) min(common$n) else NA, common_n_max = if (nrow(common)) max(common$n) else NA,
    restricted_n = sum(restricted_rows(d, rs$exclude_flag, rs$na_action)),
    restricted_same_sign = if (nrow(restr)) sum(restr$same_sign) else NA,
    restriction_change_median = if (nrow(restr)) median(restr$abs_pct_change, na.rm = TRUE) else NA,
    n_negative_rho = sum(spear$rho < 0), max_abs_rho = max(abs(spear$rho)),
    legacy_counts = if (is.null(legacy)) NULL else legacy$counts,
    n_legacy_overlap = if (is.null(legacy)) NA else length(intersect(legacy_final, sel$assay))
  )
}

git_capture <- function(root, args) {
  out <- tryCatch(suppressWarnings(system2("git", c("-C", root, args), stdout = TRUE, stderr = FALSE)),
                  error = function(e) structure(character(), status = 127L))
  status <- attr(out, "status")
  list(ok = is.null(status) || identical(as.integer(status), 0L), out = out)
}

repo_relative <- function(path, root) {
  p <- normalizePath(path, mustWork = TRUE); r <- normalizePath(root, mustWork = TRUE)
  prefix <- paste0(r, .Platform$file.sep)
  if (startsWith(p, prefix)) substring(p, nchar(prefix) + 1L) else NA_character_
}

# Files whose exact content a FROZEN result depends on: plan, data map, code lists, manifests, R code, environment.
freeze_files <- function(root, config_paths, map) {
  aux <- unlist(lapply(map$variables, function(v) v[["codes_1_file"]]))
  aux <- c(aux, map$proteins$columns$manifest, map$proteins$columns$check_names)
  f <- c(config_paths, file.path(root, aux[!is.na(aux)]), file.path(root, "environment.yml"),
         list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE), file.path(root, "run_all.R"))
  unique(f[file.exists(f)])
}

# FROZEN only if: plan status FROZEN, ISO frozen_date, frozen_by named, every decision resolved, every file the run
# depends on tracked by git, and a clean working tree. Everything else is DRAFT (or SYNTHETIC for simulated data).
provenance <- function(plan, root, synthetic, config_paths, map = list(), unresolved = character()) {
  head <- git_capture(root, c("rev-parse", "HEAD"))
  status <- git_capture(root, c("status", "--porcelain"))
  commit <- if (head$ok && length(head$out) == 1L) head$out else "nogit"
  clean <- head$ok && status$ok && !length(status$out)
  files <- freeze_files(root, config_paths, map)
  rel <- vapply(files, repo_relative, "", root = root)
  tracked <- length(rel) > 0 && !anyNA(rel) &&
    all(vapply(rel, function(p) git_capture(root, c("ls-files", "--error-unmatch", p))$ok, logical(1)))
  tag <- git_capture(root, c("describe", "--tags", "--exact-match", "HEAD"))
  date_ok <- !is.null(plan$frozen_date) && grepl("^\\d{4}-\\d{2}-\\d{2}$", format(plan$frozen_date)) &&
    !is.na(as.Date(format(plan$frozen_date), optional = TRUE))
  frozen <- identical(plan$status, "FROZEN") && date_ok && length(unlist(plan$frozen_by)) > 0 && clean && tracked &&
    !length(unresolved)
  kind <- if (synthetic) "SYNTHETIC" else if (frozen) "FROZEN" else "DRAFT"
  id <- sprintf("%s@%s", plan$plan_id, plan$version)
  line <- switch(kind,
    SYNTHETIC = sprintf("SYNTHETIC %s %s - simulated data, pipeline test only", id, substr(commit, 1, 12)),
    FROZEN = sprintf("FROZEN %s %s", id, commit),
    DRAFT = sprintf("DRAFT %s %s%s - plan not frozen, or code/config not committed; not for publication", id,
                    substr(commit, 1, 12), if (clean) "" else "-dirty"))
  list(kind = kind, line = line, commit = commit, tag = if (tag$ok) tag$out else NA, clean = clean,
       dependencies_tracked = tracked, unresolved_decisions = unresolved,
       plan_status = plan$status, frozen_date = if (is.null(plan$frozen_date)) NA else format(plan$frozen_date),
       frozen_by = unlist(plan$frozen_by), relatedness_seed = map$relatedness$seed,
       dependencies_sha256 = as.list(setNames(vapply(files, function(p) digest::digest(file = p, algo = "sha256"), ""),
                                              ifelse(is.na(rel), basename(files), rel))))
}

write_package <- function(out, res, plan, root, synthetic = FALSE, config_paths = character(), overwrite = FALSE,
                          map = list()) {
  if (dir.exists(out)) {
    if (length(list.files(out, all.files = TRUE, no.. = TRUE)) && !overwrite)
      stop("output directory is not empty: ", out, " (choose a new run id)", call. = FALSE)
    unlink(out, recursive = TRUE)
  }
  tmp <- paste0(out, ".partial")
  unlink(tmp, recursive = TRUE)
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  w <- function(x, f) fwrite(x, file.path(tmp, f))
  min_cell <- plan$export$min_cell
  prim <- res$primary
  avail <- if (nrow(res$secondary)) res$secondary[sample == "available"] else res$secondary
  # Figure 3: Basic and Extended 1-3 on the common sample (same rows, so differences reflect the added covariates);
  # the restricted sample on its own rows. Falls back to available-case rows if no common sample was fitted.
  fig_rows <- if (nrow(res$secondary)) {
    cm <- res$secondary[sample == "common"]
    rbind(if (nrow(cm)) cm else avail[model != "restricted"], avail[model == "restricted"])
  } else res$secondary

  figure_primary <- prim[, .(protein = assay, estimate, ci_low, ci_high, p, neg_log10_p, q, n, threshold, bonferroni, fdr)]
  figure_models <- if (nrow(fig_rows)) fig_rows[, .(model = model_label, protein = assay, beta = estimate, ci_low, ci_high, p, n, sample)] else
    data.table(model = character(), protein = character(), beta = numeric(), ci_low = numeric(), ci_high = numeric(),
               p = numeric(), n = integer())
  table_s2 <- merge(prim, res$spearman, by = "assay", all.x = TRUE)
  if (!is.null(res$unadjusted)) table_s2 <- merge(table_s2, res$unadjusted, by = "assay", all.x = TRUE)
  table_s2 <- merge(table_s2, res$qc[, .(assay, missing_fraction = round(missing_fraction, 3))], by = "assay", all.x = TRUE)
  qc_out <- res$qc[, .(assay, missing_fraction = round(missing_fraction, 3), analysed)]
  f <- copy(res$flow)
  f[, n_shown := fifelse(is.na(n_excluded), format(n_remaining, big.mark = ","),
                         vapply(n_excluded, small, "", min_cell = min_cell))]
  legacy_table <- if (is.null(res$legacy)) data.table() else res$legacy$table
  legacy_counts <- if (is.null(res$legacy)) data.table() else res$legacy$counts
  tables <- list(primary = figure_primary, models = figure_models, table_s2 = table_s2,
                 table_s4 = res$secondary, table_s5 = legacy_table, sensitivity = res$sensitivity,
                 qc = qc_out, table1 = res$table1, flow = f[, .(step, label, n_shown)])
  assert_aggregate(tables, nrow(res$qc), ids = res$ids)

  w(figure_primary, "primary.csv")
  w(figure_models, "models.csv")
  w(table_s2, "table_s2_primary_all_assays.csv")
  w(res$secondary, "table_s4_secondary_models.csv")
  w(legacy_table, "table_s5_legacy_staged.csv")
  w(legacy_counts, "table_s5_legacy_counts.csv")
  w(res$sensitivity, "table_s3_sensitivity.csv")
  w(qc_out, "figure_s1_protein_missingness.csv")
  w(res$table1, "table1_characteristics.csv")
  w(res$covariate_missing, "table_s1_covariate_missingness.csv")

  n_range <- function(x) if (length(x)) sprintf("n = %s–%s", format(min(x), big.mark = ","), format(max(x), big.mark = ",")) else "n = –"
  ext_n <- if (nrow(avail)) avail[model %in% unlist(plan$secondary$models), n] else integer()
  restr_n <- if (nrow(avail)) avail[model == "restricted", n] else integer()
  k <- res$keys
  flow <- list(
    primary_xlabel = sprintf("β (NPX per %s higher BD-PRS), basic model", res$beta_unit),
    figure3_note = "Basic and Extended 1-3: common sample per assay; Sensitivity: restricted sample",
    beta_unit = res$beta_unit,
    bonferroni_alpha = plan$primary$alpha, fdr_q = plan$primary$fdr_q,
    participants = c(
      lapply(seq_len(nrow(f)), function(i) list(f$label[i], paste("n =", f$n_shown[i]))),
      list(list("Basic model, per protein", n_range(prim$n)),
           list("Extended Models 1–3, per selected protein", n_range(ext_n)),
           list("Excluded: self-reported mental illness", paste("n =", small(k$cohort_n - k$restricted_n, min_cell))),
           list("Restricted (sensitivity) sample, per selected protein", n_range(restr_n)))),
    proteins = list(
      list("Assay variables in dataset", as.character(nrow(res$qc)), ""),
      list(sprintf("Excluded: >%g%% missing values", 100 * plan$proteins$max_missing_fraction),
           as.character(sum(!res$qc$analysed)), ""),
      list("Analysed in basic model [analysed]", as.character(nrow(prim)), ""),
      list(sprintf("P < %g/%d (primary) [bonferroni]", plan$primary$alpha, nrow(prim)), as.character(sum(prim$bonferroni)), ""),
      list(sprintf("FDR < %g%% (secondary) [fdr]", 100 * plan$primary$fdr_q), as.character(sum(prim$fdr)), ""))
  )
  jsonlite::write_json(flow, file.path(tmp, "flow.json"), auto_unbox = TRUE, pretty = TRUE)
  jsonlite::write_json(k, file.path(tmp, "result_keys.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
  pv <- provenance(plan, root, synthetic, config_paths, map, res$unresolved)
  writeLines(c(pv$line, sprintf("run_at: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
               paste0("input_sha256 ", names(res$inputs), ": ", res$inputs)), file.path(tmp, "PROVENANCE.txt"))
  pv$inputs_sha256 <- as.list(res$inputs)
  pv$rng <- RNGkind()
  writeLines(capture.output(sessionInfo()), file.path(tmp, "sessionInfo.txt"))
  pv$outputs_sha256 <- as.list(vapply(list.files(tmp, full.names = TRUE), function(p) digest::digest(file = p, algo = "sha256"), ""))
  names(pv$outputs_sha256) <- basename(names(pv$outputs_sha256))   # every file except provenance.json itself
  jsonlite::write_json(pv, file.path(tmp, "provenance.json"), auto_unbox = TRUE, pretty = TRUE)
  if (!file.rename(tmp, out)) stop("could not move the finished package to ", out, call. = FALSE)
  invisible(out)
}
