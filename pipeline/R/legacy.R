# The dissertation's staged procedure, reproduced from the original code (UKB.R, 2025) for Supplementary Table S5.
#   0. sample: the supplied dataset (no requirement for complete core covariates); missingness as in the plan
#   1. Spearman screen of each assay against the BD-PRS, pairwise complete, kept if P < screen_threshold
#      (1e-4 in UKB.R line 302, although the dissertation text says 0.05/365)
#   2-5. basic model, then Extended Models 1-3: each refits the assays carried forward on complete cases for its
#      own variables and keeps those with P < alpha / (number carried into that step)
#   6. restricted sample with the Extended Model 3 formula (UKB.R line 733), same threshold rule. The sample was built
#      by merge(data, supplement) (UKB.R line 697, a default inner join, so participants absent from the supplement
#      were dropped) and then anyMentalIllness != 1 | is.na(anyMentalIllness) (line 705). NOTE: the script writes
#      "sensitve.csv" but reads "sensitive.csv", and the merge step is commented out, so this definition is the most
#      likely one, not a verified one (cross-validation review, 2026-10-20).
# The primary analysis replaced this procedure because it does not control the error rate across assays.

legacy_staged <- function(d, assays, plan) {
  lg <- plan$legacy_staged
  d <- copy(d)
  for (v in names(lg$recode_override)) {            # e.g. smoking/alcohol coded as in UKB.R
    src <- lg$recode_override[[v]]
    if (src %in% names(d)) d[, (v) := get(src)]
  }
  lit <- lg$literal_thresholds
  thr_for <- function(stage, n_in) if (!is.null(lit[[stage]])) lit[[stage]] else alpha / n_in
  pass <- function(p, thr, screen = FALSE) {
    ok <- !is.na(p) & p < thr
    if (isTRUE(lg$drop_zero_p)) ok <- ok & p != 0 & !(screen & p == 1)
    ok
  }
  alpha <- lg$alpha
  exposure <- plan$exposure$variable
  covs <- model_covariates(list(covariates = list(core = plan$covariates$core,
                                                   extended = if (is.null(lg$extended)) plan$covariates$extended else lg$extended)))
  al <- isTRUE(lg$allow_aliased)
  screen_on <- if ("prs_bd_raw" %in% names(d)) "prs_bd_raw" else exposure   # UKB.R screened on the raw score
  sp <- spearman_many(d, assays, screen_on)
  thr <- if (!is.null(lit[["spearman"]])) lit[["spearman"]] else lg$screen_threshold
  ok <- pass(sp$p_spearman, thr, screen = TRUE)
  rows <- list(data.table(stage = "spearman", assay = sp$assay, estimate = sp$rho, p = sp$p_spearman,
                          n = sp$n_spearman, threshold = thr, threshold_rule = alpha / length(assays),
                          passed = ok, aliased = ""))
  carried <- sp$assay[ok]
  for (s in unlist(lg$stages)) {
    if (!length(carried)) break
    thr <- thr_for(s, length(carried))
    r <- if (s == "restricted") {
      keep <- restricted_rows(d, plan$secondary$restricted$exclude_flag, lg$restricted_na_action)
      req <- lg$restricted_require   # UKB.R merge() = inner join; only for data maps that define this variable
      if (!is.null(req) && req %in% names(d)) keep <- keep & d[[req]] %in% 1
      fit_many(d[keep], carried, exposure, covs[[lg$restricted_model]], allow_aliased = al)
    } else {
      fit_many(d, carried, exposure, covs[[s]], allow_aliased = al)
    }
    r <- r[, .(stage = s, assay, estimate, p, n, threshold = thr, threshold_rule = alpha / length(carried),
               passed = pass(p, thr), aliased)]
    rows[[length(rows) + 1]] <- r
    carried <- r[passed == TRUE, assay]
  }
  out <- rbindlist(rows)
  counts <- out[, .(n_tested = .N, n_passed = sum(passed), threshold = threshold[1]), by = stage]
  list(table = out, counts = counts, spearman = sp)
}
