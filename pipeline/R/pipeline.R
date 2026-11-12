# The full analysis, from input files to the results list that write_package() exports.

check_plan <- function(plan) {
  fixed <- list(correction = "bonferroni", secondary = "BH",
                change_metric = "100 * |beta_model - beta_basic| / |beta_basic|")
  for (k in names(fixed)) if (!identical(plan$primary[[k]], fixed[[k]]))
    stop("primary$", k, " must be '", fixed[[k]], "' (the only option implemented)", call. = FALSE)
  if (!plan$exposure$supplied_unit %in% c("standard_deviation", "raw_unit"))
    stop("exposure$supplied_unit must be standard_deviation or raw_unit", call. = FALSE)
}

unresolved_decisions <- function(plan) {
  dec <- plan[["decisions"]]
  names(dec)[vapply(dec, function(x) !identical(x[["status"]], "resolved") || is.null(x[["value"]]), logical(1))]
}

run_analysis <- function(plan, map, root, synthetic = FALSE) {
  check_plan(plan)
  inp <- load_inputs(map, root, synthetic = synthetic)
  coh <- build_cohort(inp, plan, map)
  d <- coh$data
  qc_sample <- switch(plan$proteins$missing_computed_in, input = inp$data, dataset = coh$dataset, cohort = d,
                      stop("proteins$missing_computed_in must be input, dataset or cohort", call. = FALSE))
  qc <- protein_qc(qc_sample, inp$assays, plan$proteins$max_missing_fraction)
  assays <- qc[analysed == TRUE, assay]
  message(sprintf("dataset n = %d; cohort n = %d; assays analysed = %d of %d",
                  nrow(coh$dataset), nrow(d), length(assays), nrow(qc)))
  prim <- primary_analysis(d, assays, plan)
  selected <- prim[bonferroni == TRUE, assay]
  message(sprintf("primary: %d assays pass P < %.3g", length(selected), prim$threshold[1]))
  spear <- spearman_many(coh$dataset, assays, plan$exposure$variable)   # descriptive; pairwise PRS-assay complete
  # unadjusted OLS (STROBE 16a) and each assay's mean and SD in the cohort, for Table S2
  unadj <- fit_many(d, assays, plan$exposure$variable, character())[
    , .(assay, unadj_estimate = estimate, unadj_ci_low = ci_low, unadj_ci_high = ci_high, unadj_p = p)]
  unadj[, `:=`(assay_mean = vapply(assay, function(a) mean(d[[a]], na.rm = TRUE), 0),
               assay_sd = vapply(assay, function(a) sd(d[[a]], na.rm = TRUE), 0))]
  sec <- secondary_analysis(d, selected, plan)
  sens <- sensitivity_analysis(d, assays, selected, plan, labels = inp$labels)
  legacy <- if (isTRUE(plan$legacy_staged$enabled)) legacy_staged(coh$dataset, assays, plan) else NULL
  if (!is.null(legacy)) message("legacy staged counts: ", paste(legacy$counts$n_passed, collapse = " -> "))
  res <- list(primary = prim, spearman = spear, unadjusted = unadj, secondary = sec, sensitivity = sens, legacy = legacy, qc = qc,
              flow = coh$flow, inputs = inp$inputs, table1 = table1(d, plan$export$min_cell),
              covariate_missing = covariate_missingness(d, plan, plan$export$min_cell),
              ids = inp$data$eid, labels = inp$labels, unresolved = unresolved_decisions(plan),
              beta_unit = if (plan$exposure$scale == "standardise_in_cohort" ||
                              identical(plan$exposure$supplied_unit, "standard_deviation")) "1-SD" else "1-unit")
  res$keys <- result_keys(coh, qc, prim, sec, legacy, spear, plan)
  res
}
