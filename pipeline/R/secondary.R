# Secondary (descriptive) analyses for the assays selected by the primary analysis, and sensitivity analyses.
# None of these select assays; they describe how the estimates change.

MODEL_LABELS <- c(core = "Basic", ext1 = "Extended 1", ext2 = "Extended 2", ext3 = "Extended 3",
                  restricted = "Sensitivity")

model_covariates <- function(plan) {
  prev <- unlist(plan$covariates$core)
  out <- list(core = prev)
  for (m in names(plan$covariates$extended)) {
    prev <- c(prev, unlist(plan$covariates$extended[[m]]))
    out[[m]] <- prev
  }
  out
}

# Rows in the restricted sample. na_action "exclude": unknown history is excluded (conservative);
# "keep": only participants flagged 1 are excluded (as in the original code).
restricted_rows <- function(d, flag, na_action) {
  x <- d[[flag]]
  if (is.null(x)) stop("restricted-sample flag '", flag, "' is not a mapped variable", call. = FALSE)
  switch(na_action, exclude = !is.na(x) & x == 0, keep = is.na(x) | x != 1,
         stop("restricted na_action must be exclude or keep", call. = FALSE))
}

pct_change <- function(est, base) {
  fifelse(is.na(base) | abs(base) <= sqrt(.Machine$double.eps), NA_real_, 100 * abs(est - base) / abs(base))
}

# Long table: model x assay x sample. "available" = complete cases for each model's own variables;
# "common" = rows complete for every covariate in every model and the assay (same rows for all models, so the
# change between models reflects the added covariates only).
secondary_analysis <- function(d, selected, plan) {
  if (!length(selected)) return(data.table())
  covs <- model_covariates(plan)
  exposure <- plan$exposure$variable
  models <- c("core", unlist(plan$secondary$models))
  if (length(setdiff(models, names(covs)))) stop("unknown model(s) in secondary$models", call. = FALSE)
  all_covs <- unique(unlist(covs[models]))
  rs <- plan$secondary$restricted
  out <- list()
  for (a in selected) {
    for (m in models)
      out[[length(out) + 1]] <- fit_assay(d, a, exposure, covs[[m]])[, `:=`(model = m, sample = "available")]
    if (isTRUE(plan$secondary$common_sample)) {
      common <- which(complete.cases(d[, c(a, exposure, all_covs), with = FALSE]))
      for (m in models)
        out[[length(out) + 1]] <- fit_assay(d, a, exposure, covs[[m]], rows = common)[, `:=`(model = m, sample = "common")]
    }
    rr <- which(restricted_rows(d, rs$exclude_flag, rs$na_action))
    out[[length(out) + 1]] <- fit_assay(d, a, exposure, covs[[rs$model]], rows = rr)[
      , `:=`(model = "restricted", sample = "available")]
  }
  res <- rbindlist(out)
  # the restricted-sample row is compared with the available-case basic model ("restriction-associated change")
  base <- res[model == "core", .(assay, sample, beta_basic = estimate)]
  res <- merge(res, base, by = c("assay", "sample"), all.x = TRUE)
  res[, `:=`(abs_pct_change = pct_change(estimate, beta_basic),
             same_sign = sign(estimate) == sign(beta_basic),
             model_label = MODEL_LABELS[model])]
  setcolorder(res, c("model", "model_label", "sample", "assay"))
  res[order(assay, sample, match(model, names(MODEL_LABELS)))]
}

# Sensitivity analyses. HC3 and centred age^2 for every analysed assay; optional analyses (eGFR, schizophrenia
# PRS, PRS x sex) for the selected assays, run only when the variable is available. If a selected assay is listed
# under sensitivity$egfr_required_for (assays whose plasma level depends strongly on kidney function), the eGFR
# analysis is required.

sensitivity_analysis <- function(d, assays, selected, plan, labels = NULL) {
  exposure <- plan$exposure$variable
  core <- unlist(plan$covariates$core)
  sens <- plan$sensitivity
  out <- list()
  if (isTRUE(sens$hc3))
    out$hc3 <- fit_many(d, assays, exposure, core, se = "hc3")[, analysis := "HC3 standard errors (t, residual df)"]
  if (isTRUE(sens$age_squared)) {
    age_c <- mean(d$age, na.rm = TRUE)
    out$age2 <- fit_many(d, assays, exposure, core, extra = sprintf("I((age - %.6f)^2)", age_c))[
      , analysis := sprintf("Quadratic age term (age centred at %.2f)", age_c)]
  }
  opt <- sens$optional
  if (length(selected)) {
    targets <- toupper(if (is.null(labels)) selected else unname(labels[selected]))
    kidney <- toupper(unlist(sens$egfr_required_for))
    if (any(targets %in% kidney) && !isTRUE(opt$egfr %in% names(d)))
      stop("a kidney-sensitive assay was selected, so the eGFR sensitivity analysis is required, but '",
           opt$egfr, "' is not available", call. = FALSE)
    if (!is.null(opt$egfr) && opt$egfr %in% names(d))
      out$egfr <- fit_many(d, selected, exposure, c(core, opt$egfr))[, analysis := "Adjusted for eGFR"]
    if (!is.null(opt$scz_prs) && opt$scz_prs %in% names(d))
      out$scz <- fit_many(d, selected, exposure, c(core, opt$scz_prs))[, analysis := "Adjusted for schizophrenia PRS"]
    if (isTRUE(opt$prs_by_sex))
      out$sex <- rbindlist(lapply(selected, function(a) interaction_test(d, a, exposure, core, "sex")))[
        , analysis := "BD-PRS x sex interaction (OLS)"]
  }
  rbindlist(out, fill = TRUE)
}

# Wald test of the exposure x modifier interaction (one row per assay; estimate = interaction coefficient).
interaction_test <- function(d, assay, exposure, covars, modifier) {
  x <- d[, unique(c(assay, exposure, covars, modifier)), with = FALSE]
  x <- x[complete.cases(x)]
  check_model_data(x, assay, exposure)
  term <- paste0(exposure, ":", modifier)
  fit <- lm(reformulate(c(bt(exposure), bt(covars), term), response = bt(assay)), data = x)
  cf <- summary(fit)$coefficients
  row <- which(startsWith(rownames(cf), term))
  if (length(row) != 1) stop("interaction term not identified for ", assay, call. = FALSE)
  data.table(assay = assay, n = nrow(x), estimate = cf[row, 1], se = cf[row, 2], p = cf[row, 4],
             term = rownames(cf)[row])
}
