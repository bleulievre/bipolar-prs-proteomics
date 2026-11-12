# Per-assay linear models. One assay = one outcome; the BD-PRS is the exposure.

# Backtick names that are not syntactic (Olink names such as "IL-6" would otherwise be read as arithmetic).
bt <- function(v) ifelse(make.names(v) == v, v, sprintf("`%s`", gsub("`", "", v, fixed = TRUE)))

check_model_data <- function(x, assay, exposure) {
  nums <- names(x)[vapply(x, is.numeric, logical(1))]
  if (length(nums) && any(!is.finite(as.matrix(x[, ..nums]))))
    stop("non-finite values in the model data for ", assay, call. = FALSE)
  bad <- names(x)[vapply(x, function(z) is.factor(z) && nlevels(droplevels(z)) < 2L, logical(1))]
  if (length(bad)) stop(assay, ": factor(s) with fewer than two observed levels: ", paste(bad, collapse = ", "), call. = FALSE)
  if (nrow(x) < 10) stop(assay, ": fewer than 10 complete observations", call. = FALSE)
  if (!(sd(x[[assay]]) > 0)) stop(assay, ": outcome has zero variance", call. = FALSE)
  if (!(sd(x[[exposure]]) > 0)) stop(assay, ": exposure has zero variance in this sample", call. = FALSE)
}

# Fit y ~ exposure + covariates [+ extra terms] by OLS on the rows given (complete cases are taken within).
# Returns one row: estimate, SE, 95% CI and P for the exposure (conventional OLS, or HC3 with the residual-df t
# distribution), -log10 P computed on the log scale (no underflow), n, residual df, partial R^2 of the exposure
# (nested models on the same rows) and the coefficient expressed in SDs of the assay (SD in the same rows).
fit_assay <- function(d, assay, exposure, covars, extra = character(), se = c("ols", "hc3"), rows = NULL,
                      allow_aliased = FALSE) {
  se <- match.arg(se)
  extra_vars <- if (length(extra)) all.vars(reformulate(extra)) else character()
  vars <- unique(c(assay, exposure, covars, extra_vars))
  x <- if (is.null(rows)) d[, ..vars] else d[rows, ..vars]
  x <- x[complete.cases(x)]
  check_model_data(x, assay, exposure)
  full <- lm(reformulate(c(bt(exposure), bt(covars), extra), response = bt(assay)), data = x)
  aliased <- names(coef(full))[is.na(coef(full))]
  if (length(aliased) && !allow_aliased)
    stop(assay, ": rank-deficient model (terms dropped: ", paste(aliased, collapse = ", "), ")", call. = FALSE)
  if (is.na(coef(full)[2])) stop(assay, ": the exposure is aliased", call. = FALSE)
  reduced <- lm(reformulate(c("1", bt(covars), extra), response = bt(assay)), data = x)
  rss_f <- sum(residuals(full)^2)
  rss_r <- sum(residuals(reduced)^2)
  term <- names(coef(full))[2]   # the exposure is the first term after the intercept
  b <- unname(coef(full)[term])
  v <- if (se == "ols") vcov(full, complete = FALSE) else sandwich::vcovHC(full, type = "HC3")
  s <- sqrt(v[term, term])
  df <- full$df.residual
  tq <- qt(0.975, df)
  tstat <- b / s
  log10p <- (pt(abs(tstat), df, lower.tail = FALSE, log.p = TRUE) + log(2)) / log(10)
  sd_y <- sd(x[[assay]])
  out <- data.table(assay = assay, n = nrow(x), df_resid = df, estimate = b, se = s,
                    ci_low = b - tq * s, ci_high = b + tq * s, p = 10^log10p, neg_log10_p = -log10p,
                    partial_r2 = (rss_r - rss_f) / rss_r,
                    # outcome SDs per unit of the exposure (per SD only if the exposure is standardised)
                    estimate_outcome_sd = b / sd_y, ci_low_outcome_sd = (b - tq * s) / sd_y, ci_high_outcome_sd = (b + tq * s) / sd_y,
                    aliased = paste(aliased, collapse = ";"))
  if (any(!is.finite(unlist(out[, .(estimate, se, ci_low, ci_high, neg_log10_p, partial_r2)]))) ||
      out$partial_r2 < -1e-12 || out$partial_r2 > 1)
    stop(assay, ": invalid model statistics", call. = FALSE)
  out
}

fit_many <- function(d, assays, exposure, covars, ...) {
  rbindlist(lapply(assays, function(a) fit_assay(d, a, exposure, covars, ...)))
}

# Spearman rank correlation, pairwise-complete (as in the original code); descriptive in the primary design.
spearman_many <- function(d, assays, exposure) {
  rbindlist(lapply(assays, function(a) {
    ok <- !is.na(d[[a]]) & !is.na(d[[exposure]])
    x <- d[[exposure]][ok]; y <- d[[a]][ok]
    if (sum(ok) < 10 || !(sd(x) > 0) || !(sd(y) > 0)) stop(a, ": Spearman correlation not estimable", call. = FALSE)
    ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
    data.table(assay = a, n_spearman = sum(ok), rho = unname(ct$estimate), p_spearman = ct$p.value)
  }))
}

# Primary analysis: basic model on every analysed assay, one Bonferroni threshold and one BH adjustment.
# Thresholds use strict "<", as in the original code.
primary_analysis <- function(d, assays, plan) {
  if (!length(assays)) stop("no assays to analyse", call. = FALSE)
  res <- fit_many(d, assays, plan$exposure$variable, unlist(plan$covariates$core), se = plan$primary$se)
  m <- length(assays)
  res[, `:=`(q = p.adjust(p, "BH"), threshold = plan$primary$alpha / m)]
  res[, `:=`(bonferroni = p < threshold, fdr = q < plan$primary$fdr_q)]
  res[]
}
