# Unit tests: recoding, input checks, cohort building, protein QC, model statistics, staged procedure, exports.

test_that("codes recoding follows the documented precedence", {
  dt <- data.table(a0 = c(2, 0, -3, NA, 5, 2), a1 = c(0, NA, NA, NA, NA, -3))
  spec <- list(type = "codes", columns_prefix = "a", codes_1 = list(2), codes_0 = list(0, 1), codes_na = list(-3))
  expect_equal(recode_codes(dt, spec, "x", "."), c(1, 0, NA, NA, NA, 1))
  spec$other <- 0
  spec$none <- 0
  expect_equal(recode_codes(dt, spec, "x", "."), c(1, 0, NA, 0, 0, 1))
})

test_that("a code list with a non-numeric line stops the run instead of flagging everyone", {
  f <- tempfile()
  writeLines(c("# header", "1286", "depression"), f)
  expect_error(read_codes(f), "non-numeric")
  writeLines("# nothing yet", f)
  expect_error(read_codes(f), "empty")
  writeLines(c("﻿1286", "1291"), f)
  expect_equal(read_codes(f), c(1286, 1291))
  # an empty cell must never match a code
  dt <- data.table(a0 = c(NA, NA, 1286), a1 = c(NA, 5, NA))
  writeLines(c("1286"), f)
  expect_equal(recode_codes(dt, list(columns_prefix = "a", codes_1_file = basename(f), other = 0, none = 0), "x", dirname(f)),
               c(0, 0, 1))
})

test_that("YAML key names cannot turn into booleans or be partially matched", {
  for (m in c("data_map.yaml", "data_map_glasgow.yaml")) {
    map <- yaml::read_yaml(file.path(ROOT, "config", m))
    for (v in map$variables) expect_false(any(c("TRUE", "FALSE") %in% names(v)))
  }
  expect_equal(spec_columns(data.table(p1_a0 = 1), list(columns_prefix = "p1_a"), "x"), "p1_a0")
})

test_that("factor recodes match both sex codings and send unknown codes to NA", {
  dt <- data.table(s = c(0, 1, 9), sex_plink = c(1, 2, 2), b = c(-3, 5, NA))
  f <- recode_variable(dt, list(type = "factor", column = "s", levels = list(`0` = "female", `1` = "male")), "sex", ".")
  expect_equal(as.character(f), c("female", "male", NA))
  g <- yaml::read_yaml(file.path(ROOT, "config", "data_map_glasgow.yaml"))$variables$sex
  expect_equal(as.character(recode_variable(dt, g, "sex", ".")), c("male", "female", "female"))   # sex_plink 1/2
  a <- recode_variable(dt, list(type = "split", column = "b", cut = 0, below = "BiLEVE", at_or_above = "Axiom"), "array", ".")
  expect_equal(as.character(a), c("BiLEVE", "Axiom", NA))
})

test_that("field 6150 recodes: stroke and IHD separately, high blood pressure alone is not either", {
  dt <- data.table(p6150_i0_a0 = c(1, 4, -7, NA, 3, -3), p6150_i0_a1 = c(NA, NA, NA, NA, 2, NA))
  map <- yaml::read_yaml(file.path(ROOT, "config", "data_map.yaml"))
  expect_equal(recode_variable(dt, map$variables$ihd, "ihd", "."), c(1, 0, 0, NA, 1, NA))
  expect_equal(recode_variable(dt, map$variables$stroke, "stroke", "."), c(0, 0, 0, NA, 1, NA))
})

test_that("waist-to-hip ratio rejects non-positive measurements", {
  dt <- data.table(w = c(90, 0, 80), h = c(100, 100, 0))
  expect_equal(recode_variable(dt, list(type = "ratio", numerator = "w", denominator = "h"), "whr", "."), c(0.9, NA, NA))
})

test_that("duplicate or missing participant ids stop the run", {
  expect_error(check_ids(data.table(eid = c(1, 1, 2)), "participants", "eid"), "duplicate")
  expect_error(check_ids(data.table(eid = c(1, NA)), "prs", "eid"), "missing")
})

test_that("an optional variable is skipped only when its column is absent", {
  dt <- data.table(x = c("a", "b"))
  expect_error(recode_variable(dt, list(type = "numeric", column = "x"), "v", "."), "non-numeric")
  expect_s3_class(tryCatch(recode_variable(dt, list(type = "numeric", column = "y"), "v", "."), error = function(e) e),
                  "missing_column")
})

test_that("relatedness removal leaves no related pair and does not depend on row order", {
  pairs <- data.table(ID1 = c("1", "2", "3", "7"), ID2 = c("2", "3", "4", "8"))
  ids <- as.character(1:10)
  r1 <- drop_related(ids, pairs, seed = 1)
  r2 <- drop_related(ids, pairs[c(4, 2, 1, 3)][, .(ID1 = ID2, ID2 = ID1)], seed = 1)
  expect_setequal(r1, r2)
  kept <- setdiff(ids, r1)
  expect_false(any(pairs$ID1 %in% kept & pairs$ID2 %in% kept))
})

test_that("protein QC keeps an assay with exactly the threshold fraction missing", {
  d <- data.table(x = c(rep(NA, 10), 1:90), y = c(rep(NA, 11), 1:89))
  expect_equal(protein_qc(d, c("x", "y"), 0.10)$analysed, c(TRUE, FALSE))
})

test_that("partial R2 from nested models equals t^2 / (t^2 + df); CI matches confint(); HC3 matches sandwich", {
  set.seed(1)
  d <- data.table(prs = rnorm(500), age = runif(500, 40, 70), sex = factor(sample(c("f", "m"), 500, TRUE)),
                  bmi = rnorm(500, 27, 4))
  d[, y := 0.2 * prs + 0.01 * age + rnorm(500)]
  d[sample(500, 20), y := NA]
  d[sample(500, 15), bmi := NA]
  r <- fit_assay(d, "y", "prs", c("age", "sex", "bmi"))
  x <- d[complete.cases(d)]
  fit <- lm(y ~ prs + age + sex + bmi, data = x)
  t <- coef(summary(fit))["prs", "t value"]
  expect_equal(r$partial_r2, t^2 / (t^2 + fit$df.residual), tolerance = 1e-10)
  expect_equal(c(r$ci_low, r$ci_high), unname(confint(fit)["prs", ]), tolerance = 1e-10)
  expect_equal(r$n, nrow(x))
  expect_equal(r$estimate_outcome_sd, r$estimate / sd(x$y), tolerance = 1e-10)   # SD in the regression sample
  expect_equal(r$p, coef(summary(fit))["prs", "Pr(>|t|)"], tolerance = 1e-10)
  h <- fit_assay(d, "y", "prs", c("age", "sex", "bmi"), se = "hc3")
  expect_equal(h$se, sqrt(sandwich::vcovHC(fit, type = "HC3")["prs", "prs"]), tolerance = 1e-10)
})

test_that("extreme associations keep a finite -log10 P", {
  set.seed(9)
  d <- data.table(prs = rnorm(5000), age = rnorm(5000))
  d[, y := 3 * prs + rnorm(5000, 0, 0.1)]
  r <- fit_assay(d, "y", "prs", "age")
  expect_true(is.finite(r$neg_log10_p) && r$neg_log10_p > 300)
})

test_that("bad model inputs stop with a clear message", {
  set.seed(2)
  d <- data.table(prs = rnorm(100), age = rnorm(100), const = 1, f1 = factor(rep("a", 100)), y = rnorm(100))
  d[, age2 := 2 * age]
  expect_error(fit_assay(d, "const", "prs", "age"), "zero variance")
  expect_error(fit_assay(d, "y", "prs", c("age", "f1")), "fewer than two observed levels")
  expect_error(fit_assay(d, "y", "prs", c("age", "age2")), "rank-deficient")
  expect_silent(fit_assay(d, "y", "prs", c("age", "age2"), allow_aliased = TRUE))
  d[1, y := Inf]
  expect_error(fit_assay(d, "y", "prs", "age"), "non-finite")
})

test_that("non-syntactic assay names such as IL-6 are modelled correctly", {
  set.seed(3)
  d <- data.table(prs = rnorm(300), age = rnorm(300))
  d[, `IL-6` := 0.5 * prs + rnorm(300)]
  r <- fit_assay(d, "IL-6", "prs", "age")
  expect_equal(r$estimate, unname(coef(lm(d$`IL-6` ~ d$prs + d$age))[2]), tolerance = 1e-10)
})

test_that("the primary analysis uses one Bonferroni threshold and BH across all assays", {
  set.seed(2)
  d <- data.table(prs = rnorm(400), age = rnorm(400))
  for (j in 1:5) d[, (paste0("a", j)) := (j == 1) * 0.5 * prs + rnorm(400)]
  plan <- list(exposure = list(variable = "prs"), covariates = list(core = "age"),
               primary = list(se = "ols", alpha = 0.05, fdr_q = 0.05))
  r <- primary_analysis(d, paste0("a", 1:5), plan)
  expect_equal(r$threshold, rep(0.01, 5))
  expect_equal(r$q, p.adjust(r$p, "BH"))
  expect_true(r[assay == "a1", bonferroni])
})

test_that("common-sample fits use the same rows in every model; change is NA when the basic estimate is ~0", {
  set.seed(4)
  d <- data.table(prs = rnorm(500), age = rnorm(500), bmi = rnorm(500), flag = rbinom(500, 1, 0.1))
  d[, a1 := 0.4 * prs + rnorm(500)]
  d[sample(500, 40), bmi := NA]
  plan <- list(exposure = list(variable = "prs"), covariates = list(core = "age", extended = list(ext1 = "bmi")),
               secondary = list(models = list("ext1"), common_sample = TRUE,
                                restricted = list(exclude_flag = "flag", na_action = "exclude", model = "core")))
  s <- secondary_analysis(d, "a1", plan)
  expect_equal(length(unique(s[sample == "common", n])), 1)
  expect_true(is.na(pct_change(1, 1e-12)))
})

test_that("restricted rows: unknown history is excluded or kept as configured", {
  d <- data.table(f = c(0, 1, NA))
  expect_equal(restricted_rows(d, "f", "exclude"), c(TRUE, FALSE, FALSE))
  expect_equal(restricted_rows(d, "f", "keep"), c(TRUE, FALSE, TRUE))
})

test_that("the staged procedure uses the fixed screen threshold, then alpha / number carried", {
  set.seed(3)
  d <- data.table(prs_bd = rnorm(3000), age = rnorm(3000), bmi = rnorm(3000), flag = rbinom(3000, 1, 0.1))
  d[sample(3000, 100), flag := NA]
  for (j in 1:6) d[, (paste0("a", j)) := (j <= 2) * 0.2 * prs_bd + rnorm(3000)]
  plan <- list(exposure = list(variable = "prs_bd"),
               covariates = list(core = "age", extended = list(ext1 = "bmi")),
               secondary = list(restricted = list(exclude_flag = "flag")),
               legacy_staged = list(alpha = 0.05, screen_threshold = 1e-4, stages = list("core", "ext1", "restricted"),
                                    restricted_model = "ext1", restricted_na_action = "keep"))
  lg <- legacy_staged(d, paste0("a", 1:6), plan)
  expect_equal(lg$counts[stage == "spearman", threshold], 1e-4)
  carried <- lg$counts[stage == "spearman", n_passed]
  expect_equal(lg$counts[stage == "core", threshold], 0.05 / carried)
  expect_equal(lg$counts[stage == "ext1", threshold], 0.05 / lg$counts[stage == "core", n_passed])
  # restricted stage keeps participants with unknown history (as UKB.R did)
  expect_equal(lg$table[stage == "restricted", n][1], sum(d$flag %in% 0 | is.na(d$flag)))
})

test_that("small counts and their percentages are suppressed", {
  d <- data.table(x = c(rep(1, 3), rep(0, 997)), y = c(rep(1, 995), rep(0, 5)), z = c(rep(1, 500), rep(0, 500)))
  TABLE1_VARS[c("x", "y", "z")] <<- c("x", "y", "z")
  t1 <- table1(d, 10)
  expect_match(t1[variable == "x", summary], "suppressed")
  expect_match(t1[variable == "y", summary], "suppressed")   # complement of 5
  expect_match(t1[variable == "z", summary], "500 \\(50.0%\\)")
  expect_false(grepl("%", t1[variable == "x", summary]))
})

test_that("exports refuse participant-level tables", {
  expect_error(assert_aggregate(list(bad = data.table(eid = 1:3, x = 1)), 10), "participant id")
  expect_error(assert_aggregate(list(bad = data.table(x = 1:1000)), 10), "more rows")
  expect_error(assert_aggregate(list(bad = data.table(v = c("1000042", "x"))), 10, ids = "1000042"), "contains participant ids")
  expect_true(assert_aggregate(list(ok = data.table(assay = "a", p = 0.1)), 10))
})

test_that("provenance outside a git repository is never FROZEN", {
  plan <- list(plan_id = "x", version = "1", status = "FROZEN", frozen_date = "2026-11-01", frozen_by = list("A"))
  pv <- provenance(plan, tempdir(), synthetic = FALSE, config_paths = character())
  expect_equal(pv$kind, "DRAFT")
  expect_match(pv$line, "^DRAFT ")
})

test_that("a shifted assay column range is caught by the Olink name check", {
  prot <- data.table(V1 = 1:3, n_eid = c("a", "b", "c"), age_0_0 = c(50, 60, 70), IL6 = c(1, 2, 3), LEP = c(2, 3, 4))
  map <- list(proteins = list(columns = list(from = 4, to = 5, check_names = "config/olink_assay_ukb.tsv")))
  allow <- TRUE   # first/last not set here; this test is about check_names
  expect_equal(select_assays(prot, map, "n_eid", ROOT, allow), c("IL6", "LEP"))
  map$proteins$columns$from <- 3
  map$proteins$columns$to <- 4
  expect_error(select_assays(prot, map, "n_eid", ROOT, allow), "not Olink assay names")
})

test_that("the Olink range guard stops a run whose first/last assay names do not match", {
  prot <- data.table(n_eid = c("a", "b"), x = c(1, 2), IL6 = c(1, 2), LEP = c(2, 3))
  map <- list(proteins = list(columns = list(from = 3, to = 4, first = "IL6", last = "LEP")))
  expect_equal(select_assays(prot, map, "n_eid", ROOT, FALSE), c("IL6", "LEP"))
  map$proteins$columns$first <- "X"
  expect_error(select_assays(prot, map, "n_eid", ROOT, FALSE), "shifted")
  map$proteins$columns$first <- NULL
  expect_error(select_assays(prot, map, "n_eid", ROOT, FALSE), "first/last")
})

test_that("a selected kidney-sensitive assay without eGFR stops the run", {
  set.seed(5)
  d <- data.table(prs = rnorm(300), age = rnorm(300))
  d[, a1 := 0.5 * prs + rnorm(300)]
  plan <- list(exposure = list(variable = "prs"), covariates = list(core = "age"),
               sensitivity = list(hc3 = FALSE, age_squared = FALSE, egfr_required_for = list("PROTK"),
                                  optional = list(egfr = "egfr")))
  expect_error(sensitivity_analysis(d, "a1", "a1", plan, labels = c(a1 = "PROTK")), "eGFR")
  expect_silent(sensitivity_analysis(d, "a1", "a1", plan, labels = c(a1 = "PROTB")))
})

test_that("unresolved decisions keep a run from being FROZEN", {
  plan <- yaml::read_yaml(file.path(ROOT, "config", "analysis_plan.yaml"))
  expect_true(length(unresolved_decisions(plan)) > 0)
  plan$decisions <- lapply(plan$decisions, function(x) list(status = "resolved", value = "x"))
  expect_length(unresolved_decisions(plan), 0)
})

test_that("an unadjusted model (no covariates) matches lm(y ~ prs)", {
  set.seed(8)
  d <- data.table(prs = rnorm(200))
  d[, y := 0.3 * prs + rnorm(200)]
  r <- fit_assay(d, "y", "prs", character())
  expect_equal(r$estimate, unname(coef(lm(y ~ prs, data = d))[2]), tolerance = 1e-12)
})
