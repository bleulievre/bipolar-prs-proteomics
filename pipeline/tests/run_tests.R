#!/usr/bin/env Rscript
# Rscript tests/run_tests.R   (from the repository root). The pipeline tests generate fresh simulated data in a
# temporary directory, so stale files in synthetic/ cannot affect them.
options(warnPartialMatchDollar = TRUE)
suppressPackageStartupMessages({ library(testthat); library(data.table) })
ROOT <- normalizePath(".")
for (f in c("io.R", "cohort.R", "models.R", "secondary.R", "legacy.R", "export.R", "pipeline.R"))
  source(file.path(ROOT, "R", f))
res <- test_dir(file.path(ROOT, "tests", "testthat"), reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("\n%d tests, %d failed, %d errors, %d skipped\n", nrow(df), sum(df$failed > 0), sum(df$error), sum(df$skipped)))
if (any(df$failed > 0) || any(df$error)) quit(status = 1)
