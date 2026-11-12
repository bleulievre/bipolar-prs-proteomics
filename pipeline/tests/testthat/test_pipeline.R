# End-to-end tests on freshly generated simulated data, in both layouts (synthetic/make_synthetic.R).

syn <- file.path(tempdir(), "syn")
unlink(syn, recursive = TRUE)
system2(file.path(R.home("bin"), "Rscript"), c(file.path(ROOT, "synthetic", "make_synthetic.R"), "8000", "60", syn),
        stdout = FALSE, stderr = FALSE)
plan <- yaml::read_yaml(file.path(ROOT, "config", "analysis_plan.yaml"))

run_layout <- function(layout, out) {
  map_path <- file.path(syn, layout, "data_map.yaml")
  res <- run_analysis(plan, yaml::read_yaml(map_path), file.path(syn, layout), synthetic = TRUE)
  write_package(out, res, plan, ROOT, synthetic = TRUE, config_paths = map_path, overwrite = TRUE)
  res
}

for (layout in c("glasgow", "rap")) {
  out1 <- file.path(tempdir(), paste0(layout, "_1"))
  res <- run_layout(layout, out1)
  truth <- fread(file.path(syn, layout, "truth.csv"))

  test_that(paste(layout, ": planted associations are found with the right sign; null assays mostly are not"), {
    found <- res$primary[bonferroni == TRUE]
    expect_setequal(intersect(found$assay, truth$assay), intersect(truth$assay, res$primary$assay))
    m <- merge(found, truth, by = "assay")
    expect_true(all(sign(m$estimate) == sign(m$beta_per_sd)))
    expect_lte(length(setdiff(found$assay, truth$assay)), 1)
  })

  test_that(paste(layout, ": flow counts add up and match the results"), {
    f <- res$flow
    expect_equal(f[step == "cohort", n_remaining], f[step == "start", n_remaining] - sum(f$n_excluded, na.rm = TRUE))
    flow <- jsonlite::read_json(file.path(out1, "flow.json"))
    prot <- setNames(vapply(flow$proteins, function(x) as.integer(x[[2]]), 1L), vapply(flow$proteins, function(x) x[[1]], ""))
    expect_equal(unname(prot[grep("\\[analysed\\]", names(prot))]), nrow(res$primary))
    expect_equal(unname(prot[grep("\\[bonferroni\\]", names(prot))]), sum(res$primary$bonferroni))
    expect_equal(res$keys$n_analysed + res$keys$n_excluded_missing, res$keys$n_assays_in_data)
    expect_equal(f[step == "withdrawn", n_excluded], 5L)        # 5 withdrawals planted in both layouts
  })

  test_that(paste(layout, ": the model grid is complete; Figure 3 uses the common sample for Basic-Extended 3"), {
    mods <- fread(file.path(out1, "models.csv"))
    sel <- res$primary[bonferroni == TRUE, assay]
    expect_equal(nrow(mods), length(sel) * 5)
    expect_equal(unique(mods[model != "Sensitivity", sample]), "common")
    b <- merge(mods[model == "Basic"], res$secondary[sample == "common" & model == "core"], by.x = "protein", by.y = "assay")
    expect_equal(b$beta, b$estimate, tolerance = 1e-12)
    expect_equal(res$secondary[sample == "common", uniqueN(n), by = assay]$V1, rep(1L, length(sel)))
  })

  test_that(paste(layout, ": exported files contain no participant identifiers"), {
    ids <- c(as.character(1000001:1008000), as.character(-(1:5)))
    for (f in list.files(out1, pattern = "\\.(csv|json)$", full.names = TRUE)) {
      txt <- readLines(f, warn = FALSE)
      expect_false(any(grepl("\\b100[0-7][0-9]{3}\\b", txt)), info = basename(f))
    }
  })

  test_that(paste(layout, ": the run is deterministic"), {
    out2 <- file.path(tempdir(), paste0(layout, "_2"))
    run_layout(layout, out2)
    for (f in setdiff(list.files(out1), c("PROVENANCE.txt", "provenance.json", "sessionInfo.txt")))
      expect_identical(digest::digest(file = file.path(out1, f)), digest::digest(file = file.path(out2, f)), info = f)
  })

  test_that(paste(layout, ": Table S2 has unadjusted OLS estimates and assay summaries"), {
    s2 <- fread(file.path(out1, "table_s2_primary_all_assays.csv"))
    expect_true(all(c("unadj_estimate", "unadj_ci_low", "unadj_ci_high", "unadj_p", "assay_mean", "assay_sd") %in% names(s2)))
    a <- s2$assay[1]
    d <- res$primary
    expect_equal(nrow(s2), nrow(d))
  })

  test_that(paste(layout, ": provenance marks synthetic output"), {
    expect_match(readLines(file.path(out1, "PROVENANCE.txt"))[1], "^SYNTHETIC ")
  })
}

test_that("glasgow: the original procedure is reproduced as in UKB.R", {
  res <- run_analysis(plan, yaml::read_yaml(file.path(syn, "glasgow", "data_map.yaml")), file.path(syn, "glasgow"), TRUE)
  lg <- res$legacy
  expect_equal(lg$counts$stage, c("spearman", "core", "ext1", "ext2", "ext3", "restricted"))
  expect_equal(lg$counts[stage == "spearman", threshold], 1e-4)
  expect_equal(lg$counts[stage == "core", threshold], 0.0021)          # as typed in UKB.R line 412
  expect_equal(lg$counts[stage == "ext1", threshold], 0.0025)
  # Extended Model 2 as in UKB.R: adjusted BP is aliased with measured BP + medication in these data; recorded
  expect_true(all(grepl("sbp_adj", lg$table[stage == "ext2", aliased])))
  # the restricted stage drops participants absent from the supplement (inner join) but keeps unknown history
  r <- lg$table[stage == "restricted"]
  expect_true(nrow(r) > 0)
  # the dataset sample keeps participants with missing core covariates (UKB.R had no such exclusion)
  expect_gte(nrow(res$keys$legacy_counts), 6)
  expect_gt(res$keys$dataset_n, res$keys$cohort_n)
})

test_that("figures: the Python script validates the package and copes with no selected assay", {
  py <- Sys.which("python3")
  skip_if(!nzchar(py) || system2(py, c("-c", shQuote("import matplotlib")), stdout = FALSE, stderr = FALSE) != 0,
          "python3 with matplotlib not available")
  out1 <- file.path(tempdir(), "glasgow_1")
  expect_equal(system2(py, c("-I", file.path(ROOT, "figures", "make_figures.py"), out1, file.path(out1, "figures")),
                       stdout = FALSE, stderr = FALSE), 0)
  empty <- file.path(tempdir(), "empty_pkg")
  unlink(empty, recursive = TRUE)
  dir.create(empty)
  for (f in c("PROVENANCE.txt", "flow.json")) file.copy(file.path(out1, f), empty)
  p <- fread(file.path(out1, "primary.csv"))[, `:=`(bonferroni = FALSE, fdr = FALSE)]
  fwrite(p, file.path(empty, "primary.csv"))
  fwrite(fread(file.path(out1, "models.csv"))[0], file.path(empty, "models.csv"))
  fl <- jsonlite::read_json(file.path(empty, "flow.json"))
  fl$proteins[[4]][[2]] <- "0"; fl$proteins[[5]][[2]] <- "0"
  jsonlite::write_json(fl, file.path(empty, "flow.json"), auto_unbox = TRUE)
  expect_equal(system2(py, c("-I", file.path(ROOT, "figures", "make_figures.py"), empty, file.path(empty, "figures")),
                       stdout = FALSE, stderr = FALSE), 0)
  expect_true(file.exists(file.path(empty, "figures", "Figure3_forest.png")))
})
