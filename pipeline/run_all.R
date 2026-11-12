#!/usr/bin/env Rscript
# Run the whole analysis and write one versioned analysis package (aggregate results only).
#
#   Rscript run_all.R --map config/data_map_glasgow.yaml --out outputs/<run_id>     # the group's dataset
#   Rscript run_all.R --synthetic [--layout glasgow|rap]                            # simulated data
#
# Options: --plan <yaml> (default config/analysis_plan.yaml), --map <yaml>, --root <dir> (where the data and code-list
# paths in the map are resolved), --out <dir> (must not exist or be empty), --overwrite (synthetic runs only).

options(warnPartialMatchDollar = TRUE, warnPartialMatchArgs = TRUE, warnPartialMatchAttr = TRUE)
args <- commandArgs(trailingOnly = TRUE)
known <- c("--plan", "--map", "--root", "--out", "--layout", "--synthetic", "--overwrite")
flags <- c("--synthetic", "--overwrite")
i <- 1
opts <- list()
while (i <= length(args)) {
  a <- args[i]
  if (!a %in% known) stop("unknown argument: ", a, call. = FALSE)
  if (!is.null(opts[[a]])) stop("argument given twice: ", a, call. = FALSE)
  if (a %in% flags) { opts[[a]] <- TRUE; i <- i + 1; next }
  if (i == length(args) || startsWith(args[i + 1], "--")) stop("missing value for ", a, call. = FALSE)
  opts[[a]] <- args[i + 1]
  i <- i + 2
}
opt <- function(name, default = NULL) if (is.null(opts[[paste0("--", name)]])) default else opts[[paste0("--", name)]]
here <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
synthetic <- isTRUE(opts[["--synthetic"]])
layout <- opt("layout", "glasgow")
syn_dir <- file.path(here, "synthetic", layout)
plan_path <- opt("plan", file.path(here, "config/analysis_plan.yaml"))
map_path <- opt("map", if (synthetic) file.path(syn_dir, "data_map.yaml") else file.path(here, "config/data_map_glasgow.yaml"))
root <- normalizePath(opt("root", if (synthetic) syn_dir else here))
out <- opt("out", file.path(here, "outputs", if (synthetic) paste0("synthetic_", layout) else format(Sys.Date(), "run_%Y%m%d")))
if (isTRUE(opts[["--overwrite"]]) && !synthetic) stop("--overwrite is only allowed for synthetic runs", call. = FALSE)

for (f in c("io.R", "cohort.R", "models.R", "secondary.R", "legacy.R", "export.R", "pipeline.R")) source(file.path(here, "R", f))

plan <- yaml::read_yaml(plan_path)
map <- yaml::read_yaml(map_path)
res <- run_analysis(plan, map, root, synthetic = synthetic)
write_package(out, res, plan, here, synthetic, config_paths = c(plan_path, map_path), overwrite = isTRUE(opts[["--overwrite"]]),
              map = map)
if (length(res$unresolved)) message("unresolved decisions (the run is DRAFT): ", paste(res$unresolved, collapse = ", "))
message("analysis package written to ", out)
