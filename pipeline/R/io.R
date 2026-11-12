# Reading the configuration and input files, and recoding variables as specified in the data map.
# Specs are read with [[ ]]: `$` would partially match (spec$columns would return columns_prefix).

suppressPackageStartupMessages(library(data.table))

read_table <- function(path) {
  if (!file.exists(path)) stop("input file not found: ", path, call. = FALSE)
  if (grepl("\\.dta$", path, ignore.case = TRUE)) {
    if (!requireNamespace("haven", quietly = TRUE)) stop("reading .dta files needs the haven package", call. = FALSE)
    return(as.data.table(haven::zap_labels(haven::read_dta(path))))
  }
  fread(path, na.strings = c("", "NA"), showProgress = FALSE, encoding = "UTF-8")
}

read_codes <- function(path) {
  if (!file.exists(path)) stop("code list not found: ", path, call. = FALSE)
  x <- trimws(sub("^﻿", "", readLines(path, warn = FALSE, encoding = "UTF-8")))
  x <- x[nzchar(x) & !startsWith(x, "#")]
  if (!length(x)) stop("code list is empty: ", path, " (fill it in before running)", call. = FALSE)
  v <- suppressWarnings(as.numeric(x))
  if (anyNA(v)) stop("non-numeric code(s) in ", path, ": ", paste(x[is.na(v)], collapse = ", "), call. = FALSE)
  v
}

num_or <- function(x, default) {
  if (is.null(x) || identical(x, "NA")) default else as.numeric(x)
}

codes_of <- function(spec, key) {
  v <- as.numeric(unlist(spec[[key]]))
  if (anyNA(v)) stop("codes in '", key, "' must be numbers", call. = FALSE)
  v
}

# Raised when a column is genuinely absent, so that only this case can be excused for optional variables.
missing_column <- function(name, cols) {
  structure(class = c("missing_column", "error", "condition"),
            list(message = sprintf("variable '%s': column(s) not found: %s", name, paste(cols, collapse = ", ")),
                 call = NULL))
}

spec_columns <- function(dt, spec, name) {
  cols <- if (!is.null(spec[["columns"]])) unlist(spec[["columns"]]) else
    grep(paste0("^", spec[["columns_prefix"]], "[0-9]+$"), names(dt), value = TRUE)
  missing <- setdiff(cols, names(dt))
  if (!length(cols)) stop(missing_column(name, spec[["columns_prefix"]]))
  if (length(missing)) stop(missing_column(name, missing))
  cols
}

as_number <- function(x, name) {
  if (is.numeric(x)) return(as.numeric(x))
  v <- suppressWarnings(as.numeric(as.character(x)))
  bad <- !is.na(x) & is.na(v)
  if (any(bad)) stop(sprintf("variable '%s': %d non-numeric value(s), e.g. '%s'", name, sum(bad), x[bad][1]), call. = FALSE)
  v
}

# Precedence when several columns hold codes: codes_1 > codes_0 > codes_na > other; all columns empty -> none.
recode_codes <- function(dt, spec, name, root) {
  cols <- spec_columns(dt, spec, name)
  true <- if (!is.null(spec[["codes_1_file"]])) read_codes(file.path(root, spec[["codes_1_file"]])) else codes_of(spec, "codes_1")
  false <- codes_of(spec, "codes_0")
  na <- codes_of(spec, "codes_na")
  if (!length(true) && !length(false)) stop(sprintf("variable '%s': no codes_1 or codes_0 given", name), call. = FALSE)
  m <- sapply(cols, function(cl) as_number(dt[[cl]], name))
  if (!is.matrix(m)) m <- matrix(m, nrow = nrow(dt))
  hit <- function(codes) {
    if (!length(codes)) return(rep(FALSE, nrow(m)))
    rowSums(matrix(!is.na(m) & m %in% codes, nrow(m))) > 0
  }
  other <- rowSums(matrix(!is.na(m) & !(m %in% c(true, false, na)), nrow(m))) > 0
  out <- rep(NA_real_, nrow(m))
  out[other] <- num_or(spec[["other"]], NA_real_)
  out[hit(na)] <- NA_real_
  out[hit(false)] <- 0
  out[hit(true)] <- 1
  out[rowSums(!is.na(m)) == 0] <- num_or(spec[["none"]], NA_real_)
  out
}

recode_variable <- function(dt, spec, name, root) {
  g <- function(k) spec[[k]]
  col <- function(x) {
    if (!x %in% names(dt)) stop(missing_column(name, x))
    as_number(dt[[x]], name)
  }
  drop_na_codes <- function(x) { x[x %in% codes_of(spec, "codes_na")] <- NA; x }
  switch(g("type"),
    numeric = drop_na_codes(col(g("column"))),
    precomputed = drop_na_codes(col(g("column"))),
    constant = rep(as.numeric(g("value")), nrow(dt)),
    present = rep(1, nrow(dt)),   # 1 for every row of that file; participants absent from it become 0 after merging
    ratio = {
      num <- col(g("numerator")); den <- col(g("denominator"))
      z <- num / den
      z[!is.finite(z) | den <= 0 | num <= 0] <- NA   # waist or hip of zero or less is not a valid measurement
      z
    },
    factor = {
      lv <- g("levels")
      x <- col(g("column"))
      factor(unlist(lv)[match(as.character(x), names(lv))], levels = unlist(lv))
    },
    factor_observed = {
      if (!g("column") %in% names(dt)) stop(missing_column(name, g("column")))
      x <- dt[[g("column")]]
      factor(x, levels = sort(unique(x[!is.na(x)])))
    },
    split = {
      x <- col(g("column"))
      factor(ifelse(x < g("cut"), g("below"), g("at_or_above")), levels = c(g("below"), g("at_or_above")))
    },
    codes = recode_codes(dt, spec, name, root),
    stop(sprintf("variable '%s': unknown type '%s'", name, g("type")), call. = FALSE)
  )
}

check_ids <- function(x, what, id) {
  if (!id %in% names(x)) stop("id column '", id, "' missing in the ", what, " file", call. = FALSE)
  ids <- as.character(x[[id]])
  if (anyNA(ids) || any(!nzchar(ids))) stop(what, " file: missing participant identifiers", call. = FALSE)
  if (anyDuplicated(ids)) stop(what, " file: duplicate participant identifiers (", sum(duplicated(ids)), ")", call. = FALSE)
  ids
}

select_assays <- function(prot, map, id, root, allow_all) {
  spec <- map[["proteins"]][["columns"]]
  if (is.null(spec)) {
    if (!allow_all) stop("the data map must list the assay columns (proteins: columns: from/to, or a manifest)", call. = FALSE)
    return(setdiff(names(prot), id))
  }
  if (!is.null(spec[["manifest"]])) {
    assays <- read_table(file.path(root, spec[["manifest"]]))[["column"]]
  } else {
    pos <- function(v) if (is.numeric(v)) as.integer(v) else match(v, names(prot))
    from <- pos(spec[["from"]]); to <- pos(spec[["to"]])
    if (anyNA(c(from, to)) || from > to || to > ncol(prot)) stop("invalid assay column range in the data map", call. = FALSE)
    assays <- names(prot)[from:to]
    for (end in c("first", "last")) {
      want <- spec[[end]]
      got <- if (end == "first") assays[1] else assays[length(assays)]
      if (!is.null(want) && !identical(toupper(got), toupper(want)))
        stop(sprintf("assay range %s column is '%s', not '%s': the column positions are shifted", end, got, want), call. = FALSE)
    }
    if (!allow_all && (is.null(spec[["first"]]) || is.null(spec[["last"]])))
      stop("set proteins: columns: first/last in the data map (names of the first and last assay columns)", call. = FALSE)
  }
  if (!length(assays) || anyDuplicated(assays)) stop("assay list is empty or has duplicates", call. = FALSE)
  if (length(setdiff(assays, names(prot)))) stop("assays missing from the protein file", call. = FALSE)
  if (id %in% assays) stop("the id column falls inside the assay range", call. = FALSE)
  if (!is.null(spec[["check_names"]])) {   # guard against a shifted column range (e.g. a row-name column)
    known <- toupper(read_table(file.path(root, spec[["check_names"]]))[["Assay"]])
    unknown <- assays[!toupper(assays) %in% known]
    if (length(unknown)) stop(length(unknown), " selected column(s) are not Olink assay names, e.g. ",
                              paste(head(unknown, 3), collapse = ", "), call. = FALSE)
  }
  nonnum <- assays[!vapply(prot[, ..assays], is.numeric, logical(1))]
  if (length(nonnum)) stop("non-numeric assay column(s): ", paste(head(nonnum), collapse = ", "), call. = FALSE)
  assays
}

# Target names for the selected assays (manifest "target" column if given, else the column names themselves).
labels_from_manifest <- function(map, root, assays) {
  m <- map[["proteins"]][["columns"]][["manifest"]]
  if (is.null(m)) return(setNames(assays, assays))
  man <- read_table(file.path(root, m))
  lab <- if ("target" %in% names(man)) man$target[match(assays, man$column)] else assays
  setNames(ifelse(is.na(lab), assays, lab), assays)
}

# Returns list(data = one row per participant with Olink data (analysis variables + assay columns),
#              assays, withdrawn ids, kinship pairs or NULL, inputs = file checksums).
load_inputs <- function(map, root, synthetic = FALSE) {
  id <- map[["id"]]
  path <- function(key) if (is.null(map$files[[key]])) NULL else file.path(root, map$files[[key]])
  cache <- list()
  get_file <- function(key) {
    p <- path(key)
    if (is.null(p)) stop("no file '", key, "' in the data map", call. = FALSE)
    if (is.null(cache[[p]])) {
      x <- read_table(p)
      check_ids(x, key, id)
      x[, (id) := as.character(get(id))]
      cache[[p]] <<- x
    }
    cache[[p]]
  }
  prot <- copy(get_file("proteins"))
  assays <- select_assays(prot, map, id, root, allow_all = synthetic)

  vars <- data.table(eid = prot[[id]])
  for (name in names(map$variables)) {
    spec <- map$variables[[name]]
    key <- if (is.null(spec[["file"]])) "participants" else spec[["file"]]
    v <- tryCatch({
      src <- get_file(key)
      data.table(eid = src[[id]], v = recode_variable(src, spec, name, root))
    }, missing_column = function(e) {
      if (!isTRUE(spec[["optional"]])) stop(e)
      warning("optional variable '", name, "' skipped: ", conditionMessage(e), call. = FALSE)
      NULL
    })
    if (is.null(v)) next
    setnames(v, "v", name)
    vars <- merge(vars, v, by = "eid", all.x = TRUE, sort = FALSE)
    if (identical(spec[["type"]], "present")) vars[is.na(get(name)), (name) := 0]
  }
  setnames(prot, id, "eid")
  d <- merge(vars, prot[, c("eid", assays), with = FALSE], by = "eid", all.x = TRUE, sort = TRUE)
  if (nrow(d) != nrow(prot)) stop("row count changed while merging inputs", call. = FALSE)

  # Withdrawals: UK Biobank's list (file) and/or negative ids in the extract. A real run must have one of them.
  withdrawn <- character()
  wp <- path("withdrawn")
  if (!is.null(wp)) {
    if (!file.exists(wp)) {
      if (!synthetic) stop("withdrawal list not found: ", wp, call. = FALSE)
    } else {
      raw <- trimws(sub("^﻿", "", fread(wp, header = FALSE, colClasses = "character")[[1]]))
      withdrawn <- unique(raw[grepl("^-?[0-9]+$", raw)])
      if (!length(withdrawn)) stop("withdrawal list contains no ids: ", wp, call. = FALSE)
    }
  }
  if (isTRUE(map[["withdrawn_negative_ids"]])) withdrawn <- union(withdrawn, d$eid[startsWith(d$eid, "-")])
  if (!synthetic && is.null(wp) && !isTRUE(map[["withdrawn_negative_ids"]]))
    stop("no withdrawal list and no negative-id rule: withdrawn participants must be removed", call. = FALSE)

  kin <- NULL
  if (identical(map$relatedness$method, "kinship_pairs")) {
    kin <- read_table(path("kinship"))
    if (anyNA(kin$Kinship)) stop("kinship file has missing kinship coefficients", call. = FALSE)
    kin <- kin[Kinship >= map$relatedness$min_kinship, .(ID1 = as.character(ID1), ID2 = as.character(ID2))]
  }
  files <- unique(unlist(lapply(names(map$files), path)))
  files <- files[file.exists(files)]
  inputs <- vapply(files, function(p) digest::digest(file = p, algo = "sha256"), "")
  names(inputs) <- basename(files)
  list(data = d, assays = assays, withdrawn = withdrawn, kinship = kin, inputs = inputs,
       labels = labels_from_manifest(map, root, assays))
}
