# 1. Protein and MRI input matrices
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

use("ukb")


# 2. Read the fixed feature list and save the matrix
omics_import <- function(layer, file) {
    info <- fread(file.path(code, "reference", paste0(tolower(layer), ".tsv")))
    x <- if (grepl("\\.rds$", file)) 
        as.data.table(readRDS(file))
    else fread(file)
    stopifnot("eid" %in% names(x), !anyDuplicated(x$eid))
    ids <- info[["Feature ID"]]
    if (layer == "MRI" && !all(ids %in% names(x))) {
        fields <- paste0(info[["UKB field ID"]], "_i2")
        if (all(paste0("img_", info[["UKB field ID"]]) %in% names(x))) 
            fields <- paste0("img_", info[["UKB field ID"]])
        stopifnot(all(fields %in% names(x)))
        setnames(x, fields, ids)
    }
    if (layer == "Protein") {
        old <- names(x)[toupper(names(x)) %in% ids]
        setnames(x, old, toupper(old))
    }
    stopifnot(all(ids %in% names(x)))
    x <- x[, c("eid", ids), with = FALSE]
    x[, `:=`(eid, as.character(eid))]
    x[, `:=`((ids), lapply(.SD, function(v) {
        v <- as.numeric(v)
        v[!is.finite(v)] <- NA_real_
        v
    })), .SDcols = ids]
    info[, `:=`(`N non-missing`, vapply(x[, ..ids], function(v) sum(is.finite(v)), 0L))]
    info[, `:=`(`Missing rate`, 1 - `N non-missing`/nrow(x))]
    info[, `:=`(Mean, vapply(x[, ..ids], mean, 0, na.rm = TRUE))]
    info[, `:=`(SD, vapply(x[, ..ids], sd, 0, na.rm = TRUE))]
    name <- if (layer == "Protein") 
        "prot"
    else "mri"
    saveRDS(x, file.path(phe, paste0(name, ".rds")), compress = "gzip")
    path <- file.path(phe, "omics_info.rds")
    meta <- if (file.exists(path)) 
        readRDS(path)
    else list()
    meta[[layer]] <- info
    saveRDS(meta, path, compress = "gzip")
}

if (sys.nframe() == 0L) {
    args <- commandArgs(TRUE)
    stopifnot(length(args) == 2L, args[1] %in% c("Protein", "MRI"))
    omics_import(args[1], args[2])
}

