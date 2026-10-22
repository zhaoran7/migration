# 1. Chromosome PGS aggregation
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

scores <- function(folder) {
    meta <- fread(file.path(code, "reference/pgs.tsv"))
    if (!file.exists(file.path(folder, "chr1.sscore"))) {
        x <- NULL
        for (trait in meta$score) {
            a <- NULL
            for (chr in 1:22) {
                d <- fread(file.path(folder, trait, paste0("chr", chr, ".sscore")))
                setnames(d,sub('^#','',names(d)))
                col <- grep("_SUM$", names(d), value = TRUE)
                stopifnot(length(col) == 1L)
                d <- d[, .(eid = as.character(IID), value = get(col))]
                if (is.null(a)) 
                  a <- d
                else {
                  stopifnot(setequal(a$eid, d$eid))
                  a$value <- a$value + d$value[match(a$eid, d$eid)]
                }
            }
            setnames(a, "value", trait)
            x <- if (is.null(x)) 
                a
            else merge(x, a, by = "eid", sort = FALSE)
        }
        x[, `:=`((meta$score), lapply(.SD, z)), .SDcols = meta$score]
        attr(x, "source") <- meta
        saveRDS(x, file.path(phe, "prs.rds"), compress = "gzip")
        return(invisible(x))
    }
    read <- function(chr) {
        x <- fread(file.path(folder, paste0("chr", chr, ".sscore")))
        setnames(x, sub("^#", "", names(x)))
        columns <- paste0(meta$score, "_SUM")
        stopifnot(all(c("IID", columns) %in% names(x)), !anyDuplicated(x$IID))
        x <- x[, c("IID", columns), with = FALSE]
        setnames(x, c("IID", columns), c("eid", meta$score))
        x[, `:=`(eid, as.character(eid))]
        x
    }
    x <- read(1)
    for (chr in 2:22) {
        a <- read(chr)
        stopifnot(setequal(x$eid, a$eid))
        a <- a[match(x$eid, a$eid)]
        for (v in meta$score) set(x, j = v, value = x[[v]] + a[[v]])
    }
    x[, `:=`((meta$score), lapply(.SD, z)), .SDcols = meta$score]
    attr(x, "source") <- meta
    saveRDS(x, file.path(phe, "prs.rds"), compress = "gzip")
}

if (sys.nframe() == 0L) scores(commandArgs(TRUE)[1])

