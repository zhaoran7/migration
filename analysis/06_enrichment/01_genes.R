# 1. MAGMA top-gene sets
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

if (sys.nframe() == 0L) {
    a <- commandArgs(TRUE)
    stopifnot(length(a) == 3L)
    md <- fread(a[1])
    tdi <- fread(a[2])
    convert <- function(x, name) {
        x <- x[,.(GENE=sub('\\.[0-9]+$','',as.character(GENE)),Z=ZSTAT)]
        stopifnot(all(startsWith(x$GENE,'ENSG')))
        x <- x[is.finite(Z),.(Z=max(Z)),by=GENE]
        setnames(x, "Z", name)
        x
    }
    x <- merge(convert(md, "MD"), convert(tdi, "DS"), by = "GENE", all = FALSE)
    fwrite(x, paste0(a[3], ".tsv"), sep = "\t")
    status <- system2("scdrs", c("munge-gs", "--zscore-file", shQuote(paste0(a[3], ".tsv")), "--out-file", shQuote(paste0(a[3], 
        ".gs")), "--weight", "zscore", "--n-max", "1000"))
    stopifnot(status == 0L)
}

