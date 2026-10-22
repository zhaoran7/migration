# 1. PC1-40 GWAS phenotype files
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

use("ukb")


# 2. Build the final phenotype and covariate files
gwas_pheno <- function() {
    x <- as.data.table(readRDS(file.path(phe, "all.rds")))
    x <- x[as.character(ethnic.c) == "White"]
    file <- file.path(root, "data/ref/townsend_2001.asc")
    h <- strsplit(trimws(readLines(file, n = 6)), "[[:space:]]+")
    h <- setNames(vapply(h, function(v) as.numeric(v[2]), 0), tolower(vapply(h, "[", "", 1)))
    grid <- matrix(scan(file, skip = 6, quiet = TRUE), nrow = h["nrows"], byrow = TRUE)
    available <- function(e, n) {
        j <- floor((e - h["xllcorner"])/h["cellsize"]) + 1L
        i <- h["nrows"] - floor((n - h["yllcorner"])/h["cellsize"])
        ok <- is.finite(i) & is.finite(j) & i >= 1 & i <= nrow(grid) & j >= 1 & j <= ncol(grid)
        v <- rep(NA_real_, length(e))
        v[ok] <- grid[cbind(i[ok], j[ok])]
        is.finite(v) & v != h["nodata_value"]
    }
    common <- available(x$birth_east, x$birth_north) & available(x$home_east, x$home_north)
    cov <- c("age", "age2", "sex", "center", "gen_batch", paste0("PC", 1:40))
    x[, `:=`(age = as.numeric(age), age2 = as.numeric(age)^2, sex = as.numeric(sex), gen_batch = fifelse(is.na(gen_batch), 
        NA_character_, paste0("b", gen_batch)), migration_distance = ukb_md(x), delta_tdi_ukb_cond = as.numeric(tdi) - 
        as.numeric(tdi.birth), tdi_birth = as.numeric(tdi.birth))]
    for (trait in c("migration_distance", "delta_tdi_ukb_cond")) {
        req <- c(cov, trait, if (trait == "delta_tdi_ukb_cond") "tdi_birth")
        a <- x[complete.cases(x[, ..req])]
        if (trait == "migration_distance") 
            a <- a[migration_distance > 0]
        else a <- a[eid %in% x$eid[common]]
        a <- a[, c("eid", req), with = FALSE]
        setnames(a, "eid", "IID")
        a[, `:=`(FID, IID)]
        setcolorder(a, c("FID", "IID", req))
        fwrite(a, file.path(phe, paste0(trait, ".pheno")), sep = "\t", na = "NA")
    }
}

if (sys.nframe() == 0L) gwas_pheno()

