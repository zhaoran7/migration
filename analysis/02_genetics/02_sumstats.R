# 1. GWAS merge and harmonization
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))


# 2. Merge completed chromosomes
regenie_merge <- function(folder, trait, name) {
    rows <- lapply(1:22, function(chr) {
        file <- file.path(folder, paste0("chr", chr, "_", trait, ".regenie"))
        x <- fread(file)
        stopifnot(all(c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "LOG10P") %in% 
            names(x)))
        if ("TEST" %in% names(x)) 
            x <- x[TEST == "ADD"]
        x <- x[is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(LOG10P) & LOG10P >= 0 & is.finite(N) & N > 
            0]
        stopifnot(nrow(x) > 0, all(x$CHROM == chr))
        x[, .(CHR = CHROM, POS = GENPOS, SNP = ID, EA = ALLELE1, NEA = ALLELE0, EAF = A1FREQ, BETA, SE, P = pmax(10^(-LOG10P), 
            1e-300), N)]
    })
    x <- rbindlist(rows)
    ids <- lapply(1:22,function(chr)fread(file.path(folder,paste0('chr',chr,'_',trait,'.regenie.ids')),header=FALSE))
    stopifnot(all(vapply(ids,function(a)identical(as.character(a[[2]]),as.character(ids[[1]][[2]])),TRUE)))
    fwrite(ids[[1]],file.path(phe,paste0(trait,'.ids')),sep='\t',col.names=FALSE)
    destination <- file.path(out, "gwas/raw")
    dir.create(destination, recursive = TRUE, showWarnings = FALSE)
    fwrite(x, file.path(destination, paste0(name, ".qt.gz")), sep = "\t")
}


# 3. Harmonize variants and alleles
sumstats_clean <- function(file, map, output) {
    x <- fread(file)
    a <- fread(map)
    stopifnot(all(c("SNP", "CHR", "POS") %in% names(a)), !anyDuplicated(a$SNP))
    x <- x[is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(P) & P > 0 & P <= 1 & is.finite(N) & N > 0]
    x <- x[is.finite(EAF) & pmin(EAF, 1 - EAF) > 0.01]
    x[, `:=`(c("CHR", "POS"), NULL)]
    x <- merge(a[, .(SNP, CHR, POS)], x, by = "SNP", sort = FALSE)
    setcolorder(x, c("CHR", "POS", "SNP", "EA", "NEA", "EAF", "BETA", "SE", "P", "N"))
    setorder(x, CHR, POS)
    fwrite(x, output, sep = "\t")
}

if (sys.nframe() == 0L) {
    a <- commandArgs(TRUE)
    if (a[1] == "regenie") 
        regenie_merge(a[2], a[3], a[4])
    else if (a[1] == "clean") 
        sumstats_clean(a[2], a[3], a[4])
    else stop("Use regenie <folder> <trait> <name>, or clean <input> <hg38 variant map> <output>")
}

