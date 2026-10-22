# 1. Independent COJO signals
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

if (sys.nframe() == 0L) {
    a <- commandArgs(TRUE)
    stage <- a[1]
    trait <- a[2]
    folder <- file.path(out, "cojo", trait)
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    if (stage == "prepare") {
        x <- fread(file.path(root, "data/gwas/clean", paste0(trait, ".gz")))
        x <- x[CHR %in% 1:22 & is.finite(EAF) & pmin(EAF, 1 - EAF) > 0.01 & is.finite(BETA) & is.finite(SE) & SE > 
            0 & is.finite(P) & P > 0 & P <= 1 & is.finite(N) & N > 0]
        for (chr in 1:22) fwrite(x[CHR == chr, .(SNP, A1 = EA, A2 = NEA, freq = EAF, b = BETA, se = SE, p = P, 
            N)], file.path(folder, paste0("chr", chr, ".ma")), sep = "\t")
    }
    else if (stage == "merge") {
        files <- file.path(folder, paste0("chr", 1:22, ".jma.cojo"))
        logs <- file.path(folder, paste0("chr", 1:22, ".log"))
        stopifnot(all(file.exists(logs)), all(vapply(logs, function(f) any(grepl("Analysis finished", readLines(f))), 
            TRUE)))
        present <- files[file.exists(files)]
        x <- if(length(present))rbindlist(lapply(present,fread),fill=TRUE) else
          data.table(Chr=integer(),SNP=character(),bp=integer(),refA=character(),freq=double(),
            b=double(),se=double(),p=double(),n=double(),freq_geno=double(),
            bJ=double(),bJ_se=double(),pJ=double(),LD_r=double())
        if(nrow(x)) {
          map <- fread(file.path(root,'data/gwas/clean',paste0(trait,'.gz')),select=c('SNP','CHR','POS'))
          map <- unique(map)
          stopifnot(!anyDuplicated(map$SNP))
          j <- match(x$SNP,map$SNP);stopifnot(!anyNA(j))
          x[,`:=`(Chr=map$CHR[j],bp=map$POS[j])]
        }
        fwrite(x, file.path(folder, paste0(trait, ".jma.cojo")), sep = "\t")
    }
    else stop("Use prepare or merge, followed by the trait name")
}

