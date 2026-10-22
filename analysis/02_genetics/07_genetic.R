# 1. GWAS diagnostics and comparisons
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))


# 2. Summarize GWAS quality control
ed3_qc <- function() {
    gwas <- c(md = "migration_distance_pc40", dtdi = "delta_tdi_ukb_cond")
    ld <- fread(file.path(root, "data/ldsc/latest/all.h2.res"))
    qq <- list()
    metrics <- list()
    for (id in names(gwas)) {
        name <- gwas[[id]]
        message("QC: ", name)
        file <- file.path(root, "data/gwas/clean", paste0(name, ".gz"))
        cmd <- sprintf(paste0("gzip -cd %s | awk 'BEGIN{FS=OFS=\"\\t\"} ", "NR==1{for(i=1;i<=NF;i++)h[$i]=i;print \"P\",\"N\";next} ", 
            "{print $h[\"P\"],$h[\"N\"]}'"), shQuote(file))
        x <- fread(cmd = cmd, showProgress = FALSE)
        valid <- is.finite(x$P) & x$P > 0 & x$P <= 1 & is.finite(x$N) & x$N > 0
        n <- sum(valid)
        samples <- unique(x$N[valid])
        stopifnot(n > 0, length(samples) == 1L)
        p <- sort(x$P[valid])
        ranks <- unique(c(seq_len(min(n, 2000L)), round(exp(seq(log(min(n, 2001L)), log(n), length.out = 1200)))))
        qq[[id]] <- data.table(rank = ranks, variants = n, p = p[ranks], expected = -log10((ranks - 0.5)/n), observed = -log10(p[ranks]), 
            null_lower = -log10(qbeta(0.975, ranks, n - ranks + 1)), null_upper = -log10(qbeta(0.025, ranks, n - 
                ranks + 1)))
        fit <- ld[GWAS == name]
        stopifnot(nrow(fit) == 1L)
        cojo <- fread(file.path(root, "data/cojo", name, paste0(name, ".jma.cojo")))
        metrics[[id]] <- data.table(exposure = trait[[id]], gwas = name, n = samples, variants = n, invalid = nrow(x) - 
            n, significant = sum(p < 5e-08), cojo_signals = sum(cojo$pJ < 5e-08, na.rm = TRUE), lambda_all = mean(qchisq(p[c(floor((n + 
            1)/2), ceiling((n + 1)/2))], df = 1, lower.tail = FALSE))/qchisq(0.5, df = 1), lambda_hm3 = fit$lambda_gc, 
            intercept = fit$intercept, intercept_se = fit$intercept_se, attenuation_ratio = fit$ratio, ratio_se = fit$ratio_se, 
            h2 = fit$h2, h2_se = fit$h2_se, covariates = paste("Age, age squared, sex, centre, genotyping batch, PC1-PC40", 
                if (id == "dtdi") 
                  "and Birthplace TDI"
                else ""), qq_definition = paste("All valid GWAS p values ranked before display thinning.", "Expected p=(rank-0.5)/number of valid variants.", 
                "Pointwise 95% beta-order-statistic envelope assumes independent uniform null p values;", "LD means this is a reference envelope, not a simultaneous calibration test."))
        rm(x, p)
        invisible(gc())
    }
    list(qq = qq, metrics = rbindlist(metrics))
}

ed3_read <- function(file, columns, keep, coordinates = FALSE) {
    ids <- tempfile()
    on.exit(unlink(ids), add = TRUE)
    fwrite(data.table(id = keep), ids, col.names = FALSE)
    key <- if (coordinates) 
        sprintf("$h[\"%s\"] \":\" $h[\"%s\"]", columns[["CHR"]], columns[["POS"]])
    else sprintf("$h[\"%s\"]", columns[["SNP"]])
    fields <- paste(sprintf("$h[\"%s\"]", columns), collapse = ",")
    cmd <- sprintf(paste0("gzip -cd %s | awk -v ids=%s 'BEGIN{FS=OFS=\"\\t\";", "while((getline k < ids)>0)keep[k]=1;close(ids)} ", 
        "NR==1{for(i=1;i<=NF;i++)h[$i]=i;print \"%s\";next} ", "{key=%s;if(key in keep)print %s}'"), shQuote(file), 
        shQuote(ids), paste(names(columns), collapse = "\\t"), key, fields)
    x <- fread(cmd = cmd, showProgress = FALSE)
    stopifnot(identical(names(x), names(columns)), nrow(x) > 0L)
    x[, `:=`(CHR, suppressWarnings(as.integer(sub("^chr", "", as.character(CHR)))))]
    x[CHR %in% 1:22]
}

ed3_lift <- function(x) {
    tmp <- tempfile()
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    bed <- file.path(tmp, "input.bed")
    mapped <- file.path(tmp, "mapped.bed")
    fwrite(x[, .(paste0("chr", CHR), POS - 1L, POS, SNP)], bed, sep = "\t", col.names = FALSE)
    code <- system2(file.path(root, "software/liftover/liftOver"), c("-minMatch=0.95", shQuote(bed), shQuote(file.path(root, 
        "software/liftover/hg19ToHg38.over.chain.gz")), shQuote(mapped), shQuote(file.path(tmp, "unmapped.bed"))))
    stopifnot(code == 0L)
    y <- fread(mapped, col.names = c("CHR", "start", "POS", "SNP"))
    y <- y[grepl("^chr([1-9]|1[0-9]|2[0-2])$", CHR)]
    y[, `:=`(CHR, as.integer(sub("^chr", "", CHR)))]
    y <- y[!duplicated(SNP) & !duplicated(SNP, fromLast = TRUE), .(CHR, POS, SNP)]
    merge(y, x[, setdiff(names(x), c("CHR", "POS")), with = FALSE], by = "SNP")
}


# 3. Align comparison alleles
ed3_align <- function(x, y) {
    prepare <- function(d) {
        d <- copy(d)
        d[, `:=`(c("EA", "NEA"), lapply(.SD, toupper)), .SDcols = c("EA", "NEA")]
        d <- d[CHR %in% 1:22 & is.finite(POS) & EA %chin% c("A", "C", "G", "T") & NEA %chin% c("A", "C", "G", "T") & 
            EA != NEA & !paste0(EA, NEA) %chin% c("AT", "TA", "CG", "GC") & is.finite(BETA) & is.finite(SE) & SE > 
            0 & is.finite(P) & P > 0 & P <= 1 & is.finite(N) & N > 0]
        pair <- pmin(paste0(d$EA, d$NEA), paste0(d$NEA, d$EA))
        comp <- chartr("ACGT", "TGCA", pair)
        comp <- pmin(comp, paste0(substr(comp, 2, 2), substr(comp, 1, 1)))
        d[, `:=`(alleles, pmin(pair, comp))]
        d <- d[!duplicated(d, by = c("CHR", "POS", "alleles")) & !duplicated(d, by = c("CHR", "POS", "alleles"), 
            fromLast = TRUE)]
        d
    }
    a <- prepare(x)
    b <- prepare(y)
    d <- merge(a, b, by = c("CHR", "POS", "alleles"), suffixes = c("_current", "_reference"))
    same <- with(d, (EA_current == EA_reference & NEA_current == NEA_reference) | (EA_current == chartr("ACGT", 
        "TGCA", EA_reference) & NEA_current == chartr("ACGT", "TGCA", NEA_reference)))
    opposite <- with(d, (EA_current == NEA_reference & NEA_current == EA_reference) | (EA_current == chartr("ACGT", 
        "TGCA", NEA_reference) & NEA_current == chartr("ACGT", "TGCA", EA_reference)))
    stopifnot(nrow(d) > 10000L, all(xor(same, opposite)))
    d[, `:=`(flipped = opposite, z_current = BETA_current/SE_current, z_reference = ifelse(same, 1, -1) * BETA_reference/SE_reference)]
    list(data = d, counts = data.table(current_input = nrow(x), reference_input = nrow(y), current_valid = nrow(a), 
        reference_valid = nrow(b), matched = nrow(d), flipped = sum(opposite)))
}


# 4. Compare published and study GWAS
ed3_comparison <- function() {
    ref <- file.path(root, "data/ref/migration.gz")
    url <- paste0("ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/Migration/Sumstats/", "Sumstats_logdistance_Furuya_2022Mar.txt.gz")
    if (!file.exists(ref)) {
        stage <- tempfile(fileext = ".gz")
        on.exit(unlink(stage), add = TRUE)
        download.file(url, stage, method = "libcurl", mode = "wb")
        stopifnot(system2("gzip", c("-t", shQuote(stage))) == 0L, file.copy(stage, ref))
    }
    hm3 <- scan(file.path(root, "data/ref/gsmap/gsMap_resource/LDSC_resource/", "hapmap3_snps/w_hm3.snpslist"), 
        what = "", quiet = TRUE)
    cols <- c(CHR = "CHR", POS = "POS", SNP = "SNP", EA = "EA", NEA = "NEA", BETA = "BETA", SE = "SE", P = "P", 
        N = "N")
    ans <- list()
    for (id in c("md", "dtdi")) {
        name <- if (id == "md") 
            "migration_distance_pc40"
        else "delta_tdi_ukb_cond"
        message("Z-Z: ", name)
        x <- ed3_read(file.path(root, "data/gwas/clean", paste0(name, ".gz")), cols, hm3)
        if (id == "md") {
            fields <- cols
            fields[c("POS", "EA", "NEA")] <- c("BP", "A1", "A2")
            y <- ed3_lift(ed3_read(ref, fields, hm3))
        }
        else {
            y <- ed3_read(file.path(root, "data/gwas/clean/tdi.gz"), cols, paste(x$CHR, x$POS, sep = ":"), coordinates = TRUE)
        }
        a <- ed3_align(x, y)
        d <- a$data
        limit <- ceiling(max(abs(c(d$z_current, d$z_reference))))
        step <- 2 * limit/90
        d[, `:=`(xbin = pmin(89L, floor((z_reference + limit)/step)), ybin = pmin(89L, floor((z_current + limit)/step)))]
        bins <- d[, .(variants = .N), by = .(xbin, ybin)]
        bins[, `:=`(x = -limit + (xbin + 0.5) * step, y = -limit + (ybin + 0.5) * step, width = step)]
        signals <- d[P_current < 5e-08 | P_reference < 5e-08, .(CHR, POS, SNP = SNP_current, EA = EA_current, NEA = NEA_current, 
            z_current, z_reference, p_current = P_current, p_reference = P_reference, significant_current = P_current < 
                5e-08, significant_reference = P_reference < 5e-08, flipped)]
        metrics <- copy(a$counts)
        metrics[, `:=`(exposure = trait[[id]], concordance = mean(sign(d$z_current) == sign(d$z_reference)), significant_current = sum(d$P_current < 
            5e-08), significant_reference = sum(d$P_reference < 5e-08), significant_both = sum(d$P_current < 5e-08 & 
            d$P_reference < 5e-08), significant_same_direction = sum((d$P_current < 5e-08 | d$P_reference < 5e-08) & 
            sign(d$z_current) == sign(d$z_reference)), current_n = min(d$N_current), reference_n_min = min(d$N_reference), 
            reference_n_max = max(d$N_reference), limit = limit, current_file = file.path("data/gwas/clean", paste0(name, 
                ".gz")), reference = if (id == "md") 
                "Furuya et al. (2023), log-transformed MD"
            else "Neale UKB, Recruitment TDI", reference_file = if (id == "md") 
                "data/ref/migration.gz"
            else "data/gwas/clean/tdi.gz", reference_url = if (id == "md") 
                url
            else "https://www.nealelab.is/uk-biobank", method = paste("Z=beta/SE. Common autosomal HapMap3 single-nucleotide variants.", 
                "GRCh38 position and unordered strand-equivalent allele-pair matching.", "Remove palindromic alleles and all duplicate locus/allele keys;", 
                "flip reference Z to the current effect allele.", if (id == "md") 
                  "Furuya hg19 positions lifted to GRCh38; use published LDSC-adjusted SE."
                else "", "Annotations use the separately estimated LDSC genetic correlation and standard error.", 
                "Studies share UKB participants; this is not independent replication. All matched Z values enter the 90-by-90 bins.", 
                "Coloured points are significant at p<5e-8 in either GWAS; no outlier truncation."))]
        if (id == "dtdi") {
            rg <- fread(file.path(root, "data/ldsc/latest/all.rg.res"))
            rg <- rg[(GWAS1 == name & GWAS2 == "tdi") | (GWAS2 == name & GWAS1 == "tdi")]
            stopifnot(nrow(rg) == 1L)
            metrics[, `:=`(rg = rg$rg, rg_se = rg$se)]
        }
        else {
            rg <- fread(file.path(root, "data/ldsc/latest/md_published.res"))
            stopifnot(nrow(rg) == 1L, is.finite(rg$rg), is.finite(rg$se))
            metrics[, `:=`(rg = rg$rg, rg_se = rg$se)]
        }
        ans[[id]] <- list(bins = bins, signals = signals, metrics = metrics)
        print(metrics[, .(exposure, matched, flipped, rg, rg_se, significant_both)])
        rm(x, y, a, d)
        invisible(gc())
    }
    ans
}

if (sys.nframe() == 0L) {
    saveRDS(list(qc = ed3_qc(), comparison = ed3_comparison()), file.path(out, "genetic.rds"))
}


# 5. Collect association and independent signals
gwas_points <- function(datdir = file.path(root, "data/gwas/clean"), cojodir = file.path(root, "data/cojo"), gwas = c(md = "migration_distance_pc40", 
    dtdi = "delta_tdi_ukb_cond")) {
    pcut <- 0.001
    sig <- 5e-08
    win <- 1e+07
    text_pt <- 72.27/25.4
    chr_len <- data.table::as.data.table(get("chr_lengths", envir = asNamespace("topr")))
    chr_len[, `:=`(CHR, gsub("chr", "", V1))]
    chr_len <- chr_len[CHR %in% as.character(1:22)]
    chr_len[, `:=`(CHR = as.integer(CHR), V2 = as.numeric(V2))]
    data.table::setorder(chr_len, CHR)
    stopifnot(identical(chr_len$CHR, 1:22))
    chr_len[, `:=`(offset = cumsum(data.table::shift(V2, fill = 0)), center = cumsum(data.table::shift(V2, fill = 0)) + 
        V2/2)]
    offset <- setNames(chr_len$offset, chr_len$CHR)
    xmax <- max(chr_len$offset + chr_len$V2)
    read_points <- function(file) {
        if (!file.exists(file)) 
            stop("Missing Manhattan GWAS input: ", file)
        if (any(!nzchar(Sys.which(c("gzip", "awk"))))) 
            stop("Manhattan requires gzip and awk")
        con <- gzfile(file, "rt")
        header <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1]]
        close(con)
        if (!all(c("CHR", "POS", "SNP", "P") %in% header)) 
            stop("Invalid GWAS header: ", file)
        cmd <- sprintf(paste0("gzip -cd %s | awk 'BEGIN{FS=OFS=\"\\t\"} ", "NR==1{for(i=1;i<=NF;i++) h[$i]=i; print \"CHR\",\"POS\",\"SNP\",\"P\"; next} ", 
            "$h[\"CHR\"]~/^[0-9]+$/ && $h[\"CHR\"]>=1 && $h[\"CHR\"]<=22 && ", "$h[\"P\"]!=\"NA\" && $h[\"P\"]!=\"\" && ($h[\"P\"]+0)<=%.10g ", 
            "{print $h[\"CHR\"],$h[\"POS\"],$h[\"SNP\"],$h[\"P\"]}'"), shQuote(file), pcut)
        x <- data.table::fread(cmd = cmd)
        x[, `:=`(CHR = as.integer(CHR), POS = as.numeric(POS), P = as.numeric(P))]
        if (any(!is.finite(x$P) | x$P <= 0 | !is.finite(x$POS) | x$POS < 1)) 
            stop("Nonpositive/nonfinite P or invalid POS in Manhattan candidates: ", file)
        if (!nrow(x)) 
            stop("No Manhattan candidates at P <= 1e-3: ", file)
        x[, `:=`(bin, floor(POS/2e+05))]
        x <- data.table::rbindlist(list(x[P < sig], x[P >= sig][order(P), .SD[1], by = .(CHR, bin)]), fill = TRUE)
        x[, `:=`(x = POS + offset[as.character(CHR)], y = -log10(P), col = factor(CHR%%2), bin = NULL)]
        x
    }
    read_loci <- function(file) {
        if (!file.exists(file)) 
            stop("Missing COJO input (annotations must not silently disappear): ", file)
        x <- data.table::fread(file)
        if (!all(c("Chr", "bp", "SNP", "p") %in% names(x))) 
            stop("Invalid COJO header: ", file)
        data.table::setnames(x, c("Chr", "bp", "SNP", "p"), c("CHROM", "POS", "ID", "P"))
        x[, `:=`(CHROM = as.integer(CHROM), POS = as.numeric(POS), P = as.numeric(P))]
        if (any(!x$CHROM %in% 1:22 | !is.finite(x$POS) | !is.finite(x$P) | x$P <= 0)) 
            stop("Invalid COJO loci")
        x
    }
    read_labels <- function(points, cojo, side) {
        loci <- read_loci(cojo)
        if (!nrow(loci)) 
            return(data.table::data.table())
        match_pos <- match(loci$ID, points$SNP)
        present <- !is.na(match_pos)
        if (any(loci$CHROM[present] != points$CHR[match_pos[present]] | loci$POS[present] != points$POS[match_pos[present]])) 
            stop("GWAS and COJO coordinates differ: ", cojo)
        loci[, `:=`(locus, .I)]
        g <- points[P < sig, .(CHROM = CHR, POS, ID = SNP, P)]
        top <- data.table::rbindlist(lapply(seq_len(nrow(loci)), function(i) {
            z <- g[CHROM == loci$CHROM[i] & abs(POS - loci$POS[i]) <= win/2]
            if (nrow(z)) 
                z[which.min(P)]
            else loci[i, .(CHROM, POS, ID, P)]
        }))
        top[, `:=`(locus = loci$locus, cojo_CHROM = loci$CHROM, cojo_POS = loci$POS, cojo_ID = loci$ID, cojo_P = loci$P)]
        gene <- data.table::as.data.table(topr::annotate_with_nearest_gene(as.data.frame(loci[, .(CHROM, POS, ID, 
            P)]), build = 38, protein_coding_only = TRUE))
        key <- function(d) paste(d$CHROM, d$POS, d$ID, sep = ":")
        j <- match(key(loci), key(gene))
        if (anyNA(j)) 
            stop("Nearest-gene annotation did not return all COJO loci")
        top[, `:=`(Gene_Symbol, gene$Gene_Symbol[j])]
        top[is.na(Gene_Symbol) | Gene_Symbol == "", `:=`(Gene_Symbol, cojo_ID)]
        top[, `:=`(x = POS + offset[as.character(CHROM)], y = -log10(P), side = side)]
        top
    }
    md <- read_points(file.path(datdir, paste0(gwas["md"], ".gz")))
    ds <- read_points(file.path(datdir, paste0(gwas["dtdi"], ".gz")))
    md_lab <- read_labels(md, file.path(cojodir, gwas["md"], paste0(gwas["md"], ".jma.cojo")), 1)
    ds_lab <- read_labels(ds, file.path(cojodir, gwas["dtdi"], paste0(gwas["dtdi"], ".jma.cojo")), -1)
    list(chromosomes = chr_len[, .(V1 = paste0("chr", CHR), V2)], points = list(md = md, dtdi = ds), labels = list(md = md_lab, 
        dtdi = ds_lab))
}

