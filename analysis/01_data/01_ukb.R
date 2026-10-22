# 1. UKB extraction and analysis variables
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

ukb_check <- function(x, columns = "eid") {
    stopifnot(!anyDuplicated(names(x)), all(columns %in% names(x)))
    x$eid <- as.character(x$eid)
    stopifnot(!anyNA(x$eid), all(nzchar(x$eid)), !anyDuplicated(x$eid))
    x
}

ukb_fields <- function(map = FALSE) {
    v <- c(p31 = "sex", p53 = "date_attend", p54 = "center", p129 = "birth_north", p130 = "birth_east", p191 = "date_lost", 
        p738 = "inc", p767 = "job_week_hours", p806 = "job_stand", p816 = "job_heavy", p826 = "job_shift", p3426 = "job_night", 
        p6138 = "edu", p20116 = "smoke_status", p20117 = "alcohol_status", p20277 = "job_code", p21000 = "ethnicity", 
        p21001 = "bmi", p21003 = "age", p22000 = "gen_batch", p22189 = "tdi", p22599 = "job_num", p22688 = "home_east", 
        p22689 = "home_north")
    if (map) 
        v
    else unname(v)
}

ukb_header <- function(path) {
    if (is.data.frame(path)) 
        return(setNames(names(path), names(path)))
    con <- if (grepl("\\.gz$", path)) 
        gzfile(path)
    else file(path)
    on.exit(close(con))
    h <- strsplit(sub("\r$", "", readLines(con, n = 1L)), "\t", fixed = TRUE)[[1]]
    key <- sub("^(?:f\\.)?([0-9]+)[.-]([0-9]+)\\.([0-9]+)$", "p\\1_i\\2_a\\3", h, perl = TRUE)
    stopifnot(!anyDuplicated(key))
    setNames(h, key)
}

ukb_read <- function(path, columns) {
    if (is.data.frame(path)) 
        return(ukb_check(path[, columns, drop = FALSE]))
    h <- ukb_header(path)
    if (!all(columns %in% names(h))) 
        stop("Missing fields in ", path, ": ", paste(setdiff(columns, names(h)), collapse = ", "))
    x <- as.data.frame(data.table::fread(path, select = unname(h[columns]), sep = "\t", na.strings = c("", "NA"), 
        showProgress = FALSE, nThread = 1L))
    names(x) <- names(h)[match(names(x), h)]
    ukb_check(x)
}

ukb_baseline <- function(h, field) {
    v <- intersect(c(field, paste0(field, "_i0"), paste0(field, "_i0_a0")), names(h))
    if (length(v) != 1L) 
        stop("Ambiguous/missing baseline field: ", field)
    v
}

ukb_array <- function(h, field) {
    v <- grep(paste0("^", field, "(_i0)?_a[0-9]+$"), names(h), value = TRUE)
    a <- as.integer(sub(".*_a", "", v))
    stopifnot(length(v) > 0L, !anyDuplicated(a))
    v[order(a)]
}


# 2. Read field-mapped extracts
ukb_vip <- function(path, map = ukb_fields(TRUE)) {
    h <- ukb_header(path)
    v <- vapply(names(map), function(f) ukb_baseline(h, f), character(1))
    q <- grep('^p6138_i0_a[0-9]+$', names(h), value = TRUE)
    q <- q[order(as.integer(sub('.*_a', '', q)))]
    x <- ukb_read(path, unique(c("eid", v, q)))
    edu <- if (length(q) > 1L) apply(x[q], 1, function(a) {
        a <- a[!is.na(a) & nzchar(a)]
        if (length(a)) paste(a, collapse = '|') else NA_character_
    }) else as.character(x[[v[map == 'edu']]])
    x <- x[, c('eid', v), drop = FALSE]
    names(x)[match(v, names(x))] <- unname(map)
    x$edu <- edu
    for (j in c("date_attend", "date_lost")) x[[j]] <- as.Date(x[[j]])
    x$edu <- as.character(x$edu)
    x$edu[!is.na(x$edu) & x$edu == "-7"] <- "7"
    for (j in setdiff(names(x), "eid")) {
        n <- suppressWarnings(as.numeric(x[[j]]))
        bad <- !is.na(n) & n < 0
        if (any(bad) && all(n[bad] %in% c(-1, -2, -3, -7, -10, -121, -818))) 
            x[[j]][bad] <- NA
    }
    dates <- as.Date(c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07"))
    for (j in c("date_attend", "date_lost")) x[[j]][x[[j]] %in% dates] <- as.Date(NA)
    x
}


# 3. Derive baseline and geographic variables
ukb_derive <- function(x) {
    x <- ukb_check(x, c("eid", ukb_fields()))
    x$center <- paste0("c", x$center)
    e <- x$ethnicity
    e <- ifelse(e == 2002, 12, ifelse(e == 2003, 13, ifelse(e %in% c(2004, 3004), 6, ifelse(grepl("^1", e), 1, 
        ifelse(grepl("^2", e), 2, ifelse(grepl("^3", e), 3, ifelse(grepl("^4", e), 4, e)))))))
    x$ethnic.c <- factor(e, c(1:6, 12:13), c("White", "Mixed", "Asian", "Black", "Chinese", "Other", "White-Black", 
        "White-Asian"))
    x$edu.b <- factor(ifelse(is.na(x$edu), NA, as.integer(grepl("1", x$edu))))
    code <- vapply(strsplit(as.character(x$edu), "|", fixed = TRUE), function(v) {
        v <- v[nzchar(v)]
        if (length(v)) 
            max(suppressWarnings(as.numeric(v)))
        else NA_real_
    }, numeric(1))
    x$edu.years <- unname(c(`1` = 20, `2` = 13, `3` = 10, `4` = 10, `5` = 19, `6` = 15, `7` = 7)[as.character(code)])
    x$edu.sco <- as.integer(findInterval(x$edu.years, c(0, 3, 6, 9, 12, 15, 17, 19, 21, 23)))
    x$smoke.c <- factor(x$smoke_status, 0:2, c("never", "previous", "current"))
    x$alcohol.c <- factor(x$alcohol_status, 0:2, c("never", "previous", "current"))
    x
}

ukb_md <- function(x) sqrt((as.double(x$home_east) - x$birth_east)^2 + (as.double(x$home_north) - x$birth_north)^2)/1000

ukb_geo <- function(x) {
    if (!requireNamespace("RANN", quietly = TRUE)) 
        stop("Install RANN for nearest-TDI reconstruction")
    d <- x[!is.na(x$home_east) & !is.na(x$birth_east), , drop = FALSE]
    stopifnot(nrow(d) > 0L, all(is.finite(as.matrix(d[, c("home_east", "home_north", "birth_east", "birth_north")]))))
    ref <- d[!is.na(d$tdi), c("home_east", "home_north", "tdi"), drop = FALSE]
    ref <- ref[!duplicated(ref[1:2]), , drop = FALSE]
    stopifnot(nrow(ref) > 0L)
    nn <- RANN::nn2(as.matrix(ref[1:2]), as.matrix(d[c("birth_east", "birth_north")]), k = 1)
    data.frame(eid = d$eid, move = ukb_md(d), tdi.birth = ref$tdi[nn$nn.idx[, 1]], tdi.birth.nn = nn$nn.dists[, 
        1])
}

ukb_pc <- function(path) {
    v <- ukb_array(ukb_header(path), "p22009")
    a <- as.integer(sub(".*_a", "", v))
    stopifnot(length(v) >= 40L, a[1] %in% 0:1, identical(a[1:40], seq.int(a[1], length.out = 40)))
    x <- ukb_read(path, c("eid", v[1:40]))
    names(x)[match(v[1:40], names(x))] <- paste0("PC", 1:40)
    x
}


# 4. Link death and hospital dates
ukb_death <- function(path) {
    h <- ukb_header(path)
    dcols <- grep("^p40000_i[01](_a0)?$", names(h), value = TRUE)
    ccols <- grep("^p40001_i[0-9]+(_a0)?$", names(h), value = TRUE)
    stopifnot(length(dcols) > 0L, length(ccols) > 0L)
    dcols <- dcols[order(as.integer(sub("^p40000_i([0-9]+).*", "\\1", dcols)))]
    x <- ukb_read(path, c("eid", dcols, ccols))
    d <- as.Date(x[[dcols[1]]])
    for (v in dcols[-1]) {
        a <- as.Date(x[[v]])
        d[is.na(d)] <- a[is.na(d)]
    }
    d[!is.na(d) & d > as.Date("2023-04-01")] <- as.Date("2023-04-01")
    cause <- apply(x[ccols], 1, function(v) {
        v <- v[!is.na(v) & v != ""]
        if (length(v)) 
            paste(v, collapse = " ")
        else NA_character_
    })
    cvd <- d
    cvd[is.na(d) | !grepl("^I", cause)] <- as.Date(NA)
    data.frame(eid = x$eid, date_death = d, death_icd10 = cause, date_cvd_death = cvd)
}

ukb_icd <- function(code_path, date_path) {
    dv <- ukb_array(ukb_header(date_path), "p41280")
    indices <- as.integer(sub(".*_a", "", dv))
    stopifnot(identical(indices, seq.int(0L, length.out = length(dv))))
    date <- ukb_read(date_path, c("eid", dv))
    h <- ukb_header(code_path)
    packed <- intersect(c("p41270", "p41270_i0"), names(h))
    if (length(packed)) {
        stopifnot(length(packed) == 1L)
        raw <- ukb_read(code_path, c("eid", packed))
        parts <- data.table::tstrsplit(as.character(raw[[packed]]), "|", fixed = TRUE, fill = NA_character_)
        if (length(parts) > length(dv)) 
            stop("ICD10 codes exceed date slots; refusing truncation")
        code <- data.frame(eid = raw$eid)
        for (i in seq_along(dv)) code[[paste0("i", i)]] <- if (i <= length(parts)) 
            parts[[i]]
        else rep(NA_character_, nrow(raw))
    }
    else {
        cv <- ukb_array(h, "p41270")
        stopifnot(identical(as.integer(sub(".*_a", "", cv)), indices))
        code <- ukb_read(code_path, c("eid", cv))
        names(code)[match(cv, names(code))] <- paste0("i", seq_along(cv))
    }
    stopifnot(setequal(code$eid, date$eid))
    date <- date[match(code$eid, date$eid), , drop = FALSE]
    for (i in seq_along(dv)) {
        code[[paste0("i", i)]][code[[paste0("i", i)]] == ""] <- NA_character_
        date[[dv[i]]] <- data.table::as.IDate(date[[dv[i]]])
    }
    names(date)[match(dv, names(date))] <- paste0("d", seq_along(dv))
    list(icd10.code0 = code, icd10.date0 = date)
}


# 5. Assemble participant tables
make_ukb_all <- function(ukb_root = ukb, tier = c("rdata", "tsv"), geography = c("recompute", "extract"), raw_tsv = NULL, 
    vip_tsv = file.path(ukb_root, "rap/vip.tab.gz"), pca_tsv = file.path(ukb_root, "rap/pca.tab.gz"), death_tsv = file.path(ukb_root, 
        "rap/death.tab.gz"), code_tsv = file.path(ukb_root, "rap/icd10.tab.gz"), date_tsv = file.path(ukb_root, 
        "rap/icd10Date.tab.gz")) {
    tier <- match.arg(tier)
    geography <- match.arg(geography)
    if (!requireNamespace("data.table", quietly = TRUE)) 
        stop("data.table is required")
    if (geography == "recompute" && !requireNamespace("RANN", quietly = TRUE)) 
        stop("RANN is required for recompute")
    sources <- character()
    record <- function(f) {
        stopifnot(file.exists(f))
        sources <<- unique(c(sources, normalizePath(f)))
        f
    }
    rds <- function(n) {
        file <- file.path(ukb_root, paste0(n, ".rds"))
        if (!file.exists(file)) 
            file <- file.path(ukb_root, "Rdata", paste0(n, ".rds"))
        ukb_check(as.data.frame(readRDS(record(file))))
    }
    if (tier == "rdata") {
        p <- rds("phe0")[, c("eid", ukb_fields())]
        pc <- rds("pca")[, c("eid", paste0("PC", 1:40))]
        death <- rds("death")[, c("eid", "date_death", "death_icd10", "date_cvd_death")]
        icd <- list(icd10.code0 = rds("icd10.code0"), icd10.date0 = rds("icd10.date0"))
    }
    else {
        if (!is.null(raw_tsv)) {
            h <- ukb_header(raw_tsv)
            fields <- c(names(ukb_fields(TRUE)), "p22009", "p40000", "p40001", "p41270", "p41280")
            cols <- names(h)[grepl(paste0("^(", paste(fields, collapse = "|"), ")(_|$)"), names(h))]
            raw <- ukb_read(record(raw_tsv), c("eid", cols))
            vip_tsv <- pca_tsv <- death_tsv <- code_tsv <- date_tsv <- raw
        }
        input <- function(f) if (is.data.frame(f)) 
            f
        else record(f)
        p <- ukb_vip(input(vip_tsv))
        pc <- ukb_pc(input(pca_tsv))
        death <- ukb_death(input(death_tsv))
        icd <- ukb_icd(input(code_tsv), input(date_tsv))
        if (!is.null(raw_tsv)) {
            rm(raw, vip_tsv, pca_tsv, death_tsv, code_tsv, date_tsv)
            invisible(gc())
        }
    }
    p <- ukb_derive(p)
    if (geography == "extract") {
        geo <- rds("move")[, c("eid", "move", "tdi.birth", "tdi.birth.nn")]
        stopifnot(setequal(geo$eid, p$eid[!is.na(p$home_east) & !is.na(p$birth_east)]))
        md <- ukb_md(p[match(geo$eid, p$eid), ])
        stopifnot(all(is.finite(md)), all(is.finite(geo$move)), all(abs(md - geo$move) < 1e-08))
    }
    else geo <- ukb_geo(p)
    x <- Reduce(function(a, b) merge(a, b, by = "eid", all = TRUE), list(p, pc, death, geo))
    x <- x[!grepl("^-", x$eid), , drop = FALSE]
    x <- x[, setdiff(names(x), c("edu", "ethnicity", "smoke_status", "alcohol_status"))]
    stopifnot(identical(icd$icd10.code0$eid, icd$icd10.date0$eid), identical(names(icd$icd10.code0)[-1], paste0("i", 
        seq_len(ncol(icd$icd10.code0) - 1L))), identical(names(icd$icd10.date0)[-1], paste0("d", seq_len(ncol(icd$icd10.code0) - 
        1L))))
    provenance <- list(tier = tier, geography = geography, sources = sources, field_map = ukb_fields(TRUE), RANN = if (geography == 
        "recompute") as.character(utils::packageVersion("RANN")) else NA_character_, recipe = "Baseline field extraction and fixed phenotype definitions", 
        notes = c("Education=max qualification code, then years lookup; OSGB Euclidean km; no lon/lat.", "TDI: source order, east-eligible population, first duplicate home, RANN k=1; cached TDI inherited in extract mode.", 
            "Relevant-field/ID projection, not the legacy 15-table union; no all.rds input.", "ICD blank/NA tokens become true missing values; legacy code caches may retain literal NA strings.", 
            "Death cutoff 2023-04-01 in TSV tier; Rdata tier inherits extract QC. Job history/wellbeing/omics remain separate."))
    attr(x, "ukb_provenance") <- provenance
    c(list(all = x, death = death), icd, list(provenance = provenance))
}

save_ukb_all <- function(bundle, phe_dir, overwrite = FALSE) {
    files <- c(all = file.path(phe_dir, "all.rds"), setNames(file.path(phe_dir, "ukb", paste0(c("death", "icd10.code0", 
        "icd10.date0"), ".rds")), c("death", "icd10.code0", "icd10.date0")))
    stopifnot(all(names(files) %in% names(bundle)))
    if (!overwrite && any(file.exists(files))) 
        stop("Output exists; choose a new directory")
    if (any(normalizePath(files, mustWork = FALSE) %in% bundle$provenance$sources)) 
        stop("Cannot overwrite source inputs")
    dir.create(file.path(phe_dir, "ukb"), recursive = TRUE, showWarnings = FALSE)
    for (n in names(files)) saveRDS(bundle[[n]], files[[n]], compress = "gzip")
    invisible(files)
}

prepare_ukb <- function() {
    raw <- Sys.getenv("MOVE_UKB_RAW", file.path(ukb, "pheno.tsv"))
    if (!file.exists(raw) && !nzchar(Sys.getenv("MOVE_UKB_RAW"))) 
        raw <- ""
    b <- make_ukb_all(ukb_root = ukb, tier = if (nzchar(raw)) 
        "tsv"
    else "rdata", raw_tsv = if (nzchar(raw)) 
        raw
    else NULL, geography = Sys.getenv("MOVE_GEOGRAPHY", "recompute"))
    save_ukb_all(b, phe, overwrite = TRUE)
}

tdi_group <- function(x) {
    c0 <- median(x$rtdi[is.finite(x$rtdi)], na.rm = TRUE)
    x[, `:=`(tdi4, fcase(btdi <= c0 & rtdi <= c0, "Low to low", btdi > c0 & rtdi <= c0, "High to low", btdi <= 
        c0 & rtdi > c0, "Low to high", btdi > c0 & rtdi > c0, "High to high"))]
    data.table::setattr(x, "tdi_cut", c0)
    x
}

ses_pca <- function(x) {
    v <- c("edu.years", "income", "camsis")
    ok <- complete.cases(x[, ..v])
    p <- prcomp(data.matrix(x[ok, ..v]), center = TRUE, scale. = TRUE)
    if (sum(p$rotation[, 1]) < 0) {
        p$rotation[, 1] <- -p$rotation[, 1]
        p$x[, 1] <- -p$x[, 1]
    }
    x[, `:=`(ses, NA_real_)]
    x[ok, `:=`(ses, as.numeric(p$x[, 1]))]
    data.table::setattr(x, "ses_info", data.table(variable = v, loading = p$rotation[, 1], variance = summary(p)$importance[2, 
        1], n = sum(ok)))
    x
}

job_status <- function() {
    x <- fread(file.path(code, "reference", "occupation.tsv"))
    setnames(x, "soc2000", "job_code")
    x[, `:=`(job_code, as.integer(job_code))]
    x
}


# 6. Prepare SES and analysis covariates
make_data <- function() {
    x <- as.data.table(readRDS(file.path(phe, "all.rds")))
    x[, `:=`(eid, as.integer(eid))]
    x[, `:=`(md, sqrt((as.numeric(home_east) - as.numeric(birth_east))^2 + (as.numeric(home_north) - as.numeric(birth_north))^2)/1000)]
    x[, `:=`(btdi = as.numeric(tdi.birth), rtdi = as.numeric(tdi), dtdi = as.numeric(tdi) - as.numeric(tdi.birth))]
    x <- tdi_group(x)
    keep <- unique(c("eid", "md", "dtdi", "btdi", "rtdi", "tdi4", "ethnic.c", "age", "sex", "center", "gen_batch", 
        paste0("PC", 1:5), "bmi", "smoke.c", "alcohol.c", "date_attend", "date_lost", "date_death", "death_icd10", 
        "inc", "edu.b", "edu.years", "edu.sco", "job_code", "job_num", "job_week_hours", "job_stand", "job_heavy", 
        "job_shift", "job_night", "home_east", "home_north", "birth_east", "birth_north", "lon.birth", "lat.birth", 
        "lon.home", "lat.home"))
    keep <- keep[keep %in% names(x)]
    d <- x[ethnic.c == "White" & is.finite(md), ..keep]
    d[, `:=`(job_code, as.integer(job_code))]
    d <- merge(d, job_status()[, .(job_code, occupation, siops, isei, egp, mcamsis, fcamsis)], by = "job_code", 
        all.x = TRUE, sort = FALSE)
    d[, `:=`(camsis, fifelse(sex == 1, mcamsis, fifelse(sex == 0, fcamsis, NA_real_)))]
    d[, `:=`(c("mcamsis", "fcamsis"), NULL)]
    setnames(d, c("inc", "smoke.c", "alcohol.c", "job_week_hours", "job_stand", "job_heavy", "job_shift", "job_night"), 
        c("income", "smoke", "alcohol", "work_hours", "work_stand", "work_heavy", "work_shift", "work_night"), 
        skip_absent = TRUE)
    d[, `:=`(college, as.integer(as.character(edu.b)))]
    d[, `:=`(sex = factor(sex), center = factor(center), gen_batch = factor(paste0("b", gen_batch)), smoke = factor(smoke), 
        alcohol = factor(alcohol), tdi4 = factor(tdi4, levels = tdi_lev), date_attend = as.IDate(date_attend), 
        date_lost = as.IDate(date_lost), date_death = as.IDate(date_death), death_icd10 = fmt(death_icd10))]
    v <- setdiff(cov1, "btdi_z")
    d <- d[complete.cases(d[, ..v])]
    setorder(d, eid)
    d <- tdi_group(d)
    d[, `:=`(migration_distance = md, delta_tdi = dtdi, tdi_birth = btdi, tdi_recruitment = rtdi, migration_distance_rint = rint(md), 
        delta_tdi_rint = rint(dtdi))]
    num <- intersect(c("md", "dtdi", "btdi", "rtdi", "income", "college", "edu.years", "edu.sco", "camsis", "siops", 
        "isei", "work_hours", "work_stand", "work_heavy", "work_shift", "job_num", "job_n", "job_change", "job_distinct", 
        "camsis_first", "camsis_last", "camsis_mean", "camsis_delta", "siops_mean", "isei_mean", "shift_n", "shift_change", 
        "shift_ever"), names(d))
    d[, `:=`(paste0(num, "_rint"), lapply(.SD, rint)), .SDcols = num]
    saveRDS(d, file.path(phe, "move.rds"), compress = "gzip")
    message("N = ", nrow(d))
    d
}

norm_dat <- function(x) {
    if (!"md" %in% names(x) && "migration_distance" %in% names(x)) 
        x[, `:=`(md, migration_distance)]
    if (!"dtdi" %in% names(x) && "delta_tdi" %in% names(x)) 
        x[, `:=`(dtdi, delta_tdi)]
    if (!"btdi" %in% names(x) && "tdi_birth" %in% names(x)) 
        x[, `:=`(btdi, tdi_birth)]
    if (!"rtdi" %in% names(x) && "tdi_recruitment" %in% names(x)) 
        x[, `:=`(rtdi, tdi_recruitment)]
    x <- tdi_group(x)
    x[, `:=`(tdi4, factor(tdi4, levels = tdi_lev))]
    x[]
}

urban_label <- function(from, to) {
    fifelse(from %in% 0:1 & to %in% 0:1, paste(ifelse(from == 1, "Urban", "Rural"), "to", ifelse(to == 1, "urban", 
        "rural")), NA_character_)
}

urban_history <- function(x, id) {
    setorderv(x, c(id, "wave"))
    x[, `:=`(urban_before, shift(urban)), by = id]
    x[, `:=`(urban_origin, {
        i <- which(!is.na(urban))[1]
        if (length(i)) 
            urban[i]
        else NA_integer_
    }), by = id]
    x[, `:=`(urban_origin_wave, {
        i <- which(!is.na(urban))[1]
        if (length(i)) 
            wave[i]
        else NA_integer_
    }), by = id]
    x[, `:=`(urban_transition, fifelse(moved == 1 & !is.na(urban_before) & !is.na(urban), urban_label(urban_before, 
        urban), NA_character_))]
    x[, `:=`(urban_moves, cumsum(fifelse(!is.na(moved) & moved == 1 & wave > urban_origin_wave, 1L, 0L))), by = id]
    x[, `:=`(urban_route, fifelse(urban_moves > 0 & !is.na(urban), urban_label(urban_origin, urban), NA_character_))]
    x[]
}

dat <- function() {
    f <- file.path(phe, "move.rds")
    if (!file.exists(f)) 
        make_data()
    else norm_dat(as.data.table(readRDS(f)))
}

base <- function() {
    x <- dat()
    tdi_cut <- attr(x, "tdi_cut")
    x[, `:=`(eid, as.character(eid))]
    if (testn > 0) 
        x <- x[seq_len(min(.N, testn))]
    death <- as.data.table(readRDS(ukb_input("death")))
    death[, `:=`(eid, as.character(eid))]
    death[, `:=`(date_death2 = as.IDate(date_death), date_cvd = as.IDate(date_cvd_death))]
    x <- merge(x, death[, .(eid, date_death2, date_cvd)], by = "eid", all.x = TRUE, sort = FALSE)
    x[is.na(date_death), `:=`(date_death, date_death2)]
    x[, `:=`(date_death2, NULL)]
    x <- ses_pca(x)
    ses_info <- attr(x, "ses_info")
    x[, `:=`(md_z = z(md), dtdi_z = z(dtdi), btdi_z = z(btdi), ses_z = z(ses), camsis_z = z(camsis), age_z = z(age), 
        md_q = cut(md, unique(quantile(md, seq(0, 1, 0.25), na.rm = TRUE)), include.lowest = TRUE, labels = paste0("Q", 
            1:4)), sex2 = factor(fifelse(as.character(sex) %in% c("1", "Male", "M"), "Male", "Female"), levels = c("Female", 
            "Male")))]
    x[, `:=`(tdi4, factor(as.character(tdi4), levels = tdi_lev))]
    x[, `:=`(sex = factor(sex), center = factor(center), gen_batch = factor(gen_batch), smoke = factor(smoke), 
        alcohol = factor(alcohol))]
    x <- x[complete.cases(x[, c("md_z", "dtdi_z", cov1, "date_attend"), with = FALSE])]
    data.table::setattr(x, "tdi_cut", tdi_cut)
    data.table::setattr(x, "ses_info", ses_info)
    x
}

if (sys.nframe() == 0L) {
    prepare_ukb()
    make_data()
    saveRDS(list(all = nrow(readRDS(file.path(phe, "all.rds"))), eligible = nrow(dat())), file.path(out, "sample.rds"))
}

