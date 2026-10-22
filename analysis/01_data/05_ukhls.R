# 1. UKHLS wave harmonization
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

ukhls_na <- function(x) {
    x <- suppressWarnings(as.numeric(haven::zap_labels(x)))
    x[x %in% c(-21, -20, -11, -10, -9, -8, -7, -2, -1)] <- NA_real_
    x
}

ukhls_id <- function(x) {
    y <- rep(NA_character_, length(x))
    i <- is.finite(x)
    y[i] <- format(x[i], scientific = FALSE, trim = TRUE)
    y
}

ukhls_get <- function(x, v) {
    if (v %in% names(x)) 
        as.numeric(x[[v]])
    else rep(NA_real_, nrow(x))
}

ukhls_pick <- function(x, v) {
    v <- intersect(unique(v), names(x))
    if (!length(v)) 
        return(rep(NA_real_, nrow(x)))
    Reduce(data.table::fcoalesce, lapply(v, function(z) as.numeric(x[[z]])))
}

ukhls_weight <- function(x, v) {
    v <- intersect(unique(v), names(x))
    if (!length(v)) 
        return(rep(NA_real_, nrow(x)))
    z <- lapply(v, function(y) {
        a <- as.numeric(x[[y]])
        a[a <= 0] <- NA_real_
        a
    })
    Reduce(data.table::fcoalesce, z)
}


# 2. Read wave-specific fields
ukhls_read <- function(zip, member, wave, keep, tmp, health = FALSE) {
    unzip(zip, files = member, exdir = tmp, junkpaths = TRUE)
    file <- file.path(tmp, basename(member))
    head <- read_dta(file, n_max = 1, encoding = "latin1")
    src <- names(head)
    var <- if (nzchar(wave)) 
        sub(paste0("^", wave, "_"), "", src)
    else src
    take <- var %in% keep | grepl("^ind.*_(lw|xw)$", var) | (health & grepl("^(hcond|mhcond)", var))
    src <- src[take]
    var <- var[take]
    meta <- data.table(wave = wave, file = basename(member), source = src, variable = var, label = vapply(src, 
        function(z) {
            y <- attr(head[[z]], "label")
            if (is.null(y)) 
                ""
            else y
        }, character(1)), codes = vapply(src, function(z) {
        y <- attr(head[[z]], "labels")
        if (is.null(y)) 
            ""
        else paste(names(y), unname(y), sep = "=", collapse = ";")
    }, character(1)))
    x <- as.data.table(read_dta(file, col_select = tidyselect::all_of(src), encoding = "latin1"))
    setnames(x, src, var)
    for (z in names(x)) if (is.numeric(x[[z]]) || inherits(x[[z]], "haven_labelled")) 
        set(x, j = z, value = ukhls_na(x[[z]]))
    unlink(file)
    list(data = x, meta = meta)
}

ukhls_flag <- function(x, meta, mode, pattern) {
    m <- meta[grepl(pattern, label, ignore.case = TRUE)]
    if (mode == "inc") {
        m <- m[(grepl("^(hcondn|mhcondn)", variable) | grepl("new health conditions?", label, ignore.case = TRUE)) & 
            !grepl("^(hcondns|mhcondns)|still has", paste(variable, label), ignore.case = TRUE)]
    }
    else {
        m <- m[!grepl("^(hcondn|mhcondn)", variable) & grepl("^(diagnosed health conditions?|health condition:|type of emotional)", 
            label, ignore.case = TRUE)]
    }
    v <- intersect(unique(m$variable), names(x))
    if (!length(v)) 
        return(rep(NA_integer_, nrow(x)))
    z <- as.matrix(x[, ..v])
    ok <- rowSums(!is.na(z)) > 0
    y <- rep(NA_integer_, nrow(x))
    y[ok] <- as.integer(rowSums(z[ok, , drop = FALSE] == 1, na.rm = TRUE) > 0)
    y
}


# 3. Harmonize exposures and outcomes
ukhls_wave <- function(zip, wave, tmp, disease) {
    dir <- "UKDA-6614-stata/stata/stata14_se/ukhls"
    ind <- c("pidp", "hidp", "intdatd_dv", "intdatm_dv", "intdaty_dv", "addrmov_dv", "distmov_dv", "adcts", "family", 
        "education", "memploy", "housing", "area", "movdir", "plnowm", "plnowy4", "origadd", "age_dv", "sex_dv", 
        "ethn_dv", "mastat_dv", "sf1", "health", "sf12pcs_dv", "sf12mcs_dv", "scghq1_dv", "scghq2_dv", "bmi_dv", 
        "smoker", "ncigs", "sceverdrnk", "evralc", "dklm", "fivealcdr", "fimngrs_dv", "fimnnet_dv", "fimnlabgrs_dv", 
        "qfhigh_dv", "hiqual_dv", "jbstat", "jbsoc00_cc", "jbsoc10_cc", "jbsoc20_cc", "jbnssec8_dv", "jbhrs", "jbot", 
        "wkphys")
    hh <- c("hidp", "fihhmngrs_dv", "fihhmnnet1_dv", "gor_dv", "urban_dv", "tenure_dv", "hhsize", "nchild_dv", 
        "nadults_dv")
    a <- ukhls_read(zip, sprintf("%s/%s_indresp.dta", dir, wave), wave, ind, tmp, TRUE)
    b <- ukhls_read(zip, sprintf("%s/%s_hhresp.dta", dir, wave), wave, hh, tmp)
    x <- merge(a$data, b$data, by = "hidp", all.x = TRUE, sort = FALSE)
    meta <- rbindlist(list(a$meta, b$meta), fill = TRUE)
    mv <- ukhls_get(x, "addrmov_dv")
    moved <- fcase(mv == 1, 1L, mv == 2, 0L, default = NA_integer_)
    dist <- ukhls_get(x, "distmov_dv")
    date <- suppressWarnings(as.IDate(sprintf("%04d-%02d-%02d", as.integer(ukhls_get(x, "intdaty_dv")), as.integer(ukhls_get(x, 
        "intdatm_dv")), as.integer(ukhls_get(x, "intdatd_dv")))))
    lw <- grep("^ind.*_lw$", names(x), value = TRUE)
    xw <- grep("^ind.*_xw$", names(x), value = TRUE)
    d <- data.table(pid = ukhls_id(ukhls_get(x, "pidp")), hid = ukhls_id(ukhls_get(x, "hidp")), wave = match(wave, 
        letters), date = date, moved = moved, md = fcase(moved == 0, 0, moved == 1, dist, default = NA_real_), 
        md_move = fifelse(moved == 1, dist, NA_real_), same_address = ukhls_get(x, "adcts"), original_address = ukhls_get(x, 
            "origadd"), move_direct = ukhls_get(x, "movdir"), reason_family = ukhls_get(x, "family"), reason_education = ukhls_get(x, 
            "education"), reason_employment = ukhls_get(x, "memploy"), reason_housing = ukhls_get(x, "housing"), 
        reason_area = ukhls_get(x, "area"), residence_month = ukhls_get(x, "plnowm"), residence_year = ukhls_get(x, 
            "plnowy4"), age = ukhls_get(x, "age_dv"), sex_code = ukhls_get(x, "sex_dv"), ethnicity = ukhls_get(x, 
            "ethn_dv"), marital = ukhls_get(x, "mastat_dv"), region = ukhls_get(x, "gor_dv"), urban_code = ukhls_get(x, 
            "urban_dv"), tenure = ukhls_get(x, "tenure_dv"), hh_size = ukhls_get(x, "hhsize"), children = ukhls_get(x, 
            "nchild_dv"), adults = ukhls_get(x, "nadults_dv"), qualification = ukhls_pick(x, c("hiqual_dv", "qfhigh_dv")), 
        employment = ukhls_get(x, "jbstat"), occupation = ukhls_pick(x, c("jbsoc10_cc", "jbsoc20_cc", "jbsoc00_cc")), 
        nssec = ukhls_get(x, "jbnssec8_dv"), work_hours = ukhls_get(x, "jbhrs"), overtime = ukhls_get(x, "jbot"), 
        job_physical = ukhls_get(x, "wkphys"), income = ukhls_get(x, "fimnnet_dv"), income_gross = ukhls_get(x, 
            "fimngrs_dv"), hh_income = ukhls_get(x, "fihhmnnet1_dv"), hh_income_gross = ukhls_get(x, "fihhmngrs_dv"), 
        health_self = ukhls_get(x, "sf1"), illness = ukhls_get(x, "health"), physical = ukhls_get(x, "sf12pcs_dv"), 
        mental = ukhls_get(x, "sf12mcs_dv"), ghq = ukhls_get(x, "scghq1_dv"), distress = ukhls_get(x, "scghq2_dv"), 
        bmi = ukhls_get(x, "bmi_dv"), smoker = ukhls_get(x, "smoker"), cigarettes = ukhls_get(x, "ncigs"), alcohol_ever = ukhls_pick(x, 
            c("sceverdrnk", "evralc")), alcohol_month = ukhls_get(x, "dklm"), alcohol_heavy = ukhls_get(x, "fivealcdr"), 
        weight_lw = ukhls_weight(x, c("inding2_lw", "indinui_lw", "indinub_lw", "indinus_lw", lw)), weight_xw = ukhls_weight(x, 
            c("inding2_xw", "indinui_xw", "indinub_xw", "indinus_xw", xw)))
    for (z in names(disease)) {
        d[, `:=`((paste0("inc_", z)), ukhls_flag(x, meta, "inc", disease[[z]]))]
        d[, `:=`((paste0("dx_", z)), ukhls_flag(x, meta, "dx", disease[[z]]))]
    }
    list(data = d, meta = meta)
}


# 4. Combine longitudinal records
make_ukhls <- function(waves = letters[1:15], save = TRUE) {
    zip <- Sys.glob(file.path(phe, "ukhls", "6614stata_*.zip"))
    if (length(zip) != 1) 
        stop("one study 6614 Stata archive is required")
    if (!all(waves %in% letters[1:15])) 
        stop("waves must be a-o")
    tmp <- tempfile("ukhls_")
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    dir <- "UKDA-6614-stata/stata/stata14_se/ukhls"
    disease <- c(cad_mi = "coronary heart disease|angina|heart attack|myocardial infarction", stroke = "stroke", 
        heart_failure = "heart failure", atrial_fibrillation = "atrial fibrillation", diabetes = "diabetes", depression = "depress", 
        anxiety = "anxiety", dementia = "dementia|alzheimer", hypertension = "high blood pressure|hypertension", 
        cancer = "cancer|malignancy", asthma = "asthma", arthritis = "arthritis")
    cross <- ukhls_read(zip, file.path(dir, "xwavedat.dta"), "", c("pidp", "hhorig", "memorig", "sex_dv", "birthy", 
        "racel_dv", "ethn_dv", "dcsedfl_dv", "dcsedw_dv"), tmp)
    death <- ukhls_get(cross$data, "dcsedfl_dv")
    dw <- ukhls_get(cross$data, "dcsedw_dv")
    dw <- fifelse(!is.na(dw) & dw >= 19 & dw <= 33, dw - 18, NA_real_)
    person0 <- data.table(pid = ukhls_id(ukhls_get(cross$data, "pidp")), sample = ukhls_pick(cross$data, c("memorig", 
        "hhorig")), sex_cross = ukhls_get(cross$data, "sex_dv"), birth_year = ukhls_get(cross$data, "birthy"), 
        ethnic_cross = ukhls_pick(cross$data, c("racel_dv", "ethn_dv")), death = fcase(death == 1, 1L, death == 
            3, 0L, default = NA_integer_), death_wave = as.integer(dw))
    w <- lapply(waves, function(z) {
        message("UKHLS wave ", toupper(z))
        ukhls_wave(zip, z, tmp, disease)
    })
    panel <- rbindlist(lapply(w, `[[`, "data"), fill = TRUE)
    meta <- rbindlist(c(list(cross$meta), lapply(w, `[[`, "meta")), fill = TRUE)
    panel <- merge(panel, person0, by = "pid", all.x = TRUE, sort = FALSE)
    panel[is.na(sex_code), `:=`(sex_code, sex_cross)]
    panel[is.na(ethnicity), `:=`(ethnicity, ethnic_cross)]
    panel[, `:=`(c("sex_cross", "ethnic_cross"), NULL)]
    panel[, `:=`(sex = fcase(sex_code == 1, "Male", sex_code == 2, "Female"), white = fifelse(is.na(ethnicity), 
        NA, ethnicity %in% 1:4), urban = fcase(urban_code == 1, 1L, urban_code == 2, 0L, default = NA_integer_), 
        degree = fifelse(is.na(qualification), NA_integer_, as.integer(qualification == 1)), employed = fifelse(is.na(employment), 
            NA_integer_, as.integer(employment %in% c(1, 2, 10, 11))), smoke = fcase(smoker == 1, 1L, smoker == 
            2, 0L, default = NA_integer_), log_income = fifelse(income > 0, log(pmax(income, .Machine$double.eps)), 
            NA_real_), log_hh_income = fifelse(hh_income > 0, log(pmax(hh_income, .Machine$double.eps)), NA_real_))]
    setorder(panel, pid, wave)
    locf <- c("bmi", "smoke", "cigarettes", "alcohol_ever", "alcohol_month", "alcohol_heavy", "job_physical")
    panel[, `:=`(paste0(locf, "_locf"), lapply(.SD, data.table::nafill, type = "locf")), by = pid, .SDcols = locf]
    for (z in names(disease)) {
        dx <- paste0("dx_", z)
        inc <- paste0("inc_", z)
        ever <- paste0("disease_", z)
        panel[, `:=`((ever), {
            a <- get(dx)
            b <- get(inc)
            v <- fifelse(a %in% 1L | b %in% 1L, 1L, fifelse(!is.na(a) | !is.na(b), 0L, NA_integer_))
            seen <- cumsum(!is.na(v)) > 0
            y <- cummax(fifelse(is.na(v), 0L, v))
            y[!seen] <- NA_integer_
            y
        }), by = pid]
    }
    panel[, `:=`(move_count = cumsum(fifelse(moved == 1, 1L, 0L)), distance_total = cumsum(fifelse(is.finite(md_move), 
        md_move, 0))), by = pid]
    urban_history(panel, "pid")
    person <- panel[, .(first_wave = min(wave), last_wave = max(wave), waves = .N, first_date = if (all(is.na(date))) 
        as.IDate(NA)
    else min(date, na.rm = TRUE), last_date = if (all(is.na(date))) 
        as.IDate(NA)
    else max(date, na.rm = TRUE), move_observed = sum(!is.na(moved)), moves = sum(moved == 1, na.rm = TRUE), distance_total = if (any(moved == 
        1, na.rm = TRUE) && !any(is.finite(md_move))) 
        NA_real_
    else sum(md_move, na.rm = TRUE), distance_mean = if (any(is.finite(md_move))) 
        mean(md_move, na.rm = TRUE)
    else NA_real_), by = pid]
    person[, `:=`(follow_years, as.numeric(last_date - first_date)/365.25)]
    person <- merge(person, person0, by = "pid", all.x = TRUE, sort = FALSE)
    episode <- copy(panel)
    pre <- intersect(c("age", "sex_code", "ethnicity", "white", "region", "urban", "marital", "qualification", 
        "degree", "employment", "employed", "occupation", "nssec", "work_hours", "job_physical", "income", "income_gross", 
        "hh_income", "hh_income_gross", "tenure", "hh_size", "children", "adults", "health_self", "illness", "physical", 
        "mental", "ghq", "distress", "bmi", "smoke", "cigarettes", "alcohol_ever", "alcohol_month", "alcohol_heavy", 
        paste0(locf, "_locf"), "log_income", "log_hh_income", "weight_lw", paste0("disease_", names(disease))), 
        names(episode))
    post <- intersect(c("health_self", "illness", "physical", "mental", "ghq", "distress", "bmi", "income", "hh_income", 
        "log_income", "log_hh_income", "employment", "nssec", paste0("inc_", names(disease)), paste0("disease_", 
            names(disease))), names(episode))
    episode[, `:=`(wave0 = shift(wave), date0 = shift(date), wave2 = shift(wave, type = "lead"), date2 = shift(date, 
        type = "lead")), by = pid]
    episode[, `:=`((paste0(pre, "0")), shift(.SD)), by = pid, .SDcols = pre]
    episode[, `:=`((paste0(post, "2")), shift(.SD, type = "lead")), by = pid, .SDcols = post]
    episode <- episode[!is.na(wave0) & !is.na(moved)]
    episode[, `:=`(wave_gap = wave - wave0, follow_years = as.numeric(date - date0)/365.25, next_gap = wave2 - 
        wave, next_years = as.numeric(date2 - date)/365.25)]
    change <- intersect(c("health_self", "physical", "mental", "ghq", "bmi", "income", "hh_income", "log_income", 
        "log_hh_income"), pre)
    for (z in change) episode[, `:=`((paste0("change_", z)), get(z) - get(paste0(z, "0")))]
    episode <- merge(episode, person[, .(pid, last_wave)], by = "pid", all.x = TRUE, sort = FALSE)
    episode[, `:=`(death_event, fcase(death == 1 & death_wave > wave, 1L, !is.na(death), 0L, default = NA_integer_))]
    episode[, `:=`(death_follow_waves, fifelse(death_event == 1, death_wave - wave, last_wave - wave))]
    wave_qc <- panel[, .(persons = uniqueN(pid), moved_known = sum(!is.na(moved)), movers = sum(moved == 1, na.rm = TRUE), 
        distance_known = sum(is.finite(md)), health_known = sum(is.finite(health_self)), physical_known = sum(is.finite(physical)), 
        income_known = sum(is.finite(income))), by = wave]
    qc <- data.table(item = c("persons", "person-waves", "migration intervals", "movers", "known migration distance", 
        "next-wave health follow-up", "deaths", "waves"), value = c(uniqueN(panel$pid), nrow(panel), nrow(episode), 
        sum(episode$moved == 1), sum(is.finite(episode$md)), sum(is.finite(episode$health_self2)), sum(person$death == 
            1, na.rm = TRUE), uniqueN(panel$wave)))
    design <- data.table(target = c("Migration distance", "Delta TDI", "Disease", "Mortality", "SES", "Survey weight"), 
        definition = c("Distance moved since the previous interview; non-movers are 0 and md_move is movers only", 
            "Unavailable in study 6614 without restricted geography or linked area deprivation", "Self-reported diagnosis since the previous interview plus cumulative diagnosed status", 
            "Death reported by survey wave; exact date and cause are unavailable", "Repeated education, income, employment, occupation, NS-SEC and housing measures", 
            "Broadest applicable positive longitudinal and cross-sectional adult interview weights"), replication = c("Longitudinal analogue", 
            "Not replicated", "Partial", "Partial", "Available", "Available"))
    ans <- list(panel = panel, episode = episode, person = person, meta = list(variable = unique(meta), disease = data.table(outcome = names(disease), 
        pattern = unname(disease)), design = design), qc = list(summary = qc, wave = wave_qc))
    if (save) 
        saveRDS(ans, file.path(phe, "ukhls.rds"), compress = "gzip")
    message("UKHLS: ", uniqueN(panel$pid), " persons; ", nrow(episode), " intervals; ", sum(episode$moved == 1), 
        " moves")
    invisible(ans)
}

if (sys.nframe() == 0L) {
    make_ukhls()
}

