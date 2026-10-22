# 1. HRS wave harmonization
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

hrs_require <- function() {
    packages <- c("data.table", "haven", "tidyselect")
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) 
        stop("Missing HRS packages: ", paste(missing, collapse = ", "), ". Install locally before building.", call. = FALSE)
    invisible(TRUE)
}


# 2. Locate raw releases
hrs_sources <- function(root = Sys.getenv("MOVE_ROOT", getwd())) {
    data.frame(source = c("move", "rand", "track"), archive = file.path(root, "data/phe/hrs", c("HRSXRegion2022v10Early.zip", 
        "randhrs1992_2022v1_STATA.zip", "trk2022v1.zip")), member = c("HRSXREGION22.dta", "randhrs1992_2022v1.dta", 
        "trk2022tr_r.dta"), stringsAsFactors = FALSE)
}

hrs_audit <- function(root = Sys.getenv("MOVE_ROOT", getwd())) {
    a <- hrs_sources(root)
    a$available <- file.exists(a$archive)
    a$archive_bytes <- file.info(a$archive)$size
    a$dta_bytes <- NA_real_
    a$status <- "missing archive"
    for (i in seq_len(nrow(a))) {
        if (!a$available[i]) 
            next
        z <- tryCatch(utils::unzip(a$archive[i], list = TRUE), error = function(e) NULL)
        if (is.null(z)) {
            a$status[i] <- "unreadable archive"
            next
        }
        hit <- which(z$Name == a$member[i])
        if (length(hit) != 1L) {
            a$status[i] <- "expected DTA member missing or duplicated"
            next
        }
        a$dta_bytes[i] <- z$Length[hit]
        a$status[i] <- if (a$dta_bytes[i] > 0) 
            "ready"
        else "empty DTA member"
    }
    a
}

hrs_num <- function(x) suppressWarnings(as.numeric(haven::zap_missing(x)))

hrs_bin <- function(x) {
    x <- hrs_num(x)
    data.table::fifelse(x == 1, 1L, data.table::fifelse(x == 0, 0L, NA_integer_))
}

hrs_date <- function(x) {
    if (inherits(x, "Date")) 
        return(data.table::as.IDate(x))
    data.table::as.IDate(as.Date(hrs_num(x), origin = "1960-01-01"))
}

hrs_id <- function(hhid, pn = NULL) {
    h <- hrs_num(hhid)
    if (is.null(pn)) {
        if (any(!is.finite(h) | h < 0 | h > 999999999 | h != floor(h))) 
            stop("Invalid RAND identifiers; values suppressed.", call. = FALSE)
        return(sprintf("%09.0f", h))
    }
    p <- hrs_num(pn)
    if (length(h) != length(p) || any(!is.finite(h) | h < 0 | h > 999999 | h != floor(h)) || any(!is.finite(p) | 
        p < 0 | p > 999 | p != floor(p))) 
        stop("Invalid HRS identifiers; values suppressed.", call. = FALSE)
    paste0(sprintf("%06.0f", h), sprintf("%03.0f", p))
}

hrs_take <- function(x, v, fun = hrs_num) {
    if (!v %in% names(x)) {
        if (identical(fun, hrs_date)) 
            return(data.table::as.IDate(rep(NA_character_, nrow(x))))
        if (identical(fun, hrs_bin)) 
            return(rep(NA_integer_, nrow(x)))
        return(rep(NA_real_, nrow(x)))
    }
    fun(x[[v]])
}

hrs_read <- function(path, vars, required = vars) {
    if (!file.exists(path)) 
        stop("Missing HRS DTA file: ", path, call. = FALSE)
    x <- data.table::as.data.table(haven::read_dta(path, col_select = tidyselect::any_of(vars)))
    missing <- setdiff(required, names(x))
    if (length(missing)) 
        stop("Missing required HRS fields in ", basename(path), ": ", paste(missing, collapse = ", "), call. = FALSE)
    x
}

hrs_unique <- function(x, cols, label) {
    if (anyDuplicated(x[, cols, with = FALSE])) 
        stop("Duplicate keys in HRS ", label, "; values suppressed.", call. = FALSE)
    invisible(TRUE)
}

hrs_move_map <- function() {
    data.table::data.table(variable = c("MOVE9294", "MOVE9496", "MOVE9698", "MOVE9395", "MOVE9598", "MOVE9800", 
        "MOVE0002", "MOVE0204", "MOVE0406", "MOVE0608", "MOVE0810", "MOVE1012", "MOVE1214", "MOVE1416", "MOVE1618", 
        "MOVE1820", "MOVE2022"), start = as.integer(c(1992, 1994, 1996, 1993, 1995, 1998, 2000, 2002, 2004, 2006, 
        2008, 2010, 2012, 2014, 2016, 2018, 2020)), end = as.integer(c(1994, 1996, 1998, 1995, 1998, 2000, 2002, 
        2004, 2006, 2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022)))
}


# 3. Read migration histories
hrs_move_data <- function(path) {
    map <- hrs_move_map()
    years <- c(1992:1996, seq(1998, 2022, 2))
    beale <- paste0("BEALE2013_", sprintf("%02d", years%%100))
    v <- c("HHID", "PN", "STUDY", "FIRSTIW", "REGIONB", "REGLIV10", map$variable, beale)
    x <- hrs_read(path, v)
    x[, `:=`(uid, hrs_id(HHID, PN))]
    hrs_unique(x, "uid", "mobility raw data")
    urban <- rbindlist(Map(function(y, field) data.table(uid = x$uid, year = y, urban = fcase(hrs_num(x[[field]]) %in% 
        1:2, 1L, hrs_num(x[[field]]) == 3, 0L, default = NA_integer_)), years, beale))
    urban <- urban[!is.na(urban)]
    for (j in map$variable) data.table::set(x, j = j, value = hrs_num(x[[j]]))
    hrs <- c("MOVE9294", "MOVE9496", "MOVE9698")
    ahead <- c("MOVE9395", "MOVE9598")
    nh <- rowSums(!is.na(x[, hrs, with = FALSE]))
    na <- rowSums(!is.na(x[, ahead, with = FALSE]))
    fi <- hrs_num(x$FIRSTIW)
    route <- data.table::fcase(fi %in% c(1992, 1994, 1996), "hrs", fi %in% c(1993, 1995), "ahead", nh >= na & nh > 
        0, "hrs", na > 0, "ahead", default = "later")
    x[route == "hrs", `:=`((ahead), NA_real_)]
    x[route == "ahead", `:=`((hrs), NA_real_)]
    meta <- x[, .(uid, study = hrs_num(STUDY), first_iw = fi, birth_region = hrs_num(REGIONB), school_region = hrs_num(REGLIV10), 
        early_overlap = nh > 0 & na > 0)]
    x <- data.table::melt(x, id.vars = "uid", measure.vars = map$variable, variable.name = "variable", value.name = "miles", 
        na.rm = TRUE)
    x <- merge(x, map, by = "variable", all.x = TRUE, sort = FALSE)
    x <- x[is.finite(miles) & miles >= 0, .(uid, start, end, km = miles * 1.609344, moved = as.integer(miles > 
        0))]
    x <- merge(x, urban[, .(uid, start = year, urban_start = urban)], by = c("uid", "start"), all.x = TRUE, sort = FALSE)
    x <- merge(x, urban[, .(uid, end = year, urban_end = urban)], by = c("uid", "end"), all.x = TRUE, sort = FALSE)
    x[, `:=`(urban_transition, fifelse(moved == 1, urban_label(urban_start, urban_end), NA_character_))]
    data.table::setorder(x, uid, end)
    summary <- x[, .(move_obs = .N, mt = sum(moved), md_km = sum(km), first_move = if (any(moved == 1)) 
        as.integer(min(end[moved == 1]))
    else NA_integer_, last_move = if (any(moved == 1)) 
        as.integer(max(end[moved == 1]))
    else NA_integer_), by = uid]
    list(long = x, sum = summary, meta = meta, urban = urban)
}


# 4. Harmonize repeated measurements
hrs_wave_data <- function(path) {
    w <- 1:16
    rvars <- unlist(lapply(w, function(i) paste0("r", i, c("iwbeg", "agey_b", "bmi", "smokev", "smoken", "drink", 
        "jphys", "jlift", "jstoop", "hearte", "stroke", "diabe", "psyche", "demene", "alzhee", "cesd"))))
    hvars <- unlist(lapply(w, function(i) paste0("h", i, c("itot", "atotb"))))
    stat <- c("hhidpn", "ragender", "raracem", "rahispan", "raedyrs", "raedegrm", "raddate", "radyear")
    optional <- c("r1cesd", unlist(lapply(1:9, function(i) paste0("r", i, c("demene", "alzhee")))))
    x <- hrs_read(path, c(stat, rvars, hvars), setdiff(c(stat, rvars, hvars), optional))
    uid <- hrs_id(x$hhidpn)
    person <- data.table::data.table(uid = uid, sex = data.table::fcase(hrs_num(x$ragender) == 1, "Male", hrs_num(x$ragender) == 
        2, "Female"), race = hrs_num(x$raracem), hispanic = hrs_num(x$rahispan), education = hrs_num(x$raedyrs), 
        degree = hrs_num(x$raedegrm), death = hrs_date(x$raddate), death_year = hrs_num(x$radyear))
    hrs_unique(person, "uid", "RAND raw data")
    wave <- data.table::rbindlist(lapply(w, function(i) data.table::data.table(uid = uid, wave = i, date = hrs_take(x, 
        paste0("r", i, "iwbeg"), hrs_date), age = hrs_take(x, paste0("r", i, "agey_b")), bmi = hrs_take(x, paste0("r", 
        i, "bmi")), smoke_ever = hrs_take(x, paste0("r", i, "smokev"), hrs_bin), smoke = hrs_take(x, paste0("r", 
        i, "smoken"), hrs_bin), alcohol = hrs_take(x, paste0("r", i, "drink"), hrs_bin), income = hrs_take(x, paste0("h", 
        i, "itot")), wealth = hrs_take(x, paste0("h", i, "atotb")), physical = hrs_take(x, paste0("r", i, "jphys")), 
        lifting = hrs_take(x, paste0("r", i, "jlift")), stooping = hrs_take(x, paste0("r", i, "jstoop")), heart = hrs_take(x, 
            paste0("r", i, "hearte"), hrs_bin), stroke = hrs_take(x, paste0("r", i, "stroke"), hrs_bin), diabetes = hrs_take(x, 
            paste0("r", i, "diabe"), hrs_bin), psychiatric = hrs_take(x, paste0("r", i, "psyche"), hrs_bin), dementia = hrs_take(x, 
            paste0("r", i, "demene"), hrs_bin), alzheimer = hrs_take(x, paste0("r", i, "alzhee"), hrs_bin), cesd = hrs_take(x, 
            paste0("r", i, "cesd")))), use.names = TRUE)
    wave[, `:=`(year, as.integer(format(date, "%Y")))]
    for (j in c("physical", "lifting", "stooping")) wave[!get(j) %in% 1:4, `:=`((j), NA_real_)]
    wave[, `:=`(depression, data.table::fifelse(is.finite(cesd), as.integer(cesd >= 4), NA_integer_))]
    wave[, `:=`(dementia, data.table::fifelse(dementia == 1 | alzheimer == 1, 1L, data.table::fifelse(dementia == 
        0 & alzheimer == 0, 0L, NA_integer_)))]
    wave <- wave[!is.na(date)]
    data.table::setorder(wave, uid, date)
    person <- merge(person, wave[, .(first = min(date), last = max(date)), by = uid], by = "uid", all.x = TRUE)
    list(person = person, wave = wave)
}

hrs_disease_endpoint <- function(wave, name, var) {
    x <- wave[!is.na(date) & get(var) %in% 0:1, .(id, date, status = get(var))]
    data.table::setorder(x, id, date)
    x[, {
        first <- date[1]
        prevalent <- isTRUE(status[1] == 1)
        k <- if (!prevalent) 
            which(status == 1 & date > first)[1]
        else NA_integer_
        event <- !is.na(k)
        exit <- if (event) 
            date[k]
        else max(date)
        .(outcome = name, entry = first, exit = exit, event = as.integer(event), prevalent = as.integer(prevalent), 
            visits = .N)
    }, by = id]
}

hrs_move_before <- function(endpoint, move) {
    endpoint <- data.table::copy(endpoint)
    endpoint[, `:=`(row, .I)]
    x <- merge(endpoint[, .(row, id, entry, exit, event)], move, by = "id", all.x = TRUE, allow.cartesian = TRUE)
    x[, `:=`(entry_year = as.integer(format(entry, "%Y")), exit_year = as.integer(format(exit, "%Y")))]
    x <- x[is.finite(km) & start >= entry_year & data.table::fifelse(event == 1, end < exit_year, end <= exit_year)]
    a <- x[, .(move_obs = .N, md_km = sum(km), mt = sum(moved)), by = row]
    endpoint <- merge(endpoint, a, by = "row", all.x = TRUE, sort = FALSE)
    endpoint[, `:=`(follow_years, as.numeric(exit - entry)/365.25)]
    endpoint[, `:=`(row, NULL)]
    endpoint
}

hrs_endpoint_data <- function(person, wave, move) {
    def <- c(`Heart problems` = "heart", Stroke = "stroke", Diabetes = "diabetes", `Psychiatric problems` = "psychiatric", 
        `Elevated depressive symptoms` = "depression", Dementia = "dementia")
    out <- data.table::rbindlist(Map(function(n, v) hrs_disease_endpoint(wave, n, v), names(def), def))
    death <- person[!is.na(first), .(id, entry = first, exit = data.table::fifelse(!is.na(death) & death > first, 
        death, last), event = as.integer(!is.na(death) & death > first), prevalent = 0L, visits = NA_integer_)]
    death[, `:=`(outcome, "All-cause mortality")]
    out <- rbind(out, death, fill = TRUE)
    out <- out[prevalent == 0 & !is.na(exit) & exit > entry]
    out <- hrs_move_before(out, move)
    base <- wave[out, on = .(id, date = entry), mult = "first", .(id = i.id, outcome = i.outcome, age, bmi, smoke, 
        alcohol, income, wealth, physical, lifting, stooping)]
    out <- merge(out, base, by = c("id", "outcome"), all.x = TRUE, sort = FALSE)
    out <- merge(out, person[, .(id, sex, race, hispanic, education)], by = "id", all.x = TRUE)
    data.table::setorder(out, outcome, id)
    out
}


# 5. Assemble follow-up
hrs_assemble <- function(m, r, tracker) {
    t <- data.table::data.table(uid = hrs_id(tracker$HHID, tracker$PN), ex_month = hrs_num(tracker$EXDEATHMO), 
        ex_year = hrs_num(tracker$EXDEATHYR), ndi_month = hrs_num(tracker$NMONTH), ndi_year = hrs_num(tracker$NYEAR))
    hrs_unique(t, "uid", "tracker raw data")
    p <- merge(r$person, t, by = "uid", all.x = TRUE)
    p <- merge(p, m$meta, by = "uid", all.x = TRUE)
    p <- merge(p, m$sum, by = "uid", all.x = TRUE)
    y <- data.table::fifelse(data.table::between(p$ex_year, 1992, 2023), p$ex_year, data.table::fifelse(data.table::between(p$ndi_year, 
        1992, 2023), p$ndi_year, NA_real_))
    mo <- data.table::fifelse(data.table::between(p$ex_month, 1, 12), p$ex_month, data.table::fifelse(data.table::between(p$ndi_month, 
        1, 12), p$ndi_month, 6))
    fallback <- data.table::as.IDate(sprintf("%04d-%02d-15", as.integer(y), as.integer(mo)))
    idx <- which(is.na(p$death) & is.finite(y))
    if (length(idx)) 
        data.table::set(p, i = idx, j = "death", value = fallback[idx])
    ids <- p[order(uid), .(uid)][, `:=`(id, .I)]
    p <- merge(p, ids, by = "uid")
    w <- merge(r$wave, ids, by = "uid")
    w[, `:=`(urban, m$urban[w, on = .(uid, year), urban])]
    mv <- merge(m$long, ids, by = "uid")
    p[, `:=`(c("uid", "ex_month", "ex_year", "ndi_month", "ndi_year"), NULL)]
    w[, `:=`(uid, NULL)]
    mv[, `:=`(uid, NULL)]
    data.table::setcolorder(p, c("id", setdiff(names(p), "id")))
    data.table::setcolorder(w, c("id", setdiff(names(w), "id")))
    data.table::setcolorder(mv, c("id", setdiff(names(mv), "id")))
    data.table::setcolorder(mv, c("id", "start", "end", "km", "moved", "urban_start", "urban_end", "urban_transition"))
    endpoint <- hrs_endpoint_data(p, w, mv)
    info <- data.table::data.table(item = c("RAND respondents", "Mobility respondents", "Tracker respondents", 
        "RAND-Mobility matched", "RAND-Tracker matched", "Early-route overlap"), n = c(nrow(r$person), data.table::uniqueN(m$meta$uid), 
        nrow(t), sum(r$person$uid %in% m$meta$uid), sum(r$person$uid %in% t$uid), sum(m$meta$early_overlap, na.rm = TRUE)))
    definition <- data.table::data.table(variable = c("md_km", "mt", "move_obs", "endpoint md_km/mt"), definition = c("Sum of observed wave-to-wave migration distances, converted from miles to km", 
        "Number of observed intervals with migration distance >0; not an exact count of all moves", "Number of intervals with a valid migration-distance value", 
        "Accumulated after endpoint entry and strictly before an event year, or through censoring"))
    list(person = p, wave = w, move = mv, endpoint = endpoint, info = info, definition = definition)
}

hrs_validate <- function(x) {
    expected <- list(person = c("id", "sex", "race", "hispanic", "education", "degree", "death", "death_year", 
        "first", "last", "study", "first_iw", "birth_region", "school_region", "early_overlap", "move_obs", "mt", 
        "md_km", "first_move", "last_move"), wave = c("id", "wave", "date", "age", "bmi", "smoke_ever", "smoke", 
        "alcohol", "income", "wealth", "physical", "lifting", "stooping", "heart", "stroke", "diabetes", "psychiatric", 
        "dementia", "alzheimer", "cesd", "year", "depression", "urban"), move = c("id", "start", "end", "km", "moved", 
        "urban_start", "urban_end", "urban_transition"), endpoint = c("id", "outcome", "entry", "exit", "event", 
        "prevalent", "visits", "move_obs", "md_km", "mt", "follow_years", "age", "bmi", "smoke", "alcohol", "income", 
        "wealth", "physical", "lifting", "stooping", "sex", "race", "hispanic", "education"), info = c("item", 
        "n"), definition = c("variable", "definition"))
    if (!identical(names(x), names(expected))) 
        stop("Unexpected HRS output tables.", call. = FALSE)
    for (n in names(expected)) {
        if (!data.table::is.data.table(x[[n]]) || !identical(names(x[[n]]), expected[[n]])) 
            stop("Unexpected HRS output schema: ", n, call. = FALSE)
    }
    hrs_unique(x$person, "id", "person table")
    hrs_unique(x$wave, c("id", "wave"), "wave table")
    hrs_unique(x$move, c("id", "start", "end"), "move table")
    hrs_unique(x$endpoint, c("id", "outcome"), "endpoint table")
    if (!identical(x$person$id, seq_len(nrow(x$person)))) 
        stop("Invalid HRS integer IDs.", call. = FALSE)
    for (n in c("wave", "move", "endpoint")) {
        if (anyNA(x[[n]]$id) || !all(x[[n]]$id %in% x$person$id)) 
            stop("Unmatched HRS IDs in ", n, call. = FALSE)
    }
    dates <- list(person = c("death", "first", "last"), wave = "date", endpoint = c("entry", "exit"))
    for (n in names(dates)) for (v in dates[[n]]) {
        if (!inherits(x[[n]][[v]], "IDate")) 
            stop("Invalid HRS date type: ", n, "/", v, call. = FALSE)
    }
    if (any(!is.finite(x$move$km) | x$move$km < 0) || anyNA(x$move$moved) || !all(x$move$moved == as.integer(x$move$km > 
        0))) 
        stop("Invalid HRS migration values.", call. = FALSE)
    invisible(TRUE)
}


# 6. Save the analysis input
make_hrs <- function(save = TRUE, root = Sys.getenv("MOVE_ROOT", getwd()), file = file.path(root, "data/phe/hrs.rds"), 
    overwrite = FALSE, audit_only = FALSE) {
    audit <- hrs_audit(root)
    if (audit_only) 
        return(audit)
    if (any(audit$status != "ready")) {
        bad <- audit$status != "ready"
        stop("HRS raw audit failed: ", paste(paste(audit$source[bad], audit$status[bad]), collapse = "; "), ". Expected ZIPs under ", 
            file.path(root, "data/phe/hrs"), call. = FALSE)
    }
    if (save && file.exists(file) && !overwrite) 
        stop("HRS output already exists; use overwrite = TRUE explicitly or save = FALSE.", call. = FALSE)
    hrs_require()
    tmp <- tempfile("move_hrs_")
    if (!dir.create(tmp, mode = "0700")) 
        stop("Cannot create HRS temporary directory.", call. = FALSE)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    paths <- setNames(file.path(tmp, audit$member), audit$source)
    for (i in seq_len(nrow(audit))) {
        utils::unzip(audit$archive[i], files = audit$member[i], exdir = tmp)
        if (!file.exists(paths[i]) || is.na(file.info(paths[i])$size) || file.info(paths[i])$size != audit$dta_bytes[i]) 
            stop("Incomplete HRS extraction: ", audit$source[i], call. = FALSE)
    }
    m <- hrs_move_data(paths[["move"]])
    r <- hrs_wave_data(paths[["rand"]])
    tracker <- hrs_read(paths[["track"]], c("HHID", "PN", "EXDEATHMO", "EXDEATHYR", "NMONTH", "NYEAR"))
    x <- hrs_assemble(m, r, tracker)
    hrs_validate(x)
    if (save) {
        if (!dir.exists(dirname(file))) 
            stop("HRS output directory does not exist.", call. = FALSE)
        pending <- tempfile(".hrs_", tmpdir = dirname(file), fileext = ".rds")
        on.exit(unlink(pending), add = TRUE)
        saveRDS(x, pending)
        if (file.exists(file) && !overwrite) 
            stop("HRS output appeared during build; not overwritten.", call. = FALSE)
        if (!file.rename(pending, file)) 
            stop("Cannot finalize HRS RDS; previous output was retained.", call. = FALSE)
        message("Saved HRS: ", file)
    }
    message("HRS: ", nrow(x$person), " persons; ", nrow(x$wave), " waves; ", nrow(x$move), " observed mobility intervals; ", 
        nrow(x$endpoint), " endpoints.")
    invisible(x)
}

if (sys.nframe() == 0L) {
    make_hrs(overwrite = TRUE)
}

