# 1. Supplementary Tables
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

use(c("assoc", "corr", "interaction", "region", "external", "omics", "describe", "ses", "genetic", "enrich"))

table_labels <- function(d) {
    d <- copy(as.data.table(d))
    drop <- intersect(c("eid", "IID", "FID", "id", "Source", "source_file"), names(d))
    if (length(drop)) 
        d[, `:=`((drop), NULL)]
    for (v in names(d)) if (is.list(d[[v]])) 
        d[, `:=`((v), NULL)]
    map <- c(outcome = "Outcome", exposure = "Exposure", feature = "Feature ID", n = "N", event = "Events", events = "Events", 
        beta = "Beta", se = "SE", HR = "HR", lo = "95% CI lower", hi = "95% CI upper", lower = "95% CI lower", 
        upper = "95% CI upper", FDR = "q", p = "p", q = "q", a = "Path A beta", b = "Path B log HR", proportion = "Mediation proportion", 
        p_bonf = "p Bonferroni", p_indirect = "p indirect", p_indirect_bonf = "p indirect Bonferroni")
    old <- names(d)
    new <- unname(map[old])
    missing <- is.na(new)
    new[missing] <- tools::toTitleCase(gsub("_", " ", old[missing]))
    i <- grepl("^[pq]($|_)", old)
    new[i] <- paste0(substr(old[i], 1, 1), substring(new[i], 2))
    setnames(d, make.unique(new, sep = " "))
    d
}

stack <- function(x) rbindlist(lapply(names(x), function(n) {
    d <- copy(as.data.table(x[[n]]))
    d[, `:=`(Analysis, n)]
    d
}), fill = TRUE)

samples <- function(omics, ses) {
    rows <- list()
    add <- function(cohort, sample, label, n, value, unit = NA_character_, missing = 0L) {
        rows[[length(rows) + 1L]] <<- data.table(Cohort = cohort, Sample = sample,
            Characteristic = label, N = n, Missing = missing, Value = as.character(value), Unit = unit)
    }
    summary <- function(cohort, sample, label, v, unit) {
        n <- length(v)
        v <- v[is.finite(v)]
        if (!length(v)) return(invisible(NULL))
        values <- c(sprintf('%.2f (%.2f)', mean(v), sd(v)),
            sprintf('%.2f (%.2f-%.2f)', median(v), quantile(v, .25), quantile(v, .75)),
            sprintf('%.2f-%.2f', min(v), max(v)))
        for (i in 1:3) add(cohort, sample, paste(label, c('mean (SD)', 'median (IQR)', 'range')[i]),
            n, values[i], unit, n - length(v))
    }
    categories <- function(cohort, sample, label, v) {
        for (level in unique(v[!is.na(v)])) add(cohort, sample, paste0(label, ': ', level),
            length(v), sprintf('%s (%.1f%%)', sum(v == level, na.rm = TRUE),
                100 * mean(v == level, na.rm = TRUE)), 'n (%)', sum(is.na(v)))
    }
    all <- as.data.table(readRDS(file.path(phe, 'all.rds')))
    add('UKB', 'Available baseline participants', 'Participants', 502316L, 502316L)
    nw <- sum(all$ethnic.c == 'White', na.rm = TRUE)
    add('UKB', 'White participants', 'Participants', nw, nw)
    rm(all); invisible(gc())
    x <- fig2_data()
    x <- x[complete.cases(x[, c('md_z', 'dtdi_z', fig2_cov1), with = FALSE])]
    add('UKB', 'Core phenotype sample', 'Participants', nrow(x), nrow(x))
    vars <- c(Age = 'age', MD = 'md', setNames('dtdi', dtdi_lab),
        'Birthplace TDI' = 'btdi', 'Recruitment TDI' = 'rtdi', BMI = 'bmi',
        'Education years' = 'edu.years', 'Income category' = 'income', CAMSIS = 'camsis', SES = 'ses_z')
    units <- c('years', 'km', rep('TDI units', 3), 'kg/m2', 'years', 'category', 'score', 'SD')
    for (i in seq_along(vars)) summary('UKB', 'Core phenotype sample', names(vars)[i], x[[vars[i]]], units[i])
    categorical <- c(Sex = 'sex', 'TDI group' = 'tdi4', 'Migration cross TTWA' = 'cross_ttwa',
        Smoking = 'smoke', 'Alcohol use' = 'alcohol', 'Occupational physical demand' = 'work_heavy')
    for (label in names(categorical)) {
        v <- categorical[[label]]
        a <- as.character(x[[v]])
        if (v == 'sex') a <- ifelse(a %in% c('1', 'Male', 'M'), 'Men', ifelse(is.na(a), NA, 'Women'))
        if (v == 'cross_ttwa') a <- ifelse(a == '1', 'Cross TTWA', ifelse(a == '0', 'Within TTWA', NA))
        if (v == 'work_heavy') a <- c('Never/rarely', 'Sometimes', 'Usually', 'Always')[as.integer(a)]
        categories('UKB', 'Core phenotype sample', label, a)
    }
    od <- defs(make_icd(x$eid))
    for (o in forest_out) {
        d <- ukb_surv(x, od, o)
        add('UKB', ed_out(o), 'Participants', nrow(d), nrow(d))
        add('UKB', ed_out(o), 'Incident events', nrow(d), sum(d$event))
        summary('UKB', ed_out(o), 'Follow-up', d$time, 'years')
    }
    rm(x, d, od); invisible(gc())
    for (cohort in c('UKHLS', 'HRS')) {
        raw <- readRDS(file.path(phe, paste0(tolower(cohort), '.rds')))
        pan <- ext_panel(cohort, raw)
        d <- ext_cross(cohort, pan, hispanic = FALSE)
        hist <- merge(pan[, .(id, wave, md, moved)], d[, .(id, end_wave)], by = 'id')
        mt <- hist[wave <= end_wave & !is.na(md), .(count = sum(moved == 1, na.rm = TRUE)), by = id]
        d <- merge(d, mt, by = 'id')
        req <- all.vars(as.formula(paste('~ final_ses + distance + mover + count +', fig5_cov(cohort, 'SES'))))
        d <- d[complete.cases(d[, ..req])]
        add(cohort, 'SES sample', 'Participants', nrow(d), nrow(d))
        summary(cohort, 'SES sample', 'Age', d$age, 'years')
        summary(cohort, 'SES sample', 'SES', d$final_ses, 'SD')
        summary(cohort, 'Positive-distance migrants', 'MD', d[mover == 1 & distance > 0, distance], 'km')
        summary(cohort, 'SES sample', 'Observed history duration', d$follow, 'years')
        categories(cohort, 'SES sample', 'Migration status', ifelse(d$mover == 1, 'Migrant', 'Non-migrant'))
        categories(cohort, 'SES sample', 'Sex', ifelse(as.character(d$sex) %in% c('1', 'Male', 'M'), 'Men', 'Women'))
        tv <- if (cohort == 'UKHLS') rep_ukhls_tv(raw) else rep_hrs_tv(raw)
        m <- rep_tv_risk(tv, 'All-cause mortality')[move_obs > 0]
        m[, `:=`(distance = cumulative, mover = ever)]
        req <- all.vars(as.formula(paste('~ tstart + tstop + event + mover +', fig5_cov(cohort, 'Mortality', 'Migrant'))))
        m <- m[complete.cases(m[, ..req])]
        follow <- m[, .(years = sum(tstop - tstart)/365.25), by = id]
        add(cohort, 'Mortality sample', 'Participants', nrow(follow), nrow(follow))
        add(cohort, 'Mortality sample', 'Deaths', nrow(follow), sum(m$event))
        summary(cohort, 'Mortality sample', 'Eligible at-risk follow-up', follow$years, 'years')
        rm(raw, pan, d, hist, mt, tv, m, follow); invisible(gc())
    }
    for (layer in c('Protein', 'MRI')) {
        a <- unique(omics[[layer]]$association[, .(exposure, n)])
        for (i in seq_len(nrow(a))) add('UKB', paste(layer, a$exposure[i]),
            'Association participants', a$n[i], a$n[i])
    }
    for (i in seq_len(nrow(ses$loading))) add('UKB', 'SES derivation sample',
        paste(c('Education years', 'Income', 'CAMSIS')[i], 'PC1 loading'),
        ses$complete, ses$loading$loading[i], 'Loading')
    add('UKB', 'SES derivation sample', 'PC1 variance explained', ses$complete,
        100 * ses$variance$variance[1], '%')
    rbindlist(rows)
}

run_tables <- function() {
    get <- function(n) readRDS(file.path(out, paste0(n, ".rds")))
    specs <- fread(file.path(code, "reference/tables.tsv"))
    tables <- setNames(vector("list", 25), paste0("Table S", 1:25))
    for (n in c(1:3, 21)) tables[[n]] <- fread(file.path(code, "reference", paste0("S", n, ".tsv")))
    r <- get("Fig2")
    s <- get("Supp")
    ed7 <- get("ED7")
    ed8 <- get("ED8")
    ed9 <- get("ED9")
    tables[[4]] <- samples(get('Fig6'), s$base$ses)
    g <- if (file.exists(file.path(out, "genetic.rds"))) 
        get("genetic")
    else list(qc = ed3_qc(), comparison = ed3_comparison())
    tables[[5]] <- stack(list(Quality = g$qc$metrics, Comparisons = rbindlist(lapply(g$comparison, `[[`, "metrics"))))
    for (trait in c("migration_distance_pc40", "delta_tdi_ukb_cond")) {
        x <- fread(file.path(root, "data/cojo", trait, paste0(trait, ".jma.cojo")))
        x[, `:=`(Analysis, paste("COJO", trait))]
        tables[[5]] <- rbindlist(list(tables[[5]], x), fill = TRUE)
    }
    c <- fig2a_data()
    tables[[6]] <- stack(list(Correlations = c$d, Heritability = c$h))
    p <- get("prs")
    tables[[7]] <- stack(list(Sources = p$definitions, Groups = p$bars, Associations = p$associations))
    tables[[8]] <- stack(list(SES = r$ses$main, TDI_groups = r$ses$tdi, SES_components = r$ses$comp, SES_density = r$ses$density, 
        SES_routes = r$ses$routes, MD_TTWA = r$forest$main, TDI = r$forest$tdi, Directions = direction_table(r$direction)))
    tables[[9]] <- stack(list(Main_strata = r$strata, Other_strata = ed8$Extra, Sensitivity = ed7$Effects, Nonlinearity = ed7$Tests, 
        TTWA = ed7$TTWA))
    tables[[10]] <- stack(r$interaction)
    z <- region_read()
    tables[[11]] <- stack(list(Variance = z$tests, Areas = z$area))
    e <- get("Fig5")
    tables[[12]] <- stack(c(e[c("Main", "Curves", "Within")], setNames(ed9[c("Effects", "Curves", "Within")], 
        c("Sensitivity", "Sensitivity_curves", "Sensitivity_within"))))
    info <- readRDS(file.path(phe, "omics_info.rds"))
    tables[[13]] <- rbindlist(info[c("Protein", "MRI")], fill = TRUE)
    o <- get("Fig6")
    tables[[14]] <- rbindlist(lapply(o[c("Protein", "MRI")], `[[`, "association"), fill = TRUE)
    tables[[15]] <- rbindlist(lapply(o[c("Protein", "MRI")], `[[`, "mediation"), fill = TRUE)
    enrich <- enrich_data()
    tables[[16]] <- stack(enrich[c("tissue", "cells")])
    tables[[17]] <- stack(list(Contrasts = s$omics$Protein$contrasts, Profiles = s$pattern$profiles))
    tables[[18]] <- readRDS(file.path(root, "analysis/enrich/extra/cells.rds"))
    tables[[19]] <- stack(s$flow[c("nodes", "links", "cuts")])
    tables[[20]] <- rbindlist(lapply(c("cs8", "dlpfc"), function(atlas) rbindlist(lapply(c("MD", "DS"), function(trait) {
        d <- fread(file.path(root, "analysis/enrich/spatial", atlas, "cauchy_combination", paste0(atlas, "_", trait, 
            ".Cauchy.csv.gz")))
        d[, `:=`(Atlas = atlas, Trait = trait)]
        d
    }))), fill = TRUE)
    tables[[22]] <- supp_region_table(readRDS(file.path(out, "supp/S6.rds")))
    tables[[23]] <- s$matching$effects
    tables[[24]] <- s$adjustment
    tables[[25]] <- get("Fig3")$points
    stopifnot(all(lengths(tables) > 0))
    tables <- lapply(tables, table_labels)
    tables <- tables[c(4,1,5,6,7,8,9,23,24,25,10,11,12,14,17,15,16,21,2,3,13,20,19,22,18)]
    names(tables) <- paste0("Table S", seq_along(tables))
    contents <- data.table(Table = names(tables), Title = specs$Title[match(names(tables), specs$Table)], Rows = vapply(tables, 
        nrow, 0L))
    save_xlsx(c(list(Contents = contents), tables), file.path(out, "Supplementary Tables.xlsx"))
}

if (sys.nframe() == 0L) run_tables()

