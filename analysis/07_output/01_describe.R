# 1. Sample and flow descriptions
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

use(c("assoc", "external", "region"))


# 2. Count selection steps
ed1_data <- function(gwas = NULL) {
    message("ED1: UKB selection counts")
    x <- as.data.table(readRDS(file.path(phe, "all.rds")))
    md0 <- ukb_md(x)
    initial <- 502316L
    white <- sum(x$ethnic.c == "White", na.rm = TRUE)
    geographic <- sum(x$ethnic.c == "White" & is.finite(md0), na.rm = TRUE)
    rm(x, md0)
    invisible(gc())
    stored <- nrow(dat())
    d <- fig2_data()
    v <- c("md_z", "dtdi_z", fig2_cov1)
    core <- d[complete.cases(d[, ..v])]
    ses_n <- sum(is.finite(core$ses_z))
    ncore <- nrow(core)
    f <- readRDS(file.path(out, "Fig2.rds"))
    health <- unique(as.data.table(f$forest$main)[, .(outcome, n, event)])
    stopifnot(all(f$ses$main$n == ses_n), max(health$n) == ncore)
    stages <- function(cohort, label, n) data.table(cohort, stage = seq_along(label), criterion = label, n, excluded = c(NA_integer_, 
        head(n, -1) - tail(n, -1)))
    if (is.null(gwas)) 
        gwas <- rbindlist(lapply(c("migration_distance_pc40", "delta_tdi_ukb_cond"), function(name) {
            file <- file.path(root, "data/gwas/clean", paste0(name, ".gz"))
            cmd <- sprintf(paste0("gzip -cd %s | awk 'BEGIN{FS=\"\\t\"} ", "NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} ", 
                "{n=$h[\"N\"]+0;if(!count||n<lo)lo=n;if(n>hi)hi=n;count++} ", "END{print count,lo,hi}'"), shQuote(file))
            a <- fread(cmd = cmd, header = FALSE)
            stopifnot(nrow(a) == 1L, a$V1 > 0, a$V2 > 0, a$V2 == a$V3)
            data.table(trait = if (name == "migration_distance_pc40") 
                "MD"
            else dtdi_lab, file = basename(file), n = a$V3, snps = a$V1, covariates = paste("Age, age squared, sex, centre, genotyping batch, PC1-PC40", 
                if (name == "delta_tdi_ukb_cond") 
                  "and Birthplace TDI"
                else ""))
        }))
    sources <- list(UKB = list(Stages = stages("UKB", c("Available baseline participants", "White ethnicity", "White with valid geographic coordinates", 
        "Prepared phenotypes with complete core covariates", "TTWA-linked model base", "Complete MD, Delta TDI and model covariates"), 
        c(initial, white, geographic, stored, nrow(d), ncore)), GWAS = gwas, SES = data.table(n = ses_n, excluded = ncore - 
        ses_n), Health = health, Definition = data.table(definition = c("Counts describe the current Fig2 phenotypic analyses. N is the number of people.", 
        "The prepared phenotype cohort requires age, sex, centre, genotyping batch, PC1-PC5, BMI, smoking and alcohol.", 
        "The linked model base additionally requires birthplace TDI and baseline date; primary models require both migration exposures.", 
        "SES additionally requires education years, household income and CAMSIS.", "Health samples exclude prevalent cases separately for each disease and require positive follow-up.", 
        "GWAS and phenotype analyses branch from the geographic sample, not from the lifestyle-complete phenotype sample.", 
        "The current PC1-PC40 GWAS sample sizes are read from every row of the final summary statistics, not the legacy GWAS ID files.", 
        "Trait-specific genotype, phenotype and genetic-covariate availability determine the two GWAS samples."))))
    display <- list(UKB = list(n = c(initial, white, geographic, stored, nrow(d), ncore), gwas = gwas, ses = ses_n, 
        health = range(health$n)))
    rm(core, d)
    invisible(gc())
    for (cohort in c("UKHLS", "HRS")) {
        message("ED1: ", cohort, " selection counts")
        x <- readRDS(file.path(phe, paste0(tolower(cohort), ".rds")))
        initial <- if (cohort == "UKHLS") 
            uniqueN(x$panel$pid)
        else nrow(x$person)
        pan <- ext_panel(cohort, x)
        tv <- if (cohort == "UKHLS") 
            rep_ukhls_tv(x)
        else rep_hrs_tv(x)
        eligible <- uniqueN(pan$id)
        history_ids <- unique(tv[move_obs > 0, id])
        history <- length(history_ids)
        req <- c("ses", "age", "sex", "wave", "date", if (cohort == "UKHLS") "region" else "race")
        ses_complete <- uniqueN(pan[complete.cases(pan[, ..req]), id])
        s <- ext_cross(cohort, pan, hispanic = FALSE)
        ses_n <- nrow(s)
        ses_migrant <- sum(s$mover == 1)
        ses_md <- sum(s$mover == 1 & s$distance > 0)
        ses_non <- sum(s$mover == 0)
        d <- rep_tv_risk(tv, "All-cause mortality")
        valid_followup <- uniqueN(d$id)
        d <- d[move_obs > 0]
        migration_followup <- uniqueN(d$id)
        d[, `:=`(distance = cumulative, mover = ever)]
        req <- unique(c("distance", "count", "mover", "tstart", "tstop", "event", all.vars(as.formula(paste("~", 
            fig4_cov(cohort, "Mortality", "MD"))))))
        d <- d[complete.cases(d[, ..req])]
        mortality <- uniqueN(d$id)
        deaths <- sum(d$event)
        mortality_md <- uniqueN(d[mover == 1 & distance > 0, id])
        deaths_md <- sum(d[mover == 1 & distance > 0, event])
        people <- d[, .(migrant = as.integer(any(mover == 1)), positive = any(distance > 0), events = sum(event)), 
            by = id]
        mortality_migrant <- sum(people$migrant == 1)
        mortality_non <- sum(people$migrant == 0)
        groups <- people[, .(n = .N, events = sum(events)), by = migrant]
        groups[, `:=`(criterion, ifelse(migrant == 1, "Ever migrated during eligible risk intervals", "No observed migration during eligible risk intervals"))]
        cohort_key <- cohort
        saved <- as.data.table(readRDS(file.path(out, "Fig4.rds"))$Main)[cohort == cohort_key]
        stopifnot(all(s$id %in% history_ids), eligible >= history, history >= max(ses_n, mortality), saved[outcome == 
            "SES" & exposure == "Migrant", n] == ses_n, saved[outcome == "SES" & exposure == "MD", n] == ses_md, 
            saved[outcome == "Mortality" & exposure == "Migrant", n] == mortality, saved[outcome == "Mortality" & 
                exposure == "Migrant", events] == deaths, saved[outcome == "Mortality" & exposure == "MD", n] == 
                mortality_md)
        display[[cohort]] <- list(n = c(initial, eligible, history), labels = c(if (cohort == "UKHLS") "Adult respondents\nWaves 1-15" else "RAND respondents\nWaves 1-16", 
            if (cohort == "UKHLS") "White participants\nAt least one visit aged 40-69" else "At least one visit\naged 50-69", 
            "Observed migration history\nat an eligible visit"), ses = ses_n, ses_migrant = ses_migrant, ses_non = ses_non, 
            ses_md = ses_md, mortality = mortality, mortality_migrant = mortality_migrant, mortality_non = mortality_non, 
            mortality_md = mortality_md)
        sources[[cohort]] <- list(Stages = stages(cohort, c("Interview participants", "Age / ethnicity eligible at one or more visits", 
            "Observed migration history by an eligible visit"), c(initial, eligible, history)), SES = data.table(criterion = c("SES and core covariates complete at one or more visits", 
            "Final SES analysis with migration history", "MD analysis in positive-distance migrants", "Non-migrants in the SES sample"), 
            n = c(ses_complete, ses_n, ses_md, ses_non)), Mortality = data.table(criterion = c("Valid survival follow-up", 
            "Migration history before a risk interval", "Complete mortality-model data", "MD analysis in positive-distance migrants"), 
            n = c(valid_followup, migration_followup, mortality, mortality_md), events = c(NA_integer_, NA_integer_, 
                deaths, deaths_md)), Groups = groups[, .(criterion, n, events)], Definition = data.table(definition = c("N is unique people, not person-waves. SES and mortality samples are separate overlapping subsets.", 
            "At least one eligible visit is required; the age window is not a restriction on age at first study recruitment.", 
            "SES: latest eligible complete SES visit with observed migration history and core covariates.", "Mortality: observed migration history at a visit predicts the next valid follow-up interval.", 
            if (cohort == "UKHLS") "White participants aged 40-69; age, sex and region adjusted." else "Participants aged 50-69 of all available races; age, sex and race adjusted.", 
            "SES migrant and non-migrant counts partition the final SES sample.", "Mortality counts partition people into ever-migrants and never-observed migrants over eligible risk intervals. A person can contribute intervals in both states.", 
            "MD analysis includes positive-distance migrants. Migrant-status and MT analyses also include non-migrants; reported moves with zero recorded distance are retained for migrant-status analyses.", 
            "BMI, smoking and alcohol are not required for these external-cohort core models.")))
        rm(x, pan, tv, s, d, people)
        invisible(gc())
    }
    list(display = display, source = sources)
}


# 3. Describe the study samples
supp_base <- function() {
    s <- list()
    message("SES construction")
    x <- dat()
    x[, `:=`(eid, as.character(eid))]
    ph <- ttwa_table()
    a <- merge(x, ph, by = "eid", sort = FALSE)
    v <- c("edu.years", "income", "camsis")
    ok <- complete.cases(a[, ..v])
    p <- prcomp(data.matrix(a[ok, ..v]), center = TRUE, scale. = TRUE)
    if (sum(p$rotation[, 1]) < 0) 
        p$rotation[, 1] <- -p$rotation[, 1]
    s$ses <- list(n = nrow(a), complete = sum(ok), available = rbindlist(lapply(v, function(q) data.table(variable = q, 
        n = sum(is.finite(a[[q]]))))), cor = as.data.table(as.table(cor(data.matrix(a[ok, ..v])))), loading = data.table(variable = v, 
        loading = as.numeric(p$rotation[, 1])), variance = data.table(component = paste0("PC", 1:3), variance = p$sdev^2/sum(p$sdev^2)), 
        patterns = a[, .(education = is.finite(edu.years), income = is.finite(income), occupation = is.finite(camsis))][, 
            .(n = .N), by = .(education, income, occupation)])
    rm(a, p)
    gc()
    message("Geographical measurement")
    g0 <- as.data.table(readRDS(file.path(phe, "all.rds")))
    nn <- grep("nn|dist.*tdi|tdi.*dist", names(g0), value = TRUE)
    if (length(nn)) {
        g0[, `:=`(eid, as.character(eid))]
        g0 <- g0[eid %in% x$eid]
        s$geo <- list(nearest = data.table(km = as.numeric(g0[[nn[1]]])/1000))
        s$geo$nearest <- s$geo$nearest[is.finite(km)]
    }
    else s$geo <- list()
    rm(g0)
    gc()
    dir <- direction_data()
    vv <- c("md_z", "dtdi_z", "density_z", "cross_ttwa")
    s$geo$cor <- as.data.table(as.table(cor(data.matrix(dir[, ..vv]), method = "spearman")))
    s$geo$n <- nrow(dir)
    s$geo$distance <- dir[, .(q = c(0.01, 0.25, 0.5, 0.75, 0.99), km = as.numeric(quantile(md, c(0.01, 0.25, 0.5, 
        0.75, 0.99)))), by = .(cross_ttwa)]
    s$geo$tdi <- dir[, .(btdi, rtdi)]
    s$geo$cut <- attr(x, "tdi_cut")
    s$geo$routes <- dir[, .(n = .N), by = .(tdi4, urban_route)]
    rm(dir, x)
    gc()
    message("Regional precision")
    r <- readRDS(region_path("res", "tests.rds"))
    message("Regional components: ", paste(names(r), collapse = ", "))
    s$regional <- list(area = as.data.table(r$area), tests = as.data.table(r$tests))
    if (!nrow(s$regional$area)) {
        stop("Regional area results are absent: ", paste(names(r), collapse = ", "))
    }
    rm(r)
    s$external <- supp_external()
    s
}


# 4. Describe external coverage
supp_external <- function() {
    message("External history coverage")
    external <- list()
    r4 <- readRDS(file.path(out, "Fig4.rds"))
    for (cohort in c("UKHLS", "HRS")) {
        cohort0 <- cohort
        raw <- readRDS(file.path(phe, paste0(tolower(cohort), ".rds")))
        pan <- ext_panel(cohort, raw)
        ses <- ext_cross(cohort, pan, hispanic = FALSE)
        hist <- merge(pan[, .(id, wave, date, md, moved)], ses[, .(id, end_wave)], by = "id")
        mt <- hist[wave <= end_wave & !is.na(md), .(count = sum(moved == 1, na.rm = TRUE)), by = id]
        ses <- merge(ses, mt, by = "id")
        f <- as.formula(paste("~ final_ses + distance + mover + count +", fig4_cov(cohort, "SES")))
        ses <- ses[complete.cases(ses[, all.vars(f), with = FALSE])]
        stopifnot(nrow(ses) == r4$Description[cohort == cohort0, n])
        hh <- hist[id %in% ses$id & wave <= end_wave]
        coverage <- hh[, .(observations = .N, recorded = sum(!is.na(md)), missing = sum(is.na(md))), by = id]
        coverage[, `:=`(fraction, recorded/observations)]
        tv <- if (cohort == "UKHLS") 
            rep_ukhls_tv(raw)
        else rep_hrs_tv(raw)
        mort <- rep_tv_risk(tv, "All-cause mortality")[move_obs > 0]
        mort[, `:=`(distance = cumulative, mover = ever)]
        f <- as.formula(paste("~ tstart + tstop + event + mover +", fig4_cov(cohort, "Mortality", "Migrant")))
        mort <- mort[complete.cases(mort[, all.vars(f), with = FALSE])]
        stopifnot(uniqueN(mort$id) == r4$Main[cohort == cohort0 & outcome == "Mortality" & exposure == "Migrant", 
            n])
        follow <- mort[, .(years = sum(tstop - tstart)/365.25, intervals = .N, death = sum(event)), by = id]
        external[[cohort]] <- list(n_ses = nrow(ses), n_mort = nrow(follow), events = sum(follow$death), intervals = ses[, 
            .(n = .N), by = .(observed = move_obs)], coverage = coverage[, .(n = .N), by = .(fraction)], follow = follow[, 
            .(years)], observed = ses[, .(years = follow)], interval_length = mort[, .(years = (tstop - tstart)/365.25)])
        rm(raw, pan, ses, hist, mt, hh, coverage, tv, mort, follow)
        gc()
    }
    external
}


# 5. Summarize migration-context flow
supp_flow <- function() {
    x <- direction_data()
    vars <- c("btdi", "rtdi", "md", "ses_z", "density_z", "tdi4", "urban_route", "cross_ttwa")
    x <- x[complete.cases(x[, ..vars])]
    quart <- function(v) as.character(cut(v, quantile(v, seq(0, 1, 0.25)), include.lowest = TRUE, labels = paste0("Q", 
        1:4)))
    stages <- list(MD = quart(x$md), `Birthplace TDI` = quart(x$btdi), `Recruitment TDI` = quart(x$rtdi), `TDI group` = as.character(x$tdi4), 
        `Density change` = quart(x$density_z), `Rural-urban route` = as.character(x$urban_route), `Cross TTWA` = ifelse(x$cross_ttwa == 
            1, "Cross TTWA", "Within TTWA"), SES = quart(x$ses_z))
    levels <- lapply(stages, function(v) if (all(v %in% paste0("Q", 1:4))) 
        paste0("Q", 1:4)
    else if (all(v %in% tdi_lev)) 
        tdi_lev
    else if (any(grepl("Urban|Rural", v))) 
        c("Urban to Urban", "Urban to Rural", "Rural to Urban", "Rural to Rural")
    else c("Within TTWA", "Cross TTWA"))
    nodes <- rbindlist(lapply(seq_along(stages), function(i) {
        a <- data.table(group = stages[[i]])[, .(N = .N), by = group]
        a <- a[match(levels[[i]], group)]
        a[, `:=`(axis = i, stage = names(stages)[i], height = 0.84 * N/nrow(x))]
        gap <- 0.16/(nrow(a) - 1)
        a[, `:=`(top, 1 - cumsum(shift(height, fill = 0)) - (seq_len(.N) - 1) * gap)]
        a[, `:=`(bottom, top - height)]
        a
    }))
    links <- rbindlist(lapply(seq_len(length(stages) - 1), function(i) {
        a <- data.table(quartile = stages[[1]], from = stages[[i]], to = stages[[i + 1]])[, .(N = .N), by = .(quartile, 
            from, to)]
        a[, `:=`(axis = i, height = 0.84 * N/nrow(x))]
        a[, `:=`(fi = match(from, levels[[i]]), ti = match(to, levels[[i + 1]]))]
        setorder(a, fi, quartile, ti)
        a[, `:=`(left, nodes[axis == i, top][fi] - cumsum(height) + height), by = from]
        setorder(a, ti, quartile, fi)
        a[, `:=`(right, nodes[axis == i + 1, top][ti] - cumsum(height) + height), by = to]
        a
    }))
    links[, `:=`(id, .I)]
    numeric <- c(`1` = "md", `2` = "btdi", `3` = "rtdi", `5` = "density_z", `8` = "ses_z")
    list(nodes = nodes, links = links, n = nrow(x), cuts = rbindlist(lapply(c(1, 2, 3, 5, 8), function(i) data.table(Stage = names(stages)[i], 
        Quantile = c("0%", "25%", "50%", "75%", "100%"), Value = as.numeric(quantile(x[[numeric[as.character(i)]]], 
            seq(0, 1, 0.25)))))))
}

if (sys.nframe() == 0L) {
    s <- if (file.exists(file.path(out, "Supp.rds"))) 
        readRDS(file.path(out, "Supp.rds"))
    else list()
    s$base <- supp_base()
    s$flow <- supp_flow()
    saveRDS(s, file.path(out, "Supp.rds"), compress = "gzip")
}

