# 1. External associations and SES trajectories
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

use(c("ukb", "ukhls", "hrs"))

suppressPackageStartupMessages({
    library(survival)
    library(lme4)
})

rep_ukhls_panel <- function(x) {
    d <- copy(as.data.table(x$panel))
    setorder(d, pid, wave)
    d[, `:=`(move_obs = cumsum(!is.na(md)), md_cum = cumsum(fifelse(is.na(md), 0, md)), mt_cum = cumsum(fifelse(!is.na(moved) & 
        moved == 1, 1, 0)), linked3 = wave == shift(wave) + 1 & shift(wave) == shift(wave, 2) + 1, md3 = frollsum(md, 
        3, align = "right")), by = pid]
    d
}

rep_ukhls_tv <- function(x) {
    d <- rep_ukhls_panel(x)
    setnames(d, "pid", "id")
    wd <- d[!is.na(date), .(wave_date = as.IDate(round(median(as.numeric(date))))), by = wave]
    d <- merge(d, wd, by.x = "death_wave", by.y = "wave", all.x = TRUE, sort = FALSE)
    d[, `:=`(death_date, fifelse(death == 1, wave_date, as.IDate(NA)))]
    setorder(d, id, wave)
    nxt <- c("date", paste0("disease_", c("cad_mi", "stroke", "diabetes", "depression", "anxiety", "dementia")))
    d[, `:=`(paste0("next_", nxt), shift(.SD, type = "lead")), by = id, .SDcols = nxt]
    d[, `:=`(recent = md, rolling3 = md3, cumulative = md_cum, count = mt_cum, distance_rate = fifelse(move_obs > 
        0, md_cum/move_obs, NA_real_), move_mean = fifelse(mt_cum > 0, md_cum/mt_cum, NA_real_), move_rate = fifelse(move_obs > 
        0, mt_cum/move_obs, NA_real_), move = as.integer(moved == 1), ever = as.integer(mt_cum > 0), alcohol = fcase(alcohol_ever_locf == 
        1, 1L, alcohol_ever_locf == 2, 0L), work = factor(fcase(job_physical_locf == 4, "None/almost none", job_physical_locf == 
        3, "Some", job_physical_locf == 2, "Most", job_physical_locf == 1, "All/almost all"), c("None/almost none", 
        "Some", "Most", "All/almost all")))]
    d[age %between% c(40, 69) & ethnicity %in% 1:4]
}

rep_hrs_tv <- function(x) {
    d <- copy(as.data.table(x$wave))
    m <- as.data.table(x$move)
    p <- as.data.table(x$person)
    a <- m[, .(recent = sum(km), moved = as.integer(any(moved == 1))), by = .(id, year = end)]
    d <- merge(d, a, by = c("id", "year"), all.x = TRUE, sort = FALSE)
    d <- merge(d, p[, .(id, sex, race, hispanic, death)], by = "id", all.x = TRUE, sort = FALSE)
    setorder(d, id, date)
    d[, `:=`(move_obs = cumsum(!is.na(recent)), cumulative = cumsum(fifelse(is.na(recent), 0, recent)), count = cumsum(fifelse(!is.na(moved) & 
        moved == 1, 1, 0))), by = id]
    d[, `:=`(rolling3, frollsum(recent, 3, align = "right")), by = id]
    nxt <- c("date", "heart", "stroke", "diabetes", "depression", "psychiatric", "dementia")
    d[, `:=`(paste0("next_", nxt), shift(.SD, type = "lead")), by = id, .SDcols = nxt]
    d[, `:=`(distance_rate = fifelse(move_obs > 0, cumulative/move_obs, NA_real_), move_mean = fifelse(count > 
        0, cumulative/count, NA_real_), move_rate = fifelse(move_obs > 0, count/move_obs, NA_real_), move = as.integer(moved == 
        1), ever = as.integer(count > 0), death_date = death, work = factor(fcase(physical == 4, "None/almost none", 
        physical == 3, "Some", physical == 2, "Most", physical == 1, "All/almost all"), c("None/almost none", "Some", 
        "Most", "All/almost all")))]
    d[age %between% c(50, 69)]
}

rep_tv_risk <- function(d, outcome, var = NULL) {
    if (is.null(var)) {
        x <- d[!is.na(date) & (is.na(death_date) | date < death_date)]
        x[, `:=`(stop, fifelse(!is.na(death_date) & (is.na(next_date) | death_date <= next_date), death_date, next_date))]
        x[, `:=`(event, as.integer(!is.na(death_date) & !is.na(stop) & stop == death_date))]
    }
    else {
        first <- d[!is.na(get(var)), .SD[1], by = id][get(var) == 0, id]
        x <- d[id %in% first & get(var) == 0 & !is.na(get(paste0("next_", var)))]
        x[, `:=`(stop = get("next_date"), event = as.integer(get(paste0("next_", var)) == 1))]
        x[!is.na(death_date) & death_date < stop, `:=`(stop = death_date, event = 0L)]
    }
    x <- x[!is.na(stop) & stop > date]
    setorder(x, id, date)
    x <- x[, {
        j <- which(event == 1)[1]
        if (is.na(j)) 
            .SD
        else .SD[seq_len(j)]
    }, by = id]
    x[, `:=`(tstart = as.numeric(date), tstop = as.numeric(stop), outcome = outcome)]
    x
}

ext_rz <- function(x) {
    x <- as.numeric(x)
    i <- is.finite(x)
    if (sum(i) < 3) 
        return(rep(NA_real_, length(x)))
    q <- quantile(x[i], c(0.01, 0.99))
    x[i] <- pmin(pmax(x[i], q[1]), q[2])
    z(x)
}


# 2. Harmonize the cohort panels
ext_panel <- function(cohort, x) {
    if (cohort == "UKHLS") {
        d <- copy(as.data.table(x$panel))[age %between% c(40, 69) & white == 1]
        d[, `:=`(edu, fcase(qualification == 1, 5, qualification == 2, 4, qualification == 3, 3, qualification == 
            4, 2, qualification == 5, 1, qualification == 9, 0, default = NA_real_))]
        d[, `:=`(inc, log1p(pmax(hh_income, 0)))]
        setnames(d, "pid", "id")
        d[, `:=`(race = NA_integer_, hispanic = NA_integer_)]
    }
    else {
        d <- merge(copy(as.data.table(x$wave)), as.data.table(x$person)[, .(id, edu = education, sex0 = sex, race, 
            hispanic)], by = "id", all.x = TRUE)
        a <- as.data.table(x$move)[, .(md = sum(km), moved = as.integer(any(moved == 1))), by = .(id, year = end)]
        d <- merge(d, a, by = c("id", "year"), all.x = TRUE)[age %between% c(50, 69)]
        d[, `:=`(sex = sex0, inc = log1p(pmax(income, 0)))]
    }
    d[, `:=`(edu_z = ext_rz(edu), inc_z = ext_rz(inc)), by = wave]
    d[, `:=`(ses, rowMeans(cbind(edu_z, inc_z), na.rm = FALSE))]
    d[, `:=`(ses, ext_rz(ses)), by = wave]
    setorder(d, id, wave)
    d
}

ext_cross <- function(cohort, d, hispanic = TRUE) {
    v <- c("ses", "age", "sex", "wave", "date", if (cohort == "UKHLS") "region" else c("race", if (hispanic) "hispanic"))
    s <- d[complete.cases(d[, ..v])]
    f <- s[, .SD[.N], by = id][, .(id, end_wave = wave, end_date = date, final_ses = ses, age, sex, region = if (cohort == 
        "UKHLS") 
        region
    else NA_real_, race, hispanic)]
    a <- merge(d[, .(id, wave, date, md, moved)], f[, .(id, end_wave)], by = "id")
    a <- a[wave <= end_wave & !is.na(md), .(distance = sum(md), move_obs = .N, mover = as.integer(any(moved == 
        1, na.rm = TRUE)), first_date = min(date)), by = id]
    x <- merge(f, a, by = "id")[move_obs > 0]
    x[, `:=`(follow, as.numeric(end_date - first_date)/365.25)]
    x[is.finite(follow) & follow >= 0]
}


# 3. Fix cohort-specific covariates
fig5_cov <- function(cohort, outcome, exposure = "MD") {
    geo <- if (cohort == "UKHLS") 
        "factor(region)"
    else "factor(race)"
    paste("age + sex +", geo, "+", if (outcome == "SES") 
        "factor(end_wave) + follow + move_obs"
    else paste("strata(wave) + cluster(id)", if (exposure != "Migrant") 
        "+ move_obs"
    else ""))
}

fig5_model <- function(d, cohort, outcome, terms, exposure = "MD") {
    f <- as.formula(paste(if (outcome == "SES") 
        "final_ses ~"
    else "Surv(tstart,tstop,event) ~", terms, "+", fig5_cov(cohort, outcome, exposure)))
    d <- d[complete.cases(d[, all.vars(f), with = FALSE])]
    w <- character()
    m <- withCallingHandlers(if (outcome == "SES") 
        lm(f, d)
    else coxph(f, d, ties = "efron"), warning = function(e) {
        w <<- c(w, conditionMessage(e))
        invokeRestart("muffleWarning")
    })
    list(b = coef(m), V = if (outcome == "SES") vcovHC(m, type = "HC3") else vcov(m), mean = if (outcome == "SES") colMeans(model.matrix(m)) else NULL, 
        n = uniqueN(d$id), observations = nrow(d), events = if (outcome == "SES") NA_integer_ else sum(d$event), 
        status = if (length(w)) paste(w, collapse = "; ") else "OK", formula = paste(deparse(f), collapse = " "))
}

fig5_basis <- function(v, exposure) {
    if (exposure == "MD") {
        k <- quantile(v, c(0.05, 0.35, 0.65, 0.95))
        if (anyDuplicated(k)) 
            k <- c(min(v), quantile(v[v > min(v)], c(1/3, 2/3)), max(v))
        splines::ns(v, knots = k[2:3], Boundary.knots = k[c(1, 4)])
    }
    else splines::ns(v, knots = c(1, 2), Boundary.knots = c(0, max(3, max(v))))
}

fig5_predict <- function(b, V, L, outcome = "SES") {
    L <- as.matrix(L)
    e <- drop(L %*% b)
    se <- sqrt(pmax(0, rowSums((L %*% V) * L)))
    d <- data.table(estimate = e, lo = e - 1.96 * se, hi = e + 1.96 * se)
    if (outcome == "Mortality") 
        d[, `:=`(names(d), lapply(.SD, exp))]
    d
}

fig5_test <- function(B, x, b, V) {
    linear <- lm.fit(cbind(1, B), x)
    stopifnot(max(abs(linear$residuals)) < 1e-07 * max(1, max(abs(x))))
    C <- t(qr.Q(qr(matrix(linear$coefficients[-1], ncol = 1)), complete = TRUE)[, -1])
    wald <- function(b, V) pchisq(drop(crossprod(b, solve(V, b))), length(b), lower.tail = FALSE)
    list(p_overall = wald(b, V), p_nonlinear = wald(drop(C %*% b), C %*% V %*% t(C)))
}


# 4. Estimate between-person associations
fig5_assoc <- function(d, cohort, outcome, exposure) {
    s <- copy(if (exposure == "MD") 
        d[mover == 1 & distance > 0]
    else d)
    f <- as.formula(paste("~ distance + count + mover +", fig5_cov(cohort, outcome, exposure)))
    v <- unique(c(all.vars(f), if (outcome == "SES") "final_ses" else c("tstart", "tstop", "event")))
    s <- s[complete.cases(s[, ..v])]
    s[, `:=`(xx, switch(exposure, Migrant = as.numeric(mover), MD = z(log1p(distance)), MT = z(count)))]
    m <- fig5_model(s, cohort, outcome, "xx", exposure)
    L <- matrix(as.numeric(names(m$b) == "xx"), nrow = 1)
    a <- fig5_predict(m$b, m$V, L, outcome)
    cbind(data.table(cohort, outcome, exposure, n = m$n, observations = m$observations, events = m$events, p = 2 * 
        pnorm(-abs(m$b["xx"]/sqrt(m$V["xx", "xx"]))), status = m$status, formula = m$formula), a)
}


# 5. Estimate nonlinear curves
fig5_curve <- function(d, cohort, outcome, exposure) {
    f <- as.formula(paste("~ distance + count + mover +", fig5_cov(cohort, outcome)))
    v <- unique(c(all.vars(f), if (outcome == "SES") "final_ses" else c("tstart", "tstop", "event")))
    d <- copy(d[complete.cases(d[, ..v])])
    if (exposure == "MD") 
        d <- d[mover == 1 & distance > 0]
    d[, `:=`(xx, if (exposure == "MD") 
        log1p(distance)
    else count)]
    nm <- paste0("s", 1:3)
    B <- fig5_basis(d$xx, exposure)
    d[, `:=`((nm), as.data.table(B))]
    lim <- quantile(d$xx, c(0.01, 0.99))
    grid <- if (exposure == "MD") 
        seq(lim[1], lim[2], length.out = 160)
    else 0:floor(lim[2])
    ref <- if (exposure == "MD") 
        median(d$xx)
    else 1
    m <- fig5_model(d, cohort, outcome, paste(nm, collapse = " + "))
    stopifnot(all(is.finite(m$b)))
    L0 <- if (outcome == "SES") 
        m$mean
    else setNames(rep(0, length(m$b)), names(m$b))
    L0[nm] <- predict(B, ref)[1, ]
    L <- matrix(L0, nrow = length(grid), ncol = length(m$b), byrow = TRUE, dimnames = list(NULL, names(m$b)))
    L[, nm] <- predict(B, grid)
    A <- sweep(L, 2, L0, "-")
    a <- cbind(data.table(x = if (exposure == "MD") 
        expm1(grid)
    else grid, reference = if (exposure == "MD") 
        expm1(ref)
    else ref), fig5_predict(m$b, m$V, A, outcome))
    test <- fig5_test(B, d$xx, m$b[nm], m$V[nm, nm])
    a[, `:=`(cohort = cohort, outcome = outcome, exposure = exposure, n = m$n, observations = m$observations, events = m$events, 
        p_overall = test$p_overall, p_nonlinear = test$p_nonlinear, status = m$status)]
    a
}


# 6. Estimate within-person SES trajectories
fig5_within <- function(pan, cohort, exposure, adjusted = FALSE) {
    d <- copy(pan)
    setorder(d, id, wave)
    d[, `:=`(obs = cumsum(!is.na(md)), distance = cumsum(fifelse(is.na(md), 0, md)), count = cumsum(fifelse(!is.na(moved) & 
        moved == 1, 1, 0))), by = id]
    d[, `:=`(xx, if (exposure == "MD") 
        distance
    else count)]
    v <- c("ses", "xx", "age", "sex", "wave", if (cohort == "UKHLS") "region" else "race")
    if (adjusted) {
        v <- c(v, "bmi", "smoke")
        d <- d[is.finite(bmi) & bmi >= 10 & bmi <= 80 & !is.na(smoke) & smoke >= 0]
    }
    d <- d[obs > 0 & complete.cases(d[, ..v])]
    ids <- d[, .(n = .N, v = uniqueN(xx), m = max(xx)), by = id][n >= 2 & v >= 2 & m > 0, id]
    d <- d[id %in% ids]
    if (exposure == "MD") 
        d[, `:=`(xx, log1p(xx))]
    B <- fig5_basis(d$xx, exposure)
    bn <- paste0("b", 1:3)
    an <- paste0("a", 1:3)
    wn <- paste0("w", 1:3)
    d[, `:=`((bn), as.data.table(B))]
    for (j in 1:3) {
        bj <- bn[j]
        aj <- an[j]
        wj <- wn[j]
        d[, `:=`((aj), mean(get(bj))), by = id]
        d[, `:=`((wj), get(bj) - get(aj))]
    }
    f <- as.formula(paste("ses ~", paste(c(wn, an), collapse = "+"), "+ age + sex + factor(wave) +", if (cohort == 
        "UKHLS") 
        "factor(region)"
    else "factor(race)", if (adjusted) 
        "+ bmi + factor(smoke)"
    else "", "+ (1|id)"))
    w <- character()
    m <- withCallingHandlers(lmer(f, d, REML = FALSE, control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)), 
        warning = function(e) {
            w <<- c(w, conditionMessage(e))
            invokeRestart("muffleWarning")
        })
    lim <- quantile(d$xx, c(0.01, 0.99))
    grid <- if (exposure == "MD") 
        seq(lim[1], lim[2], length.out = 160)
    else 0:floor(lim[2])
    ref <- if (exposure == "MD") 
        median(d[distance > 0, xx])
    else 1
    L <- sweep(predict(B, grid), 2, predict(B, ref)[1, ], "-")
    b <- fixef(m)[wn]
    V <- as.matrix(vcov(m)[wn, wn])
    test <- fig5_test(B, d$xx, b, V)
    a <- cbind(data.table(cohort, exposure, x = if (exposure == "MD") 
        expm1(grid)
    else grid, reference = if (exposure == "MD") 
        expm1(ref)
    else ref), fig5_predict(b, V, L))
    a[, `:=`(n = uniqueN(d$id), observations = nrow(d), p_overall = test$p_overall, p_nonlinear = test$p_nonlinear, 
        status = if (length(w)) 
            paste(w, collapse = "; ")
        else "OK", singular = isSingular(m), formula = paste(deparse(f), collapse = " "))]
    a
}

fig5_waves <- function(pan, ids, cohort) {
    h <- copy(pan[id %in% ids])
    setorder(h, id, wave)
    h[, `:=`(obs = cumsum(!is.na(md)), cumulative = cumsum(fifelse(is.na(md), 0, md))), by = id]
    w <- h[, {
        v <- cumulative[obs > 0 & cumulative > 0]
        q <- if (length(v)) 
            quantile(v, c(0.25, 0.5, 0.75))
        else rep(NA_real_, 3)
        .(cohort, n = uniqueN(id), available = sum(obs > 0), migrants = length(v), recorded = sum(!is.na(md)), 
            interval_km = sum(md, na.rm = TRUE), mean_md = if (any(obs > 0)) 
                mean(cumulative[obs > 0])
            else NA_real_, md_q1 = q[1], md_median = q[2], md_q3 = q[3])
    }, by = wave]
    setorder(w, wave)
    w[, `:=`(total_km, fifelse(cumsum(recorded) > 0, cumsum(interval_km), NA_real_))]
    w
}


# 7. Save cohort results
run_fig5_fit <- function() {
    setDTthreads(1)
    main <- curves <- within <- status <- counts <- densities <- description <- waves <- list()
    for (cohort in c("UKHLS", "HRS")) {
        message(cohort, ": core models; log MD")
        x <- readRDS(file.path(phe, paste0(tolower(cohort), ".rds")))
        pan <- ext_panel(cohort, x)
        ses <- ext_cross(cohort, pan, hispanic = FALSE)
        hist <- merge(pan[, .(id, wave, md, moved)], ses[, .(id, end_wave)], by = "id")
        mt <- hist[wave <= end_wave & !is.na(md), .(count = sum(moved == 1, na.rm = TRUE)), by = id]
        ses <- merge(ses, mt, by = "id")
        mort <- rep_tv_risk(if (cohort == "UKHLS") 
            rep_ukhls_tv(x)
        else rep_hrs_tv(x), "All-cause mortality")[move_obs > 0]
        mort[, `:=`(distance = cumulative, mover = ever)]
        for (o in c("SES", "Mortality")) {
            d <- if (o == "SES") 
                ses
            else mort
            for (e in c("Migrant", "MD", "MT")) main[[paste(cohort, o, e)]] <- fig5_assoc(d, cohort, o, e)
            for (e in c("MD", "MT")) curves[[paste(cohort, o, e)]] <- fig5_curve(d, cohort, o, e)
        }
        for (e in c("MD", "MT")) {
            message(cohort, ": individual SES change, ", e)
            within[[paste(cohort, e)]] <- fig5_within(pan, cohort, e)
        }
        f <- as.formula(paste("~ final_ses + distance + mover + count +", fig5_cov(cohort, "SES")))
        s <- ses[complete.cases(ses[, all.vars(f), with = FALSE])]
        status[[cohort]] <- s[, .(n = .N), by = mover][, `:=`(cohort, cohort)]
        counts[[cohort]] <- s[, .(n = .N), by = count][, `:=`(cohort = cohort, fraction = n/sum(n))]
        v <- s[mover == 1 & distance > 0, distance]
        den <- density(log1p(v), n = 256, from = 0, to = log1p(max(v)))
        densities[[cohort]] <- data.table(cohort, x = den$x, y = den$y)
        waves[[cohort]] <- fig5_waves(pan, s$id, cohort)
        description[[cohort]] <- data.table(cohort, n = nrow(s), migrants = sum(s$mover == 1), nonmigrants = sum(s$mover == 
            0), md_median = median(v), md_q1 = quantile(v, 0.25), md_q3 = quantile(v, 0.75), md_max = max(v), mt_max = max(s$count))
    }
    r <- list(Main = rbindlist(main), Curves = rbindlist(curves), Within = rbindlist(within), Description = rbindlist(description), 
        Waves = rbindlist(waves), Status = rbindlist(status), Counts = rbindlist(counts), Density = rbindlist(densities))
    r$Main[, `:=`(q, p.adjust(ifelse(status == "OK", p, NA_real_), "BH")), by = cohort]
    r$Methods <- data.table(item = c("Sample", "SES", "MD", "MT", "Migrant", "Core adjustment", "Mortality", "Curves", 
        "Individual SES change", "Missingness", "Multiplicity", "Display", "Waves", "Version"), definition = c("UKHLS: White, age 40-69; HRS: age 50-69. Latest eligible SES wave with complete core covariates. No eligibility requirement for BMI, smoking, alcohol or Hispanic ethnicity.", 
        "Mean of education and log household-income standardized within wave, then SES standardized within wave; relative socioeconomic position.", 
        "Observed cumulative km in positive-distance migrants only. Natural log(1+km), then unit SD for forest coefficients, computed once in each cohort/outcome complete analysis set, not separately by wave. Cox uses interval observations as analysis units. Splines use log(1+km) directly.", 
        "Observed cumulative move count; per SD, including zero. Missing intervals are not counted as confirmed non-moves.", 
        "At least one observed residential move, not necessarily migration across a city or labour-market boundary.", 
        "Age, sex, region (UKHLS) or race (HRS). SES additionally wave, observation duration and observed intervals. Cox uses wave strata, participant-cluster robust SE; MD/MT additionally observed intervals. No BMI, smoking, alcohol or Hispanic ethnicity.", 
        "Time-updated Cox: observed migration history at an interview predicts the next follow-up interval. UKHLS death dates approximate death-wave median interview date. HRS uses recorded death date.", 
        "Natural cubic splines: log(1+MD) among migrants, knots at 5/35/65/95 percentiles. If ties duplicate knots, use full-range boundaries and tertiles above the tied minimum. MT includes zero, internal knots 1/2, boundary 0/max count. SES shows adjusted contrasts, not absolute predicted SES. SES/Cox reference: cohort/outcome-specific median MD among migrants or MT=1; no non-migrant reference point.", 
        "Mixed model with person-centred spline bases, person-mean bases, age/sex/wave/geography and random intercept. >=2 complete waves and changing exposure. MD uses log(1+km), including pre-move zero states of eventual movers. Reference: median positive-distance person-wave MD or MT=1. Relative SES change, not a person's absolute SES or a causal trajectory. Overall/nonlinear tests concern within-person coefficients conditional on person-mean spline bases.", 
        "No new imputation. Cumulative histories sum observed intervals only; wholly unobserved histories excluded. Observed non-moves distinguished from missing records.", 
        "BH-FDR separately within each cohort over 3 exposures x 2 outcomes (6 primary tests). Filled forest points q<0.05; open otherwise. Curve p overall is a 3-df Wald test of spline coefficients; p non-linear is a 2-df Wald test of departure from linearity in log(1+MD) or MT, using the same fitted model and covariance. HC3 for cross-sectional SES, participant-robust for Cox, model-based for mixed models. Curve tests are exploratory and unadjusted.", 
        "All complete observations fitted. MD curves show 1st-99th percentiles, km labels on a log(1+km) axis; MT predictions at integer counts through 99th percentile, joined by lines. Panel a describes core-model SES sample. Intervals/bands are pointwise 95% CI. Descriptive wave/count lines do not receive inferential p values.", 
        fig5_wave_note, "2026-09-21 core adjustment, log MD, reference-centred curves and Wald tests, 3 rows x 4 columns"))
    saveRDS(r, file.path(out, "Fig5.rds"))
    print(r$Main[, .(cohort, outcome, exposure, n, events, estimate, lo, hi, p, q, status)])
    invisible(r)
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

urban_label <- function(from, to) {
    fifelse(from %in% 0:1 & to %in% 0:1, paste(ifelse(from == 1, "Urban", "Rural"), "to", ifelse(to == 1, "urban", 
        "rural")), NA_character_)
}

urban_path_data <- function() {
    paths <- coverage <- list()
    for (cohort in c("UKHLS", "HRS")) {
        x <- readRDS(file.path(phe, paste0(tolower(cohort), ".rds")))
        if (cohort == "UKHLS") {
            d <- as.data.table(x$panel)[, .(id = pid, wave, urban, moved)]
        }
        else {
            d <- as.data.table(x$wave)[, .(id, wave, year, urban)]
            m <- as.data.table(x$move)[, .(moved = as.integer(any(moved == 1))), by = .(id, year = end)]
            d <- merge(d, m, by = c("id", "year"), all.x = TRUE, sort = FALSE)
        }
        setorder(d, id, wave)
        switch <- d[!is.na(urban), .(id, wave, urban, moved)]
        switch[, `:=`(prior, shift(urban)), by = id]
        switch <- switch[!is.na(prior) & urban != prior]
        a <- d[!is.na(urban), .(observed = .N, sequence = paste(rle(urban)$values, collapse = "")), by = id]
        a <- merge(a, d[, .(move_recorded = any(moved == 1, na.rm = TRUE)), by = id], by = "id", sort = FALSE)
        a <- a[observed >= 2]
        a[, `:=`(path, fcase(sequence == "1", "Always urban", sequence == "0", "Always rural", sequence == "01", 
            "Rural to urban", sequence == "10", "Urban to rural", default = "Multiple switches"))]
        paths[[cohort]] <- a[, .(people = .N, with_recorded_move = sum(move_recorded), without_recorded_move = sum(!move_recorded)), 
            by = path][, `:=`(cohort, cohort)]
        coverage[[cohort]] <- data.table(cohort, people = uniqueN(d$id), waves = nrow(d), waves_with_urban = sum(!is.na(d$urban)), 
            urban_coverage = mean(!is.na(d$urban)), people_with_two_urban_waves = nrow(a), people_with_observed_switch = a[!path %in% 
                c("Always urban", "Always rural"), .N], switch_without_recorded_move = a[!path %in% c("Always urban", 
                "Always rural") & !move_recorded, .N], switch_intervals = nrow(switch), switch_with_move = sum(switch$moved == 
                1, na.rm = TRUE), switch_without_move = sum(switch$moved == 0, na.rm = TRUE), switch_move_missing = sum(is.na(switch$moved)))
    }
    list(Paths = rbindlist(paths), Coverage = rbindlist(coverage))
}

urban_sample <- function(x, cohort, outcome) {
    if (outcome == "SES") {
        pan <- urban_history(ext_panel(cohort, x), "id")
        d <- ext_cross(cohort, pan, hispanic = FALSE)
        route <- pan[, .(id, end_wave = wave, urban_route)]
        d <- merge(d, route, by = c("id", "end_wave"), all.x = TRUE, sort = FALSE)
    }
    else {
        pan <- if (cohort == "UKHLS") 
            rep_ukhls_tv(x)
        else rep_hrs_tv(x)
        pan <- urban_history(pan, "id")
        d <- rep_tv_risk(pan, "All-cause mortality")
        d[, `:=`(distance, cumulative)]
    }
    d <- d[urban_route %in% urban_levels & is.finite(distance) & distance > 0]
    d[, `:=`(route = factor(urban_route, urban_levels), md_log = log1p(distance))]
    f <- as.formula(paste(if (outcome == "SES") 
        "final_ses ~"
    else "Surv(tstart,tstop,event) ~", "route + md_log +", fig5_cov(cohort, outcome, "MD")))
    d[complete.cases(d[, all.vars(f), with = FALSE])]
}

urban_fit <- function(d, cohort, outcome) {
    m <- fig5_model(d, cohort, outcome, "route + md_log", "MD")
    pair <- data.table(comparison = c("Urban to rural vs Urban to urban", "Rural to urban vs Urban to urban", "Rural to rural vs Urban to urban", 
        "Rural to urban vs Urban to rural"), a = c("Urban to rural", "Rural to urban", "Rural to rural", "Rural to urban"), 
        b = c("Urban to urban", "Urban to urban", "Urban to urban", "Urban to rural"))
    pair[, `:=`(c("estimate", "lower", "upper", "p"), {
        L <- setNames(rep(0, length(m$b)), names(m$b))
        if (a != urban_levels[1]) 
            L[paste0("route", a)] <- 1
        if (b != urban_levels[1]) 
            L[paste0("route", b)] <- -1
        beta <- sum(L * m$b)
        se <- sqrt(drop(crossprod(L, m$V %*% L)))
        v <- beta + c(0, -1.96, 1.96) * se
        if (outcome == "Mortality") 
            v <- exp(v)
        as.list(c(v, 2 * pnorm(-abs(beta/se))))
    }), by = .(comparison, a, b)]
    pair[, `:=`(cohort = cohort, outcome = outcome, n = m$n, observations = m$observations, events = m$events, 
        status = m$status, formula = m$formula)]
    pair[, `:=`(c("a", "b"), NULL)]
    pair
}

urban_interaction <- function(d, cohort, outcome) {
    m <- fig5_model(d, cohort, outcome, "route * md_log", "MD")
    terms <- grep(":md_log$", names(m$b), value = TRUE)
    v <- m$V[terms, terms, drop = FALSE]
    stat <- drop(crossprod(m$b[terms], solve(v, m$b[terms])))
    p_int <- pchisq(stat, length(terms), lower.tail = FALSE)
    rbindlist(lapply(urban_levels, function(route) {
        L <- setNames(rep(0, length(m$b)), names(m$b))
        L["md_log"] <- log(2)
        if (route != urban_levels[1]) 
            L[paste0("route", route, ":md_log")] <- log(2)
        beta <- sum(L * m$b)
        se <- sqrt(drop(crossprod(L, m$V %*% L)))
        a <- beta + c(0, -1.96, 1.96) * se
        if (outcome == "Mortality") 
            a <- exp(a)
        data.table(cohort, outcome, route, n = m$n, events = m$events, estimate = a[1], lower = a[2], upper = a[3], 
            p = 2 * pnorm(-abs(beta/se)), p_interaction = p_int, status = m$status)
    }))
}

if (sys.nframe() == 0L) {
    run_fig5_fit()
}

