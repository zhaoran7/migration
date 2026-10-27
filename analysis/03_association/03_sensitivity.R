# 1. Sensitivity and extra stratified models
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

use(c("assoc", "interaction", "external"))

ed_coef <- function(m, vars, outcome, group) {
    b <- coef(m)[vars]
    se <- sqrt(diag(vcov(m)))[vars]
    data.table(outcome, exposure = unname(med_exp[vars]), group, n = m$n, event = m$nevent, beta = unname(b), se = unname(se), 
        estimate = exp(unname(b)), lower = exp(b - 1.96 * se), upper = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b/se)))
}

ed7_join <- function(parts) {
    stopifnot(all(lengths(parts) > 0), !any(vapply(parts, inherits, logical(1), "try-error")))
    r <- list(Effects = rbindlist(lapply(parts, `[[`, "Effects")), TTWA = rbindlist(lapply(parts, `[[`, "TTWA")), 
        Curves = rbindlist(lapply(parts, `[[`, "Curves")), Tests = rbindlist(lapply(parts, `[[`, "Tests")), Methods = parts[[1]]$Methods)
    r$Effects[, `:=`(q, p.adjust(p, "BH")), by = group]
    r$TTWA[, `:=`(q, p.adjust(p, "BH"))]
    r$Tests[, `:=`(q_overall = p.adjust(p_overall, "BH"), q_nonlinear = p.adjust(p_nonlinear, "BH"))]
    r
}


# 2. Run phenotype sensitivity models
ed7_fit <- function(outcomes = forest_out) {
    if (.Platform$OS.type != "windows" && length(outcomes) > 1 && core > 1) 
        return(ed7_join(mclapply(outcomes, ed7_fit, mc.cores = min(core, 4))))
    x <- fig2_data()
    od <- defs(make_icd(x$eid))
    setDTthreads(1)
    ranks <- dat()[, .(eid = as.character(eid), value = z(rint(md)))]
    x[, `:=`(md_rint_z, ranks$value[match(eid, ranks$eid)])]
    effects <- curves <- tests <- cross <- list()
    for (o in outcomes) {
        message("ED7: ", o)
        original <- ukb_surv(x, od, o)
        base <- c("md_z", "dtdi_z", fig2_cov1)
        required <- c(base, "tdi4")
        d <- original[complete.cases(original[, ..required])]
        main <- coxph(reformulate(base, "Surv(time,event)"), d, ties = "efron")
        main_tdi <- main
        effects[[paste(o, 0, "main")]] <- rbind(ed_coef(main, "md_z", o, "Baseline"), ed_coef(main_tdi, "dtdi_z", 
            o, "Baseline"))
        dr <- copy(d)
        dr[, `:=`(md_z, md_rint_z)]
        mr <- coxph(reformulate(base, "Surv(time,event)"), dr, ties = "efron")
        mtr <- mr
        effects[[paste(o, "RINT")]] <- rbind(ed_coef(mr, "md_z", o, "RINT"), ed_coef(mtr, "dtdi_z", o, "RINT"))
        rm(dr, mr, mtr)
        for (id in c("md", "dtdi", "ttwa")) {
            rhs <- switch(id, md = base, dtdi = c("md_z", "dtdi_z", tdi_cov1), ttwa = c("cross_ttwa", "dtdi_z", 
                cov1))
            s <- if (id == "ttwa") 
                original[complete.cases(original[, ..rhs])]
            else d
            vars <- switch(id, md = "md_z", dtdi = "dtdi_z", ttwa = "cross_ttwa")
            m <- switch(id, md = main, dtdi = main_tdi, ttwa = coxph(reformulate(rhs, "Surv(time,event)"), s, ties = "efron"))
            if (id == "ttwa") 
                effects[[paste(o, 0, id)]] <- ed_coef(m, vars, o, "Baseline")
            if (id == "ttwa") {
                adjusted <- coxph(reformulate(c(rhs, "md_z"), "Surv(time,event)"), s[is.finite(md_z)], ties = "efron")
                cross[[o]] <- ed_coef(adjusted, "cross_ttwa", o, "MD adjusted")
            }
            for (lag in c(2, 5)) {
                u <- copy(s[time > lag])
                u[, `:=`(time, time - lag)]
                mm <- coxph(reformulate(rhs, "Surv(time,event)"), u, ties = "efron")
                effects[[paste(o, lag, id)]] <- ed_coef(mm, vars, o, paste0(lag, " years"))
                rm(mm, u)
                invisible(gc())
            }
            if (id == "ttwa") 
                rm(m, s)
        }
        for (v in c("md_z", "dtdi_z")) {
            k <- as.numeric(quantile(d[[v]], c(0.05, 0.35, 0.65, 0.95)))
            B <- splines::ns(d[[v]], knots = k[2:3], Boundary.knots = k[c(1, 4)])
            bn <- paste0("spline", 1:3)
            d[, `:=`((bn), as.data.table(B))]
            rhs <- c(bn, setdiff(if (v == "md_z") base else c("md_z", "dtdi_z", tdi_cov1), v))
            m <- coxph(reformulate(rhs, "Surv(time,event)"), d, ties = "efron")
            non <- anova(if (v == "md_z") 
                main
            else main_tdi, m, test = "LRT")
            sm <- coef(m)[bn]
            V <- vcov(m)[bn, bn, drop = FALSE]
            p0 <- pchisq(drop(crossprod(sm, solve(V, sm))), 3, lower.tail = FALSE)
            pc <- grep("Pr|P\\(", names(non), value = TRUE)[1]
            t <- data.table(outcome = o, exposure = med_exp[v], p_overall = p0, p_nonlinear = as.numeric(non[[pc]][2]), 
                n = m$n, event = m$nevent, knot1 = k[1], knot2 = k[2], knot3 = k[3], knot4 = k[4])
            tests[[paste(o, v)]] <- t
            g <- seq(quantile(d[[v]], 0.01), quantile(d[[v]], 0.99), length.out = 101)
            L <- predict(B, g) - matrix(predict(B, 0), length(g), 3, byrow = TRUE)
            b <- drop(L %*% sm)
            se <- sqrt(rowSums((L %*% V) * L))
            curves[[paste(o, v)]] <- data.table(outcome = o, exposure = med_exp[v], value = g, estimate = exp(b), 
                lower = exp(b - 1.96 * se), upper = exp(b + 1.96 * se))
            rm(m, B)
            invisible(gc())
        }
        rm(main, main_tdi, d, original)
        invisible(gc())
    }
    e <- rbindlist(effects)
    e[, `:=`(q, p.adjust(p, "BH")), by = group]
    t <- rbindlist(tests)
    t[, `:=`(q_overall = p.adjust(p_overall, "BH"), q_nonlinear = p.adjust(p_nonlinear, "BH"))]
    cc <- rbindlist(cross)
    cc[, `:=`(q, p.adjust(p, "BH"))]
    list(Effects = e, TTWA = cc, Curves = rbindlist(curves), Tests = t, Methods = data.table(item = c("Early events", 
        "Scaling", "Adjustment", "Multiplicity"), definition = c("Condition on remaining event-free and observed at 2 or 5 years; restart follow-up at that time.", 
        "Primary Fig2 uses raw-distance Z scores, not its separately stored RINT variable. RINT sensitivity ranks MD across the full cleaned UKB cohort; exposure SDs are frozen before outcome exclusions.", 
        "MD and Delta TDI are mutually adjusted with fixed disease covariates; TTWA: Delta TDI, Birthplace TDI and fixed disease covariates; additional linear MD in the TTWA sensitivity model.", 
        "BH across 30 exposure-outcome estimates per exclusion, 20 coefficients for MD-RINT sensitivity, and 20 spline tests per test type.")))
}


# 3. Extend strata and interaction results
ed8_fit <- function() {
    r <- readRDS(file.path(out, "Fig2.rds"))
    x <- fig2_data()
    od <- defs(make_icd(x$eid))
    x[, `:=`(sex2, factor(fifelse(as.character(sex) %in% c("1", "Male", "M"), "Male", "Female"), levels = c("Female", 
        "Male")))]
    x[, `:=`(age_group, cut(age, c(40, 50, 60, 70, Inf), right = FALSE, labels = c("40-49", "50-59", "60-69", "70+")))]
    bounds <- rbindlist(lapply(c("md_z", "dtdi_z"), function(v) data.table(exposure = v, lower = quantile(x[[v]], 
        0.01, na.rm = TRUE), upper = quantile(x[[v]], 0.99, na.rm = TRUE))))
    fit <- function(o) {
        message("ED8 strata: ", o)
        d <- ukb_surv(x, od, o)
        rbindlist(lapply(c("sex", "age"), function(id) {
            levels <- if (id == "sex") 
                c("Female", "Male")
            else c("40-49", "50-59", "60-69", "70+")
            v <- if (id == "sex") 
                "sex2"
            else "age_group"
            cv <- if (id == "sex") 
                setdiff(fig2_cov1, "sex")
            else fig2_cov1
            rbindlist(lapply(levels, function(g) {
                a <- copy(d[as.character(get(v)) == g])
                rhs <- c("md_z", "dtdi_z", cv)
                a <- a[complete.cases(a[, ..rhs])]
                f <- intersect(rhs, names(a)[vapply(a, is.factor, logical(1))])
                a[, `:=`((f), lapply(.SD, droplevels)), .SDcols = f]
                rhs <- setdiff(rhs, f[vapply(a[, ..f], nlevels, integer(1)) < 2])
                warn <- character()
                m <- try(withCallingHandlers(coxph(reformulate(rhs, "Surv(time,event)"), a, ties = "efron"), warning = function(w) {
                  warn <<- c(warn, conditionMessage(w))
                  invokeRestart("muffleWarning")
                }), silent = TRUE)
                ok <- !inherits(m, "try-error") && is.finite(coef(m)["md_z"])
                b <- if (ok) 
                  unname(coef(m)["md_z"])
                else NA_real_
                se <- if (ok) 
                  sqrt(vcov(m)["md_z", "md_z"])
                else NA_real_
                data.table(strat = id, outcome = o, group = g, n = nrow(a), event = sum(a$event), HR = exp(b), 
                  lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b/se)), status = if (!ok) 
                    "Fit failed"
                  else if (length(warn)) 
                    paste(unique(warn), collapse = "; ")
                  else "OK", formula = paste(deparse(reformulate(rhs, "Surv(time,event)")), collapse = " "))
            }))
        }))
    }
    if (.Platform$OS.type == "windows" && core > 1) {
        cl <- parallel::makePSOCKcluster(min(core, 2), outfile = "")
        on.exit(parallel::stopCluster(cl), add = TRUE)
        parallel::clusterCall(cl, function(path) {
            Sys.setlocale("LC_CTYPE", "English_United States.utf8")
            Sys.setenv(MOVE_CODE=path)
            source(file.path(path, "00_common/01_base.R"), encoding = "UTF-8")
            use("sensitivity")
        }, code)
        parallel::clusterExport(cl, c("x", "od"), envir = environment())
        environment(fit) <- .GlobalEnv
        extra <- parallel::parLapply(cl, forest_out, fit)
    }
    else extra <- parallel::mclapply(forest_out, fit, mc.cores = max(1, min(core, 2)))
    stopifnot(all(lengths(extra) > 0), !any(vapply(extra, inherits, logical(1), "try-error")))
    a <- rbindlist(extra)
    a[, `:=`(FDR, p.adjust(p, "BH")), by = strat]
    c(r$interaction, list(Strata = r$strata, Extra = a, Bounds = bounds, Methods = data.table(item = c("Stratification", 
        "Age groups", "Multiplicity", "Curves", "Main figure overlap"), definition = c("Separate Cox models within sex/age strata, retaining the Fig2 mutual MD/Delta TDI adjustment and other covariates; sex omitted when constant within a sex stratum. MD SD fixed in the full cohort.", 
        "40-49, 50-59, 60-69 and 70+ years; age retained as a continuous covariate within each stratum.", "BH across each complete stratification family; all 93 omnibus interaction tests corrected together and all 155 individual interaction coefficients corrected together.", 
        "Conditional HR curves derived from the existing linear-exposure interaction models, not newly fitted spline interactions; shown within the 1st-99th percentile exposure range.", 
        "Only FDR-significant interaction/outcome combinations absent from the eight Fig2d curves are displayed. Source data contain plotted panels only."))))
}


# 4. Prepare external sensitivity samples
ed9_base <- function() {
    effects <- weights <- list()
    for (cohort in c("UKHLS", "HRS")) {
        message("ED9: ", cohort)
        raw <- readRDS(file.path(phe, paste0(tolower(cohort), ".rds")))
        pan <- ext_panel(cohort, raw)
        if (cohort == "UKHLS") 
            pan[, `:=`(bmi = bmi_locf, smoke = smoke_locf)]
        ses <- ext_cross(cohort, pan, hispanic = FALSE)
        hist <- merge(pan[, .(id, wave, md, moved)], ses[, .(id, end_wave)], by = "id")
        mt <- hist[wave <= end_wave & !is.na(md), .(count = sum(moved == 1, na.rm = TRUE)), by = id]
        ses <- merge(ses, mt, by = "id")
        lc <- pan[, .(id, end_wave = wave, bmi, smoke)]
        ses <- merge(ses, lc, by = c("id", "end_wave"), all.x = TRUE, sort = FALSE)
        mort <- rep_tv_risk(if (cohort == "UKHLS") 
            rep_ukhls_tv(raw)
        else rep_hrs_tv(raw), "All-cause mortality")[move_obs > 0]
        if (cohort == "UKHLS") 
            mort[, `:=`(bmi = bmi_locf, smoke = smoke_locf)]
        mort[, `:=`(distance = cumulative, mover = ever)]
        weights[[cohort]] <- data.table(cohort, field = if (cohort == "UKHLS") 
            "weight_xw / weight_lw"
        else "Not in cleaned data", available = cohort == "UKHLS", used = FALSE, reason = if (cohort == "UKHLS") 
            "Fallback weights mix samples and wave-specific eligibility; longitudinal design strata/PSU and a consistent target population must be verified before weighted inference."
        else "Cleaned HRS wave data contain no respondent survey weights; no ad hoc weights imputed.")
        for (o in c("SES", "Mortality")) for (e in c("Migrant", "MD", "MT")) {
            d <- copy(if (o == "SES") 
                ses
            else mort)
            if (e == "MD") 
                d <- d[mover == 1 & distance > 0]
            cv <- fig5_cov(cohort, o, e)
            f <- as.formula(paste("~ distance + count + mover +", cv))
            v <- unique(c(all.vars(f), if (o == "SES") "final_ses" else c("tstart", "tstop", "event")))
            d <- d[complete.cases(d[, ..v])]
            d[, `:=`(xx, switch(e, Migrant = as.numeric(mover), MD = z(log1p(distance)), MT = z(count)))]
            matched <- d[is.finite(bmi) & bmi >= 10 & bmi <= 80 & !is.na(smoke) & smoke >= 0]
            for (g in c("Core", "Matched", "Adjusted")) {
                s <- if (g == "Core") 
                  d
                else matched
                lhs <- if (o == "SES") 
                  "final_ses"
                else "Surv(tstart,tstop,event)"
                formula <- as.formula(paste(lhs, "~ xx +", cv, if (g == "Adjusted") 
                  "+ bmi + factor(smoke)"
                else ""))
                warn <- character()
                m <- withCallingHandlers(if (o == "SES") 
                  lm(formula, s)
                else coxph(formula, s, ties = "efron"), warning = function(w) {
                  warn <<- c(warn, conditionMessage(w))
                  invokeRestart("muffleWarning")
                })
                V <- if (o == "SES") 
                  vcovHC(m, type = "HC3")
                else vcov(m)
                b <- unname(coef(m)["xx"])
                se <- sqrt(V["xx", "xx"])
                effects[[paste(cohort, o, e, g)]] <- data.table(cohort, outcome = o, exposure = e, group = g, n = uniqueN(s$id), 
                  observations = nrow(s), event = if (o == "SES") 
                    NA_integer_
                  else sum(s$event), beta = b, se, p = 2 * pnorm(-abs(b/se)), estimate = if (o == "SES") 
                    b
                  else exp(b), lower = if (o == "SES") 
                    b - 1.96 * se
                  else exp(b - 1.96 * se), upper = if (o == "SES") 
                    b + 1.96 * se
                  else exp(b + 1.96 * se), status = if (length(warn)) 
                    paste(unique(warn), collapse = "; ")
                  else "OK", formula = paste(deparse(formula), collapse = " "), scale_mean = if (e == "MD") 
                    mean(log1p(d$distance))
                  else if (e == "MT") 
                    mean(d$count)
                  else NA_real_, scale_sd = if (e == "MD") 
                    sd(log1p(d$distance))
                  else if (e == "MT") 
                    sd(d$count)
                  else NA_real_)
                rm(m)
                invisible(gc())
            }
        }
        rm(raw, pan, ses, mort, hist)
        invisible(gc())
    }
    a <- rbindlist(effects)
    a[, `:=`(q, p.adjust(ifelse(status == "OK", p, NA_real_), "BH")), by = .(cohort, group)]
    ref <- readRDS(file.path(out, "Fig5.rds"))$Main
    check <- merge(a[group == "Core"], ref, by = c("cohort", "outcome", "exposure"), suffixes = c("", "_main"))
    check[, `:=`(n_difference = n - n_main, beta_difference = estimate - estimate_main)]
    if (any(check$n_difference != 0) || any(abs(check$beta_difference) > 1e-06)) 
        stop("ED9 core estimates do not reproduce Fig5.")
    list(Effects = a, Weights = rbindlist(weights), Check = check, Methods = data.table(item = c("Core", "Matched", 
        "Adjusted", "Scaling", "Intervals", "Correction"), definition = c("Exact Fig5 core models; HRS race, not Hispanic origin.", 
        "Core model restricted to BMI and smoking complete cases; BMI 10-80 kg/m2. UKHLS uses latest available prior values because BMI is only observed at wave 1 and smoking begins at wave 6; HRS uses current wave values.", 
        "Same matched observations as Matched; add BMI and categorical smoking, not alcohol.", "Freeze full core-sample means and SDs before restricting observations.", 
        "SES HC3 standard errors; time-updated Cox with participant-cluster robust standard errors.", "BH across six exposure-outcome tests within each cohort/model; models with convergence warnings excluded from significance classification.")))
}

ed9_curve <- function(d, cohort, outcome, exposure) {
    cv <- fig5_cov(cohort, outcome, exposure)
    v <- unique(c(all.vars(as.formula(paste("~ distance + count + mover +", cv))), if (outcome == "SES") "final_ses" else c("tstart", 
        "tstop", "event")))
    d <- copy(d[complete.cases(d[, ..v])])
    if (exposure == "MD") 
        d <- d[mover == 1 & distance > 0]
    d[, `:=`(xx, if (exposure == "MD") 
        log1p(distance)
    else count)]
    B <- fig5_basis(d$xx, exposure)
    nm <- paste0("s", 1:3)
    d[, `:=`((nm), as.data.table(B))]
    ref <- if (exposure == "MD") 
        median(d$xx)
    else 1
    d <- d[is.finite(bmi) & bmi >= 10 & bmi <= 80 & !is.na(smoke) & smoke >= 0]
    lim <- quantile(d$xx, c(0.01, 0.99))
    grid <- if (exposure == "MD") 
        seq(lim[1], lim[2], length.out = 160)
    else 0:floor(lim[2])
    m <- fig5_model(d, cohort, outcome, paste(c(nm, "bmi", "factor(smoke)"), collapse = " + "), exposure)
    stopifnot(all(is.finite(m$b[nm])), all(is.finite(m$V[nm, nm])))
    L <- predict(B, grid) - matrix(predict(B, ref), length(grid), 3, byrow = TRUE)
    test <- fig5_test(as.matrix(d[, ..nm]), d$xx, m$b[nm], m$V[nm, nm])
    a <- cbind(data.table(cohort, outcome, exposure, group = "Adjusted", x = if (exposure == "MD") 
        expm1(grid)
    else grid, reference = if (exposure == "MD") 
        expm1(ref)
    else ref), fig5_predict(m$b[nm], m$V[nm, nm], L, outcome))
    a[, `:=`(n = m$n, observations = m$observations, events = m$events, p_overall = test$p_overall, p_nonlinear = test$p_nonlinear, 
        status = m$status, formula = m$formula)]
    a
}


# 5. Fit external sensitivity models
ed9_fit <- function() {
    file <- file.path(out, "ED9.rds")
    r <- if (file.exists(file) && !stale(file, file.path(phe, c("ukhls.rds", "hrs.rds")))) 
        readRDS(file)
    else ed9_base()
    curves <- within <- list()
    for (cohort in c("UKHLS", "HRS")) {
        message("ED9 curves: ", cohort)
        raw <- readRDS(file.path(phe, paste0(tolower(cohort), ".rds")))
        pan <- ext_panel(cohort, raw)
        if (cohort == "UKHLS") 
            pan[, `:=`(bmi = bmi_locf, smoke = smoke_locf)]
        ses <- ext_cross(cohort, pan, hispanic = FALSE)
        hist <- merge(pan[, .(id, wave, md, moved)], ses[, .(id, end_wave)], by = "id")
        mt <- hist[wave <= end_wave & !is.na(md), .(count = sum(moved == 1, na.rm = TRUE)), by = id]
        ses <- merge(ses, mt, by = "id")
        ses <- merge(ses, pan[, .(id, end_wave = wave, bmi, smoke)], by = c("id", "end_wave"), all.x = TRUE, sort = FALSE)
        mort <- rep_tv_risk(if (cohort == "UKHLS") 
            rep_ukhls_tv(raw)
        else rep_hrs_tv(raw), "All-cause mortality")[move_obs > 0]
        if (cohort == "UKHLS") 
            mort[, `:=`(bmi = bmi_locf, smoke = smoke_locf)]
        mort[, `:=`(distance = cumulative, mover = ever)]
        for (o in c("SES", "Mortality")) for (e in c("MD", "MT")) {
            message(cohort, " ", o, " ", e)
            curves[[paste(cohort, o, e)]] <- ed9_curve(if (o == "SES") 
                ses
            else mort, cohort, o, e)
        }
        for (e in c("MD", "MT")) {
            message(cohort, " individual SES ", e)
            within[[paste(cohort, e)]] <- fig5_within(pan, cohort, e, adjusted = TRUE)
        }
        rm(raw, pan, ses, mort, hist)
        invisible(gc())
    }
    r$Curves <- rbindlist(curves)
    r$Within <- rbindlist(within)
    c0 <- unique(r$Curves[, .(cohort, outcome, exposure, n, observations, events)])
    c0 <- merge(c0, r$Effects[group == "Adjusted", .(cohort, outcome, exposure, n_main = n, observations_main = observations, 
        events_main = event)], by = c("cohort", "outcome", "exposure"))
    stopifnot(nrow(c0) == 8, all(c0$n == c0$n_main), all(c0$observations == c0$observations_main))
    r$CurveCheck <- c0
    r$Methods <- rbind(r$Methods, data.table(item = c("Displayed models", "Curve adjustment", "Within-person sensitivity", 
        "Point markers"), definition = c("Panels b-e show the BMI/smoking-adjusted model. Core and Matched effect estimates remain in source data; panel a compares Core versus Adjusted sample sizes.", 
        "Same Fig5 spline specification, pointwise intervals, core-sample knots and reference; fitted on BMI/smoking-complete observations with BMI and smoking added. No alcohol or Hispanic-origin adjustment.", 
        "Same Fig5 person-centred spline mixed model with BMI and categorical smoking added; require at least two complete visits and changing exposure. Within-person bases/reference are computed in this eligible sample; SES scaling stays frozen from the full wave data.", 
        "Filled forest circles use nominal p<0.05, open otherwise. BH q values are retained in source data and are not the point-shape criterion.")))
    r
}


# 6. Repeat SES and health models
birth_fit <- function() {
    x <- fig2_data()
    od <- defs(make_icd(x$eid))
    parts <- mclapply(c("SES", forest_out), function(o) {
        message("Birthplace TDI sensitivity: ", o)
        cv <- if (o == "SES") 
            fig2_cov0
        else fig2_cov1
        d <- if (o == "SES") 
            x[is.finite(ses_z)]
        else ukb_surv(x, od, o)
        req <- c("md_z", "dtdi_z", cv, "btdi_z")
        d <- d[complete.cases(d[, ..req])]
        rbindlist(lapply(c("Primary", "Birthplace TDI adjusted"), function(label0) {
            rhs <- c("md_z", "dtdi_z", cv, if (label0 != "Primary") "btdi_z")
            m <- if (o == "SES") 
                lm(reformulate(rhs, "ses_z"), d)
            else coxph(reformulate(rhs, "Surv(time,event)"), d, ties = "efron")
            a <- ukb_coef(m, o, if (o == "SES") 
                "Linear"
            else "Cox", c("md", "dtdi"))
            a[, .(outcome, exposure = label, adjustment = label0, parameter = if (o == "SES") 
                "Beta"
            else "HR", n, event, beta, se, estimate, lower, upper, p)]
        }))
    }, mc.cores = min(core, 2))
    stopifnot(length(parts) == 11L, all(lengths(parts) > 0), !any(vapply(parts, inherits, logical(1), "try-error")))
    a <- rbindlist(parts)
    a[, `:=`(q, p.adjust(p, "BH")), by = .(adjustment, family = ifelse(outcome == "SES", "SES", "Health"))]
    a
}


# 7. Repeat models in the matched subset
tdi_match <- function(limit = 1) {
    x <- fig2_data()
    g <- as.data.table(readRDS(ukb_input("move")))
    stopifnot("tdi.birth.nn" %in% names(g), !anyDuplicated(g$eid))
    x[, `:=`(match_km, as.numeric(g$tdi.birth.nn[match(eid, as.character(g$eid))])/1000)]
    x[, `:=`(tdi4f, factor(as.character(tdi4), levels = c("High to low", "Low to low", "High to high", "Low to high")))]
    pair <- data.table(a = c("Low to high", "High to high", "High to low", "Low to high"), b = c("High to low", 
        "Low to low", "High to high", "Low to low"))
    pair[, `:=`(term, paste(a, "vs", b))]
    od <- defs(make_icd(x$eid))
    extract <- function(m, variables, o, sample) {
        iscox <- inherits(m, "coxph")
        b <- coef(m)[variables]
        se <- sqrt(diag(vcov(m)))[variables]
        data.table(outcome = o, sample, term = unname(med_exp[variables]), n = if (iscox) 
            m$n
        else nobs(m), event = if (iscox) 
            m$nevent
        else NA_integer_, beta = unname(b), se = unname(se), estimate = if (iscox) 
            exp(b)
        else b, lower = if (iscox) 
            exp(b - 1.96 * se)
        else b - 1.96 * se, upper = if (iscox) 
            exp(b + 1.96 * se)
        else b + 1.96 * se, p = 2 * pnorm(-abs(b/se)))
    }
    fit <- function(o) {
        message("TDI matching: ", o)
        iscox <- o != "SES"
        d <- if (iscox) 
            ukb_surv(x, od, o)
        else x[is.finite(ses_z)]
        cv <- if (iscox) 
            tdi_cov1
        else tdi_cov0
        required <- c("md_z", "dtdi_z", "tdi4f", cv)
        d <- d[complete.cases(d[, ..required])]
        rbindlist(lapply(c("All participants", paste0("Match <= ", limit, " km")), function(label) {
            a <- if (label == "All participants") 
                d
            else d[is.finite(match_km) & match_km <= limit]
            lhs <- if (iscox) 
                "Surv(time,event)"
            else "ses_z"
            model <- function(rhs) if (iscox) 
                coxph(reformulate(rhs, lhs), a, ties = "efron")
            else lm(reformulate(rhs, lhs), a)
            m <- model(c("md_z", "dtdi_z", cv))
            cont <- extract(m, "dtdi_z", o, label)
            m <- model(c("tdi4f", "md_z", cv))
            b <- coef(m)
            V <- vcov(m)
            groups <- rbindlist(lapply(seq_len(nrow(pair)), function(i) {
                w <- setNames(rep(0, length(b)), names(b))
                for (v in c("a", "b")) if (pair[[v]][i] != levels(a$tdi4f)[1]) 
                  w[paste0("tdi4f", pair[[v]][i])] <- if (v == "a") 
                    1
                  else -1
                est <- sum(w * b)
                se <- sqrt(drop(w %*% V %*% w))
                data.table(outcome = o, sample = label, term = pair$term[i], n = nrow(a), event = if (iscox) 
                  sum(a$event)
                else NA_integer_, beta = est, se, estimate = if (iscox) 
                  exp(est)
                else est, lower = if (iscox) 
                  exp(est - 1.96 * se)
                else est - 1.96 * se, upper = if (iscox) 
                  exp(est + 1.96 * se)
                else est + 1.96 * se, p = 2 * pnorm(-abs(est/se)))
            }))
            rbind(cont, groups)
        }))
    }
    a <- rbindlist(mclapply(c("SES", forest_out), fit, mc.cores = min(core, 3)))
    a[, `:=`(family, ifelse(outcome == "SES", "SES", "Health"))]
    a[, `:=`(q, p.adjust(p, "BH")), by = .(sample, family)]
    required <- c("md_z", "dtdi_z", "tdi4f", tdi_cov1)
    list(effects = a, coverage = x[complete.cases(x[, ..required]), .(Participants = .N, Matched = sum(is.finite(match_km) & 
        match_km <= limit), Missing = sum(!is.finite(match_km)))], limit_km = limit, definition = "Nearest birthplace TDI reference within 1 km; full-cohort scaling and fixed primary covariates.")
}

if (sys.nframe() == 0L) {
    a <- commandArgs(TRUE)
    ids <- if (length(a)) 
        a
    else c("ED7", "ED8", "ED9", "birth", "matching")
    for (id in intersect(ids, c("ED7", "ED8", "ED9"))) {
        x <- switch(id, ED7 = ed7_fit(), ED8 = ed8_fit(), ED9 = ed9_fit())
        x$model <- fig2_model
        saveRDS(x, file.path(out, paste0(id, ".rds")))
    }
    s <- if (file.exists(file.path(out, "Supp.rds"))) 
        readRDS(file.path(out, "Supp.rds"))
    else list()
    if ("birth" %in% ids) 
        s$adjustment <- birth_fit()
    if ("matching" %in% ids) 
        s$matching <- tdi_match()
    saveRDS(s, file.path(out, "Supp.rds"), compress = "gzip")
}

