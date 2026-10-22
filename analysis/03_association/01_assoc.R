# 1. SES and health associations and strata
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

use(c("ukb", "geo"))

suppressPackageStartupMessages({
    library(survival)
    library(parallel)
})


# 2. Prepare complete analysis samples
ttwa_data <- function(covariates = cov1, cross_required = TRUE) {
    x <- dat()
    x[, `:=`(eid, as.character(eid))]
    if (testn > 0) 
        x <- x[seq_len(min(.N, testn))]
    ph <- ttwa_table()
    x <- merge(x, ph, by = "eid", all = FALSE, sort = FALSE)
    death <- as.data.table(readRDS(ukb_input("death")))
    death[, `:=`(eid, as.character(eid))]
    death[, `:=`(date_death2 = as.IDate(date_death), date_cvd = as.IDate(date_cvd_death))]
    x <- merge(x, death[, .(eid, date_death2, date_cvd)], by = "eid", all.x = TRUE, sort = FALSE)
    x[is.na(date_death), `:=`(date_death, date_death2)]
    x[, `:=`(date_death2, NULL)]
    x <- ses_pca(x)
    x[, `:=`(md_z = z(md), dtdi_z = z(dtdi), ses_z = z(ses), btdi_z = z(btdi), cross_ttwa = as.integer(cross_ttwa), 
        sex = factor(sex), center = factor(center), gen_batch = factor(gen_batch), smoke = factor(smoke), alcohol = factor(alcohol))]
    req <- c(covariates, "date_attend", if (cross_required) "cross_ttwa")
    x <- x[complete.cases(x[, req, with = FALSE])]
    x[, `:=`(c("sex", "center", "gen_batch", "smoke", "alcohol"), lapply(.SD, droplevels)), .SDcols = c("sex", 
        "center", "gen_batch", "smoke", "alcohol")]
    x
}

ukb_data <- function(covariates = cov1) {
    x <- ttwa_data(covariates)
    x <- x[complete.cases(x[, unique(c(ukb_var, covariates, "date_attend", "origin", "destination")), with = FALSE])]
    x[, `:=`(c("sex", "center", "gen_batch", "smoke", "alcohol"), lapply(.SD, droplevels)), .SDcols = c("sex", 
        "center", "gen_batch", "smoke", "alcohol")]
    x
}

fig2_data <- function() ttwa_data(cov1, FALSE)

make_icd <- function(ids) {
    code <- as.data.table(readRDS(ukb_input("icd10.code0")))
    date <- as.data.table(readRDS(ukb_input("icd10.date0")))
    code[, `:=`(eid, as.character(eid))]
    date[, `:=`(eid, as.character(eid))]
    ids <- as.character(ids)
    code <- melt(code[eid %in% ids], "eid", variable.name = "idx", value.name = "code", na.rm = TRUE)
    date <- melt(date[eid %in% ids], "eid", variable.name = "idx", value.name = "date", na.rm = TRUE)
    code[, `:=`(idx, sub("^i", "", idx))]
    date[, `:=`(idx, sub("^d", "", idx))]
    x <- merge(code, date, by = c("eid", "idx"))
    x[, `:=`(code, fmt(code))]
    unique(x[!is.na(code) & !is.na(date), .(eid, code, date = as.IDate(date))])
}

defs <- function(icd) {
    d <- list(`CAD/MI` = c("I20", "I21", "I22", "I23", "I24", "I25"), Stroke = c("I60", "I61", "I62", "I63", "I64"), 
        `Heart failure` = "I50", `Atrial fibrillation` = "I48", `Type 2 diabetes` = "E11", Depression = c("F32", 
            "F33"), `Anxiety disorder` = c("F40", "F41"), Dementia = c("F00", "F01", "F02", "F03", "G30"))
    c(list(`All-cause mortality` = NULL), lapply(d, function(p) icd[hit(code, p), .(event_date = min(date)), by = eid]))
}


# 3. Define endpoint-specific follow-up
survdat <- function(x, ed = NULL, mort = FALSE, cvd = FALSE) {
    d <- copy(x)
    if (mort && !cvd) {
        d[, `:=`(end0 = fifelse(!is.na(date_death) & date_death <= end, date_death, end), event = !is.na(date_death) & 
            date_death > date_attend & date_death <= end, prev = FALSE)]
    }
    else if (mort && cvd) {
        d[, `:=`(end0 = fifelse(!is.na(date_death) & date_death <= end, date_death, end), event = !is.na(date_cvd) & 
            date_cvd > date_attend & date_cvd <= end, prev = FALSE)]
    }
    else {
        d <- merge(d, ed, by = "eid", all.x = TRUE, sort = FALSE)
        d[, `:=`(censor, fifelse(!is.na(date_death) & date_death <= end, date_death, end))]
        d[, `:=`(event = !is.na(event_date) & event_date > date_attend & event_date <= censor, prev = !is.na(event_date) & 
            event_date <= date_attend, end0 = fifelse(!is.na(event_date) & event_date > date_attend & event_date <= 
            censor, event_date, censor))]
    }
    d[, `:=`(time, as.numeric(as.IDate(end0) - date_attend)/365.25)]
    d[!prev & is.finite(time) & time > 0]
}

ukb_surv <- function(x, od, outcome) {
    if (outcome == "All-cause mortality") 
        return(survdat(x, mort = TRUE))
    if (outcome == "CVD mortality") 
        return(survdat(x, mort = TRUE, cvd = TRUE))
    survdat(x, od[[outcome]])
}

ukb_coef <- function(m, outcome, type, exposures = names(ukb_var)) {
    sm <- summary(m)
    tab <- sm$coefficients
    pc <- grep("^Pr", colnames(tab), value = TRUE)[1]
    rbindlist(lapply(exposures, function(e) {
        v <- ukb_var[[e]]
        b <- coef(m)[v]
        se <- sqrt(vcov(m)[v, v])
        data.table(outcome, exposure = e, label = ukb_lab[[e]], model = type, n = if (type == "Cox") 
            m$n
        else nobs(m), event = if (type == "Cox") 
            m$nevent
        else NA_integer_, beta = b, se, estimate = if (type == "Cox") 
            exp(b)
        else b, lower = if (type == "Cox") 
            exp(b - 1.96 * se)
        else b - 1.96 * se, upper = if (type == "Cox") 
            exp(b + 1.96 * se)
        else b + 1.96 * se, p = tab[v, pc])
    }))
}

ukb_profile <- function(d) {
    x <- d[1]
    for (v in names(x)) {
        z0 <- d[[v]]
        if (is.factor(z0)) {
            u <- names(which.max(table(z0)))
            x[[v]] <- factor(u, levels = levels(z0))
        }
        else if (is.numeric(z0) || is.integer(z0)) {
            x[[v]] <- mean(z0, na.rm = TRUE)
        }
    }
    x
}


# 4. Estimate SES associations
fig2b_data <- function() {
    d <- fig2_data()[is.finite(ses_z)]
    d1 <- d[complete.cases(d[, c("md_z", "dtdi_z", fig2_cov0), with = FALSE])]
    d2 <- d[complete.cases(d[, c("cross_ttwa", "dtdi_z", cov0), with = FALSE])]
    m1 <- lm(as.formula(paste("ses_z ~", paste(c("md_z", "dtdi_z", fig2_cov0), collapse = " + "))), d1)
    m2 <- lm(as.formula(paste("ses_z ~", paste(c("cross_ttwa", "dtdi_z", cov0), collapse = " + "))), d2)
    a <- rbind(ukb_coef(m1, "SES", "Linear", "md"), ukb_coef(m1, "SES", "Linear", "dtdi"), ukb_coef(m2, "SES", 
        "Linear", "ttwa"))
    a[, `:=`(FDR = p.adjust(p, "BH"), filled = p.adjust(p, "BH") < 0.05)]
    q <- quantile(d1$md, c(0.01, 0.99), na.rm = TRUE)
    nd <- ukb_profile(d1)[rep(1, 101)]
    nd[, `:=`(md, seq(q[1], q[2], length.out = 101))]
    fm <- as.formula(paste("ses_z ~ splines::ns(md, df = 3) +", paste(c("dtdi_z", fig2_cov0), collapse = " + ")))
    pr <- predict(lm(fm, d1), nd, se.fit = TRUE)
    md <- data.table(x = nd$md, fit = as.numeric(pr$fit), lo = pr$fit - 1.96 * pr$se.fit, hi = pr$fit + 1.96 * 
        pr$se.fit)
    td <- copy(d1[!is.na(tdi4)])
    td[, `:=`(tdi4, factor(as.character(tdi4), levels = tdi_lev))]
    ft <- as.formula(paste("ses_z ~", paste(c("tdi4", "md_z", tdi_cov0), collapse = " + ")))
    mt <- lm(ft, td)
    nt <- ukb_profile(td)[rep(1, 4)]
    nt[, `:=`(tdi4, factor(tdi_lev, levels = tdi_lev))]
    pt <- predict(mt, nt, se.fit = TRUE)
    tdi <- data.table(group = tdi_lev, fit = as.numeric(pt$fit), lo = pt$fit - 1.96 * pt$se.fit, hi = pt$fit + 
        1.96 * pt$se.fit, origin = c("Low", "High", "Low", "High"), destination = c("Low", "Low", "High", "High"))
    comp <- rbindlist(lapply(c("edu.years", "income", "camsis"), function(v) {
        z1 <- copy(d1)
        z1[, `:=`(y, z(get(v)))]
        z2 <- copy(d2)
        z2[, `:=`(y, z(get(v)))]
        m1 <- lm(as.formula(paste("y ~", paste(c("md_z", "dtdi_z", fig2_cov0), collapse = " + "))), z1)
        m2 <- lm(as.formula(paste("y ~", paste(c("cross_ttwa", "dtdi_z", cov0), collapse = " + "))), z2)
        rbind(ukb_coef(m1, v, "Linear", "md"), ukb_coef(m1, v, "Linear", "dtdi"), ukb_coef(m2, v, "Linear", "ttwa"))[, 
            .(component = v, exposure, estimate, lower, upper, p)]
    }))
    comp[, `:=`(FDR = p.adjust(p, "BH"), filled = p.adjust(p, "BH") < 0.05)]
    c(list(main = a, md = md, tdi = tdi, comp = comp), fig2b_direction())
}

fig2b_direction <- function(x = direction_data()) {
    d <- x[is.finite(ses_z)]
    fit <- function(v) lm(reformulate(c(v, "md_z", "dtdi_z", cov0), "ses_z"), d)
    m <- fit("urban_route")
    nd <- ukb_profile(d)[rep(1, 4)]
    nd[, `:=`(urban_route, factor(levels(d$urban_route), levels = levels(d$urban_route)))]
    pr <- predict(m, nd, se.fit = TRUE)
    urban <- data.table(group = levels(d$urban_route), fit = as.numeric(pr$fit), lo = pr$fit - 1.96 * pr$se.fit, 
        hi = pr$fit + 1.96 * pr$se.fit, origin = c("Urban", "Urban", "Rural", "Rural"), destination = c("Urban", 
            "Rural", "Urban", "Rural"), n = nobs(m))
    co <- summary(m)$coefficients
    i <- startsWith(rownames(co), "urban_route")
    routes <- data.table(outcome = "SES", term = rownames(co)[i], n = nobs(m), beta = co[i, "Estimate"], se = co[i, 
        "Std. Error"], p = co[i, "Pr(>|t|)"])
    routes[, `:=`(lower = beta - 1.96 * se, upper = beta + 1.96 * se, FDR = p.adjust(p, "BH"))]
    m <- fit("density_z")
    b <- unname(coef(m)["density_z"])
    se <- sqrt(vcov(m)["density_z", "density_z"])
    density <- data.table(outcome = "SES", exposure = "density", label = "Density change", model = "Linear", n = nobs(m), 
        event = NA_integer_, beta = b, se, estimate = b, lower = b - 1.96 * se, upper = b + 1.96 * se, p = summary(m)$coefficients["density_z", 
            "Pr(>|t|)"])
    density[, `:=`(FDR = p, filled = p < 0.05)]
    list(urban = urban, density = density, routes = routes)
}


# 5. Estimate health associations
fig2c_res <- function() {
    pair <- data.table(term = c("Low to high vs High to low", "High to high vs Low to low", "High to low vs High to high", 
        "Low to high vs Low to low"), a = c("Low to high", "High to high", "High to low", "Low to high"), b = c("High to low", 
        "Low to low", "High to high", "Low to low"))
    x <- fig2_data()
    x[, `:=`(md_q = cut(md, unique(quantile(md, seq(0, 1, 0.25), na.rm = TRUE)), include.lowest = TRUE, labels = paste0("Q", 
        1:4)), tdi4f = factor(as.character(tdi4), levels = c("High to low", "Low to low", "High to high", "Low to high")))]
    icd <- make_icd(x$eid)
    od <- defs(icd)
    r <- mclapply(forest_out, function(o) {
        d <- ukb_surv(x, od, o)
        d1 <- d[complete.cases(d[, c("md_z", "dtdi_z", fig2_cov1, "tdi4f"), with = FALSE])]
        d2 <- d[complete.cases(d[, c("cross_ttwa", "dtdi_z", cov1), with = FALSE])]
        m1 <- coxph(as.formula(paste("Surv(time,event) ~", paste(c("md_z", "dtdi_z", fig2_cov1), collapse = " + "))), 
            d1, ties = "efron")
        m2 <- coxph(as.formula(paste("Surv(time,event) ~", paste(c("cross_ttwa", "dtdi_z", cov1), collapse = " + "))), 
            d2, ties = "efron")
        z1 <- ukb_coef(m1, o, "Cox", c("md", "dtdi"))
        z2 <- ukb_coef(m2, o, "Cox", "ttwa")
        main <- rbind(z1[exposure == "md"], z2)[, .(outcome, term = unname(ukb_lab[exposure]), n, event, HR = estimate, 
            lo = lower, hi = upper, p)]
        mq <- coxph(as.formula(paste("Surv(time,event) ~", paste(c("md_q", "dtdi_z", fig2_cov1), collapse = " + "))), 
            d1, ties = "efron")
        quartile <- rbind(data.table(outcome = o, term = "MD Q1", n = mq$n, event = mq$nevent, HR = 1, lo = NA_real_, 
            hi = NA_real_, p = NA_real_), rbindlist(lapply(2:4, function(i) {
            v <- paste0("md_qQ", i)
            b <- coef(mq)[v]
            se <- sqrt(vcov(mq)[v, v])
            data.table(outcome = o, term = paste0("MD Q", i), n = mq$n, event = mq$nevent, HR = exp(b), lo = exp(b - 
                1.96 * se), hi = exp(b + 1.96 * se), p = summary(mq)$coefficients[v, "Pr(>|z|)"])
        })))
        main <- rbind(quartile, main, fill = TRUE)
        cont <- ukb_coef(m1, o, "Cox", "dtdi")[, .(outcome, term = paste("Per SD", dtdi_lab), n, event, HR = estimate, 
            lo = lower, hi = upper, p)]
        mg <- coxph(as.formula(paste("Surv(time,event) ~", paste(c("tdi4f", "md_z", tdi_cov1), collapse = " + "))), 
            d1, ties = "efron")
        comp <- rbindlist(lapply(seq_len(nrow(pair)), function(i) fig2d_contrast(mg, pair$a[i], pair$b[i], pair$term[i], 
            o, nrow(d1), sum(d1$event))))
        list(main = main, tdi = rbind(cont, comp, fill = TRUE))
    }, mc.cores = core)
    stopifnot(all(lengths(r) > 0), !any(vapply(r, inherits, logical(1), "try-error")))
    ans <- list(main = rbindlist(lapply(r, `[[`, "main")), tdi = rbindlist(lapply(r, `[[`, "tdi")))
    for (nm in names(ans)) {
        ans[[nm]][, `:=`(FDR, NA_real_)]
        ans[[nm]][is.finite(p), `:=`(FDR, p.adjust(p, "BH"))]
        ans[[nm]][, `:=`(sig, is.finite(FDR) & FDR < 0.05)]
    }
    ans
}

fig2d_contrast <- function(m, a, b, label, outcome, n, event) {
    cf <- coef(m)
    vc <- vcov(m)
    v <- setNames(rep(0, length(cf)), names(cf))
    ref <- levels(model.frame(m)$tdi4f)[1]
    trm <- function(g) if (g == ref) 
        NA_character_
    else paste0("tdi4f", g)
    ta <- trm(a)
    tb <- trm(b)
    if (!is.na(ta) && ta %in% names(v)) 
        v[ta] <- 1
    if (!is.na(tb) && tb %in% names(v)) 
        v[tb] <- v[tb] - 1
    beta <- sum(v * cf)
    se <- sqrt(as.numeric(t(v) %*% vc %*% v))
    data.table(outcome, term = label, n, event, model = "cox", HR = exp(beta), lo = exp(beta - 1.96 * se), hi = exp(beta + 
        1.96 * se), p = 2 * pnorm(-abs(beta/se)))
}


# 6. Analyze migration directions
direction_data <- function() {
    file <- file.path(root, "data/ref/population_2011.asc")
    header <- strsplit(trimws(readLines(file, n = 6)), "[[:space:]]+")
    h <- setNames(vapply(header, function(x) as.numeric(x[2]), 0), tolower(vapply(header, `[`, "", 1)))
    g <- matrix(scan(file, skip = 6, quiet = TRUE), nrow = h["nrows"], byrow = TRUE)
    density <- function(east, north) {
        col <- floor((east - h["xllcorner"])/h["cellsize"]) + 1L
        row <- h["nrows"] - floor((north - h["yllcorner"])/h["cellsize"])
        ok <- is.finite(row) & is.finite(col) & row >= 1 & row <= nrow(g) & col >= 1 & col <= ncol(g)
        x <- rep(NA_real_, length(east))
        x[ok] <- g[cbind(row[ok], col[ok])]
        x[x < 0 | x == h["nodata_value"]] <- NA_real_
        x/(h["cellsize"]/1000)^2
    }
    u <- as.data.table(readRDS(file.path(root, "data/ref/urban/ukb.rds")))
    x <- merge(fig2_data(), u[, .(eid, urban_route)], by = "eid", sort = FALSE)
    x[, `:=`(density_shift, log1p(density(home_east, home_north)) - log1p(density(birth_east, birth_north)))]
    need <- c("density_shift", "urban_route", "md_z", "dtdi_z", cov1)
    x <- x[complete.cases(x[, ..need]) & is.finite(density_shift)]
    x[, `:=`(density_z = z(density_shift), urban_route = factor(urban_route, levels = c("Urban to Urban", "Urban to Rural", 
        "Rural to Urban", "Rural to Rural")))]
    x
}

direction_fit <- function() {
    x <- direction_data()
    od <- defs(make_icd(x$eid))
    results <- rbindlist(lapply(forest_out, function(outcome) {
        d <- ukb_surv(x, od, outcome)
        rbindlist(lapply(c("density_z", "urban_route"), function(v) {
            f <- as.formula(paste("Surv(time,event) ~", paste(c(v, "md_z", "dtdi_z", cov1), collapse = " + ")))
            m <- coxph(f, d, ties = "efron")
            s <- summary(m)$coefficients
            i <- startsWith(rownames(s), v)
            data.table(outcome, term = rownames(s)[i], n = nrow(d), event = sum(d$event), beta = s[i, "coef"], 
                se = s[i, "se(coef)"], HR = exp(s[i, "coef"]), lo = exp(s[i, "coef"] - 1.96 * s[i, "se(coef)"]), 
                hi = exp(s[i, "coef"] + 1.96 * s[i, "se(coef)"]), p = s[i, "Pr(>|z|)"])
        }))
    }))
    results[, `:=`(FDR = p.adjust(p, "BH"), sig = p.adjust(p, "BH") < 0.05)]
    list(results = results, counts = x[, .(N = .N), urban_route], density_sd = sd(x$density_shift), definition = "log(1 + recruitment population/km2) - log(1 + birthplace population/km2)", 
        reference = "Urban to Urban", covariates = c("md_z", "dtdi_z", cov1))
}

direction_table <- function(r) {
    ref <- unique(r$results[, .(outcome, n, event)])
    ref[, `:=`(term = "urban_reference", HR = 1, lo = NA_real_, hi = NA_real_, p = NA_real_, FDR = NA_real_, sig = FALSE)]
    rbindlist(list(r$results, ref), fill = TRUE)
}


# 7. Fit prespecified strata
strat_spec <- function() list(tdi = list(title = "MD stratified by TDI group", var = "tdi_g", levels = tdi_lev, 
    colors = c("#79295E", "#A33D78", "#C65D94", "#DF8DB4"), drop = character()), work = list(title = "MD stratified by occupational physical demand", 
    var = "work_g", levels = c("Never/rarely", "Sometimes", "Usually", "Always"), colors = c("#176873", "#3D8F99", 
        "#72AFB7", "#A2CCD0"), drop = character()))

strat_fit <- function(d, outcome, id, s) {
    v <- s$var
    lev <- s$levels
    cv <- setdiff(fig2_cov1, s$drop)
    z <- copy(d[!is.na(get(v))])
    z[, `:=`(group0, factor(get(v), levels = lev))]
    z <- z[complete.cases(z[, c("md_z", "dtdi_z", "group0", cv), with = FALSE])]
    m0 <- try(coxph(as.formula(paste("Surv(time,event) ~ md_z + dtdi_z + group0 +", paste(cv, collapse = " + "))), 
        z, ties = "efron"), silent = TRUE)
    mi <- try(coxph(as.formula(paste("Surv(time,event) ~ md_z * group0 + dtdi_z +", paste(cv, collapse = " + "))), 
        z, ties = "efron"), silent = TRUE)
    pint <- NA_real_
    if (!inherits(m0, "try-error") && !inherits(mi, "try-error")) {
        q <- anova(m0, mi, test = "LRT")
        pc <- grep("Pr|P\\(", names(q), value = TRUE)
        if (length(pc)) 
            pint <- as.numeric(q[2, pc[1]])
    }
    rbindlist(lapply(lev, function(g) {
        a <- z[group0 == g]
        m <- try(coxph(as.formula(paste("Surv(time,event) ~ md_z + dtdi_z +", paste(cv, collapse = " + "))), a, 
            ties = "efron"), silent = TRUE)
        if (inherits(m, "try-error")) 
            return(data.table(strat = id, panel = s$title, outcome, group = g, n = nrow(a), event = sum(a$event), 
                HR = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, p_interaction = pint))
        q <- summary(m)
        data.table(strat = id, panel = s$title, outcome, group = g, n = nrow(a), event = sum(a$event), HR = q$conf.int["md_z", 
            "exp(coef)"], lo = q$conf.int["md_z", "lower .95"], hi = q$conf.int["md_z", "upper .95"], p = q$coefficients["md_z", 
            "Pr(>|z|)"], p_interaction = pint)
    }))
}

fig2c_strat_res <- function() {
    x <- fig2_data()
    x[, `:=`(tdi_g = as.character(tdi4), work_g = fcase(work_heavy == 1, "Never/rarely", work_heavy == 2, "Sometimes", 
        work_heavy == 3, "Usually", work_heavy == 4, "Always"))]
    sp <- strat_spec()
    icd <- make_icd(x$eid)
    od <- defs(icd)
    r <- rbindlist(mclapply(forest_out, function(o) {
        d <- survdat(x, od[[o]], o %in% c("All-cause mortality", "CVD mortality"), o == "CVD mortality")
        rbindlist(lapply(names(sp), function(id) strat_fit(d, o, id, sp[[id]])))
    }, mc.cores = core), fill = TRUE)
    r[, `:=`(FDR, p.adjust(p, "BH")), by = strat]
    q <- unique(r[, .(strat, outcome, p_interaction)])
    q[, `:=`(interaction_FDR, p.adjust(p_interaction, "BH")), by = strat]
    merge(r, q, by = c("strat", "outcome", "p_interaction"), all.x = TRUE, sort = FALSE)
}


# 8. Save aggregate results
fig2_results <- function(part = NULL, fresh = FALSE) {
    file <- file.path(out, "Fig2.rds")
    r <- if (!stale(file, ukb_inputs())) 
        readRDS(file)
    else list()
    if (!identical(r$model, fig2_model)) 
        r <- list(model = fig2_model)
    build <- list(ses = fig2b_data, forest = fig2c_res, strata = fig2c_strat_res, direction = direction_fit)
    for (nm in if (is.null(part)) 
        names(build)
    else part) {
        if (fresh || is.null(r[[nm]])) {
            message("UKB: ", nm)
            r[[nm]] <- build[[nm]]()
            saveRDS(r, file, compress = "gzip")
        }
    }
    if (is.null(part) || part == "ses") {
        if (is.null(r$ses$urban)) {
            r$ses <- c(r$ses, fig2b_direction())
            saveRDS(r, file, compress = "gzip")
        }
    }
    if (is.null(part)) 
        r
    else r[[part]]
}

if (sys.nframe() == 0L) {
    a <- commandArgs(TRUE)
    fig2_results(if (length(a)) 
        a[1]
    else NULL, fresh = TRUE)
}

