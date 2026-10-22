# 1. Interaction tests and curves
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

use("assoc")

# 2. Fit nested models
int_fit <- function(d, outcome, type, name, base, full) {
    lhs <- if (type == "Cox") 
        "Surv(time,event)"
    else "ses_z"
    fit <- function(rhs) {
        f <- as.formula(paste(lhs, "~", paste(rhs, collapse = " + ")))
        if (type == "Cox") 
            coxph(f, d, ties = "efron", model = TRUE)
        else lm(f, d)
    }
    m0 <- fit(base)
    m1 <- fit(full)
    a <- anova(m0, m1, test = if (type == "Cox") 
        "LRT"
    else "F")
    pc <- grep("Pr|P\\(", names(a), value = TRUE)[1]
    nm <- names(coef(m1))[grepl(":", names(coef(m1)), fixed = TRUE)]
    b <- coef(m1)[nm]
    se <- sqrt(diag(vcov(m1)))[nm]
    ip <- if (type == "Cox") 
        2 * pnorm(-abs(b/se))
    else 2 * stats::pt(-abs(b/se), df = df.residual(m1))
    term <- data.table(outcome, interaction = name, model = type, term = nm, beta = unname(b), se = unname(se), 
        ratio = if (type == "Cox") 
            exp(unname(b))
        else NA_real_, lower = if (type == "Cox") 
            exp(b - 1.96 * se)
        else b - 1.96 * se, upper = if (type == "Cox") 
            exp(b + 1.96 * se)
        else b + 1.96 * se, p = unname(ip))
    list(test = data.table(outcome, interaction = name, model = type, df = a$Df[2], n = if (type == "Cox") m1$n else nobs(m1), 
        event = if (type == "Cox") m1$nevent else NA_integer_, interaction_beta = if (length(b) == 1) unname(b) else NA_real_, 
        interaction_ratio = if (type == "Cox" && length(b) == 1) exp(unname(b)) else NA_real_, p = as.numeric(a[[pc]][2])), 
        term = term, model = m1)
}

int_effect <- function(m, outcome, name, type, exposure, modifier) {
    mf <- model.frame(m)
    mf <- mf[, -attr(terms(m), "response"), drop = FALSE]
    ref <- as.data.frame(lapply(mf, function(v) {
        if (is.factor(v)) 
            factor(levels(v)[1], levels = levels(v))
        else 0
    }))
    rbindlist(lapply(seq_len(nrow(modifier)), function(i) {
        a <- ref
        for (v in intersect(names(modifier), names(a))) {
            a[[v]] <- if (is.factor(mf[[v]])) 
                factor(modifier[[v]][i], levels = levels(mf[[v]]))
            else modifier[[v]][i]
        }
        a0 <- a1 <- a
        a0[[exposure]] <- 0
        a1[[exposure]] <- 1
        tr <- delete.response(terms(m))
        x0 <- model.matrix(tr, a0, contrasts.arg = m$contrasts)
        x1 <- model.matrix(tr, a1, contrasts.arg = m$contrasts)
        b <- coef(m)
        w <- x1[, names(b), drop = FALSE] - x0[, names(b), drop = FALSE]
        est <- drop(w %*% b)
        se <- sqrt(drop(w %*% vcov(m) %*% t(w)))
        data.table(outcome, interaction = name, model = type, exposure, level = modifier$level[i], order = modifier$order[i], 
            beta = est, se, estimate = if (type == "Cox") 
                exp(est)
            else est, lower = if (type == "Cox") 
                exp(est - 1.96 * se)
            else est - 1.96 * se, upper = if (type == "Cox") 
                exp(est + 1.96 * se)
            else est + 1.96 * se)
    }))
}

int_run <- function(d, outcome, type, name, base, full, exposure, modifier) {
    z0 <- int_fit(d, outcome, type, name, base, full)
    list(test = z0$test, term = z0$term, effect = int_effect(z0$model, outcome, name, type, exposure, modifier))
}


# 3. Test health interactions
int_health <- function(x, od, outcome) {
    d <- ukb_surv(x, od, outcome)
    cv <- fig2_cov1
    age <- data.table(age_z = predict(lm(age_z ~ age, d), newdata = data.frame(age = seq(40, 70, 10))), level = as.character(seq(40, 
        70, 10)), order = 1:4)
    sex <- data.table(sex2 = c("Female", "Male"), level = c("Women", "Men"), order = 1:2)
    tdi <- data.table(tdi4 = tdi_lev, level = tdi_lev, order = 1:4)
    work <- data.table(work_g = c("Never/rarely", "Sometimes", "Usually", "Always"), level = c("Never/rarely", 
        "Sometimes", "Usually", "Always"), order = 1:4)
    core <- c("md_z", "dtdi_z", cv)
    z0 <- list(int_run(d, outcome, "Cox", "MD x dTDI", core, c("md_z * dtdi_z", cv), "md_z", data.table(dtdi_z = c(-1, 
        0, 1), level = c("-1 SD", "0", "+1 SD"), order = 1:3)), int_run(d, outcome, "Cox", "MD x age", c("md_z", 
        "age_z", "dtdi_z", setdiff(cv, "age")), c("md_z * age_z", "dtdi_z", setdiff(cv, "age")), "md_z", age), 
        int_run(d, outcome, "Cox", "dTDI x age", c("dtdi_z", "age_z", "md_z", setdiff(tdi_cov1, "age")), c("dtdi_z * age_z", 
            "md_z", setdiff(tdi_cov1, "age")), "dtdi_z", age), int_run(d, outcome, "Cox", "MD x sex", c("md_z", 
            "sex2", "dtdi_z", setdiff(cv, "sex")), c("md_z * sex2", "dtdi_z", setdiff(cv, "sex")), "md_z", sex), 
        int_run(d, outcome, "Cox", "dTDI x sex", c("dtdi_z", "sex2", "md_z", setdiff(tdi_cov1, "sex")), c("dtdi_z * sex2", 
            "md_z", setdiff(tdi_cov1, "sex")), "dtdi_z", sex), int_run(d, outcome, "Cox", "MD x TDI group", c("md_z", 
            "tdi4", cv), c("md_z * tdi4", cv), "md_z", tdi))
    dw <- copy(d)[, `:=`(work_g, factor(fcase(work_heavy == 1, "Never/rarely", work_heavy == 2, "Sometimes", work_heavy == 
        3, "Usually", work_heavy == 4, "Always"), levels = work$level))][!is.na(work_g)]
    z0[[length(z0) + 1]] <- int_run(dw, outcome, "Cox", "MD x physical demand", c("md_z", "work_g", "dtdi_z", cv), 
        c("md_z * work_g", "dtdi_z", cv), "md_z", work)
    list(test = rbindlist(lapply(z0, `[[`, "test")), term = rbindlist(lapply(z0, `[[`, "term")), effect = rbindlist(lapply(z0, 
        `[[`, "effect")))
}


# 4. Test SES interactions
int_ses <- function(x) {
    d <- x[is.finite(ses_z)]
    cv <- fig2_cov0
    core <- c("md_z", "dtdi_z", cv)
    z0 <- list(int_run(d, "SES", "Linear", "MD x dTDI", core, c("md_z * dtdi_z", cv), "md_z", data.table(dtdi_z = c(-1, 
        0, 1), level = c("-1 SD", "0", "+1 SD"), order = 1:3)))
    list(test = rbindlist(lapply(z0, `[[`, "test")), term = rbindlist(lapply(z0, `[[`, "term")), effect = rbindlist(lapply(z0, 
        `[[`, "effect")))
}

int_direction <- function(x = direction_data()) {
    q <- quantile(x$density_shift, seq(0, 1, 0.25))
    x[, `:=`(density_q, cut(density_shift, q, include.lowest = TRUE, labels = paste0("Q", 1:4)))]
    density <- x[, .(density_z = median(density_z)), by = density_q][order(density_q)]
    od <- defs(make_icd(x$eid))
    parts <- mclapply(c("SES", forest_out), function(o) {
        type <- if (o == "SES") 
            "Linear"
        else "Cox"
        cv <- if (o == "SES") 
            fig2_cov0
        else fig2_cov1
        d <- if (o == "SES") 
            copy(x[is.finite(ses_z)])
        else ukb_surv(x, od, o)
        r <- lapply(c("urban_route", "density_z"), function(v) {
            name <- if (v == "urban_route") 
                "MD x rural-urban route"
            else "MD x density change"
            modifier <- if (v == "urban_route") 
                data.table(urban_route = levels(d$urban_route), level = levels(d$urban_route), order = 1:4)
            else density[, .(density_z, level = as.character(density_q), order = 1:4)]
            int_run(d, o, type, name, c("md_z", v, "dtdi_z", cv), c(paste0("md_z * ", v), "dtdi_z", cv), "md_z", 
                modifier)
        })
        list(Tests = rbindlist(lapply(r, `[[`, "test")), Terms = rbindlist(lapply(r, `[[`, "term")), Effects = rbindlist(lapply(r, 
            `[[`, "effect")))
    }, mc.cores = core)
    stopifnot(all(lengths(parts) > 0), !any(vapply(parts, inherits, logical(1), "try-error")))
    setNames(lapply(c("Tests", "Terms", "Effects"), function(v) rbindlist(lapply(parts, `[[`, v))), c("Tests", 
        "Terms", "Effects"))
}


# 5. Correct the testing family
int_adjust <- function(r) {
    r$Tests[, `:=`(FDR, p.adjust(p, "BH"))]
    r$Terms[, `:=`(FDR, p.adjust(p, "BH"))]
    drop <- intersect(c("p", "FDR", "p_interaction", "q_interaction"), names(r$Effects))
    if (length(drop)) 
        r$Effects[, `:=`((drop), NULL)]
    r$Effects <- merge(r$Effects, r$Tests[, .(outcome, interaction, p, FDR)], by = c("outcome", "interaction"), 
        all.x = TRUE, sort = FALSE)
    r
}

int_main <- function() data.table(interaction = c("MD x rural-urban route", "dTDI x age", "MD x dTDI", "MD x physical demand", 
    "dTDI x age", "dTDI x sex", "MD x age", "MD x sex"), outcome = c("SES", "CVD mortality", "SES", "Anxiety disorder", 
    "All-cause mortality", "All-cause mortality", "Anxiety disorder", "CAD/MI"), type = c("urban", "age", "ses", 
    "work", "age", "sex", "age", "sex"))


# 6. Save estimates and curves
run_interaction <- function() {
    x <- fig2_data()
    x[, `:=`(age_z = z(age), sex2 = factor(fifelse(as.character(sex) %in% c("1", "Male", "M"), "Male", "Female"), 
        levels = c("Female", "Male")))]
    icd <- make_icd(x$eid)
    od <- defs(icd)
    z0 <- mclapply(forest_out, function(o) {
        message("interaction ", o)
        int_health(x, od, o)
    }, mc.cores = core)
    stopifnot(all(lengths(z0) > 0), !any(vapply(z0, inherits, logical(1), "try-error")))
    z0 <- c(z0, list(int_ses(x)))
    test <- rbindlist(lapply(z0, `[[`, "test"))
    term <- rbindlist(lapply(z0, `[[`, "term"))
    effect <- rbindlist(lapply(z0, `[[`, "effect"))
    result <- list(Tests = test, Terms = term, Effects = effect)
    extra <- int_direction()
    result <- int_adjust(setNames(lapply(names(result), function(v) rbindlist(list(result[[v]], extra[[v]]), fill = TRUE)), 
        names(result)))
    if (testn == 0) {
        file <- file.path(out, "Fig2.rds")
        cache <- if (!stale(file, ukb_inputs())) 
            readRDS(file)
        else list()
        if (!identical(cache$model, fig2_model)) 
            cache <- list(model = fig2_model)
        cache$interaction <- result
        saveRDS(cache, file, compress = "gzip")
    }
    invisible(result)
}

int_bounds <- function() {
    x <- fig2_data()
    rbindlist(lapply(c("md_z", "dtdi_z"), function(v) data.table(exposure = v, lower = quantile(x[[v]], 0.01, na.rm = TRUE), 
        upper = quantile(x[[v]], 0.99, na.rm = TRUE), minimum = min(x[[v]], na.rm = TRUE), maximum = max(x[[v]], 
            na.rm = TRUE))))
}

interaction_data <- function() {
    file <- file.path(out, "Fig2.rds")
    cache <- if (!stale(file, ukb_inputs())) 
        readRDS(file)
    else list()
    r <- if (identical(cache$model, fig2_model)) 
        cache$interaction
    else NULL
    if (is.null(r)) 
        return(run_interaction())
    if (!"MD x density change" %chin% r$Tests$interaction) {
        extra <- int_direction()
        r <- int_adjust(setNames(lapply(names(r), function(v) rbindlist(list(r[[v]], extra[[v]]), fill = TRUE)), 
            names(r)))
        cache$interaction <- r
        saveRDS(cache, file, compress = "gzip")
    }
    r
}

if (sys.nframe() == 0L) run_interaction()

