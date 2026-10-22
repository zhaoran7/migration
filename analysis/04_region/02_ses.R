# 1. Regional SES heterogeneity
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

use(c("assoc", "region"))

suppressPackageStartupMessages({
    library(lme4)
})


# 2. Fit regional SES models
supp_region_fit <- function() {
    x <- fig2_data()
    cv <- c("md_z", "dtdi_z", setdiff(fig2_cov0, "btdi_z"))
    x <- x[is.finite(ses_z) & complete.cases(x[, c(cv, "origin", "destination"), with = FALSE])]
    x[, `:=`(origin = factor(origin), destination = factor(destination))]
    for (v in names(x)) if (is.factor(x[[v]])) 
        set(x, j = v, value = droplevels(x[[v]]))
    fit <- function(random) {
        m <- lme4::lmer(reformulate(c(cv, random), "ses_z"), x, REML = FALSE, control = lme4::lmerControl(optimizer = "bobyqa", 
            optCtrl = list(maxfun = 1e+05)))
        stopifnot(m@optinfo$conv$opt == 0L, is.null(m@optinfo$conv$lme4$messages))
        m
    }
    m <- fit(c("(1 | origin)", "(1 | destination)"))
    vc <- as.data.table(lme4::VarCorr(m))
    variance <- setNames(vc$vcov, vc$grp)
    areas <- rbindlist(lapply(c("origin", "destination"), function(v) {
        a <- lme4::ranef(m, condVar = TRUE)[[v]]
        data.table(context = v, code = rownames(a), beta = a[, 1], se = sqrt(attr(a, "postVar")[1, 1, ]))
    }))
    count <- rbindlist(lapply(c("origin", "destination"), function(v) x[, .(n = .N), by = .(code = as.character(get(v)))][, 
        `:=`(context, v)]))
    areas <- merge(areas, count, by = c("context", "code"))
    areas[, `:=`(lower = beta - 1.96 * se, upper = beta + 1.96 * se)]
    tests <- rbindlist(lapply(c("origin", "destination"), function(v) {
        reduced <- fit(paste0("(1 | ", setdiff(c("origin", "destination"), v), ")"))
        lr <- 2 * as.numeric(logLik(m) - logLik(reduced))
        stopifnot(lr > -0.05)
        data.table(context = v, n = nrow(x), variance = variance[v], SD = sqrt(variance[v]), VPC = variance[v]/sum(variance), 
            LR = lr, p = pchisq(max(0, lr), 1, lower.tail = FALSE))
    }))
    tests[, `:=`(q, p.adjust(p, "BH"))]
    slopes <- lapply(c("md_z", "dtdi_z"), function(v) {
        x[, `:=`(slope = get(v), origin_s = origin, destination_s = destination)]
        z <- fit(c("(1 | origin)", "(1 | destination)", "(0 + slope | origin_s)", "(0 + slope | destination_s)"))
        s <- as.data.table(lme4::VarCorr(z))[grp %chin% c("origin_s", "destination_s")]
        lr <- 2 * as.numeric(logLik(z) - logLik(m))
        stopifnot(lr > -0.05)
        metrics <- data.table(outcome = "SES", exposure = v, context = sub("_s$", "", s$grp), slope_variance = s$vcov, 
            slope_sd = sqrt(s$vcov), mean_beta = unname(lme4::fixef(z)[v]))
        metrics[, `:=`(beta_low = mean_beta - 1.96 * slope_sd, beta_high = mean_beta + 1.96 * slope_sd)]
        list(test = data.table(outcome = "SES", exposure = v, model = "Linear", n = nrow(x), events = NA_integer_, 
            LR_raw = lr, df = 2L, p = pchisq(max(0, lr), 2, lower.tail = FALSE)), metrics = metrics, fixed = lme4::fixef(z), 
            variance = lme4::VarCorr(z), formula = deparse(formula(z)))
    })
    names(slopes) <- c("md_z", "dtdi_z")
    list(tests = tests, areas = areas, slopes = slopes, N = nrow(x), fixed = lme4::fixef(m), variance = variance, 
        formula = deparse(formula(m)), covariates = cv)
}


# 3. Collect regional estimates
supp_region_table <- function(r) {
    label <- c(origin = "Birthplace", destination = "Recruitment")
    first <- r$tests[, .(Analysis = "Regional variance", Exposure = NA_character_, Location = unname(label[context]), 
        Region = NA_character_, Participants = n, `Regional variance` = variance, `Regional SD` = SD, `Variance partition (%)` = 100 * 
            VPC, `Likelihood-ratio statistic` = LR, `Degrees of freedom` = 1L, p, q)]
    areas <- r$areas[, .(Analysis = "Conditional regional deviation", Location = unname(label[context]), Region = code, 
        Participants = n, Beta = beta, `Conditional 95% CI lower` = lower, `Conditional 95% CI upper` = upper)]
    tests <- rbindlist(lapply(r$slopes, `[[`, "test"))[, `:=`(q, p.adjust(p, "BH"))]
    last <- tests[, .(Analysis = "Regional slope test", Exposure = fifelse(exposure == "md_z", "MD", "Delta TDI"), 
        Participants = n, `Likelihood-ratio statistic` = LR_raw, `Degrees of freedom` = df, p, q)]
    ranges <- rbindlist(lapply(r$slopes, `[[`, "metrics"))[, .(Analysis = "Between-area slope distribution", Exposure = fifelse(exposure == 
        "md_z", "MD", "Delta TDI"), Location = unname(label[context]), Participants = r$N, `Regional variance` = slope_variance, 
        `Regional SD` = slope_sd, Beta = mean_beta, `Between-area 95% range lower` = beta_low, `Between-area 95% range upper` = beta_high)]
    d <- rbindlist(list(first, areas, last, ranges), fill = TRUE)
    for (v in c("p", "q")) set(d, j = v, value = ifelse(is.na(d[[v]]), NA_character_, ifelse(d[[v]] == 0, "<1e-300", 
        formatC(d[[v]], format = "e", digits = 3))))
    d
}

if (sys.nframe() == 0L) {
    supp_region_fit()
}

