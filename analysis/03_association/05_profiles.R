# 1. Twenty-four joint migration profiles
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

use(c("assoc", "urban"))


# 2. Define outcome-independent profiles
profile_cluster <- function(x) {
    x <- copy(x)
    x[, `:=`(stratum, interaction(urban_route, cross_ttwa, drop = TRUE))]
    Z <- scale(cbind(log1p(x$md), x$dtdi, x$density_shift))
    set.seed(20261020)
    offset <- 0L
    models <- list()
    x[, `:=`(profile, NA_character_)]
    for (st in levels(x$stratum)) {
        ii <- which(x$stratum == st)
        m <- kmeans(Z[ii, , drop = FALSE], centers = 3, nstart = 30, algorithm = "Lloyd", iter.max = 200)
        ord <- order(m$centers[, 1])
        ids <- sprintf("G%02d", offset + 1:3)
        x[ii, `:=`(profile, ids[match(m$cluster, ord)])]
        models[[st]] <- list(centres = m$centers[ord, , drop = FALSE], ids = ids)
        offset <- offset + 3L
    }
    stopifnot(offset == 24L, !anyNA(x$profile))
    list(data = x, clusters = models)
}


# 3. Link geographic context
profile_geography <- function(x) {
    locations <- as.data.table(readRDS(file.path(phe, "all.rds")))[, .(eid = as.character(eid), birth_east, birth_north, 
        home_east, home_north)]
    g <- locations[match(x$eid, eid)]
    areas <- urban_polygons()
    b <- urban_assign(g$birth_east, g$birth_north, areas)
    h <- urban_assign(g$home_east, g$home_north, areas, grid = TRUE)
    file <- file.path(root, "data/ref/population_2011.asc")
    a <- strsplit(trimws(readLines(file, n = 6)), "[[:space:]]+")
    header <- setNames(vapply(a, function(v) as.numeric(v[2]), 0), tolower(vapply(a, "[", "", 1)))
    grid <- matrix(scan(file, skip = 6, quiet = TRUE), nrow = header["nrows"], byrow = TRUE)
    density <- function(e, n) {
        a <- sapply(seq_len(9), function(i) {
            offsets <- CJ(e = -1:1, n = -1:1)
            col <- floor((e - header["xllcorner"])/header["cellsize"]) + 1L + offsets$e[i]
            row <- header["nrows"] - floor((n - header["yllcorner"])/header["cellsize"]) - offsets$n[i]
            ok <- is.finite(row) & is.finite(col) & row >= 1 & row <= nrow(grid) & col >= 1 & col <= ncol(grid)
            v <- rep(NA_real_, length(e))
            v[ok] <- grid[cbind(row[ok], col[ok])]
            v[v < 0 | v == header["nodata_value"]] <- NA_real_
            v/(header["cellsize"]/1000)^2
        })
        rowMeans(a, na.rm = TRUE)
    }
    x[, `:=`(urban_route = factor(paste(b, "to", h), levels = c("Urban to Urban", "Urban to Rural", "Rural to Urban", 
        "Rural to Rural")), density_shift = log1p(density(g$home_east, g$home_north)) - log1p(density(g$birth_east, 
        g$birth_north)))]
    stopifnot(!anyNA(x$urban_route), all(is.finite(x$density_shift)))
    x[, `:=`(density_z, z(density_shift))]
    x
}


# 4. Estimate adjusted SES and ten-year risks
profile_fit <- function() {
    original <- direction_data()
    keep <- unique(c("eid", "md", "md_z", "dtdi", "dtdi_z", "cross_ttwa", "urban_route", "density_shift", "density_z", 
        "ses_z", cov1, "date_attend", "date_death", "date_cvd"))
    original <- survdat(original[, ..keep], mort = TRUE)
    old <- profile_cluster(original)$data
    current <- profile_cluster(profile_geography(copy(original)))
    x <- current$data
    counts <- x[, .(N = .N), by = profile]
    reference <- counts[which.max(N), profile]
    od <- defs(make_icd(x$eid))
    fit <- function(outcome) {
        ses <- outcome == "SES"
        d <- if (ses) 
            copy(x[is.finite(ses_z)])
        else ukb_surv(x, od, outcome)
        ref <- if (ses) 
            copy(old[is.finite(ses_z)])
        else ukb_surv(old, od, outcome)
        supported <- function(a) a[, .(n = .N, events = if (ses) 
            NA_integer_
        else sum(event)), by = profile][, `:=`(supported, n >= 750L & (if (ses) 
            TRUE
        else events >= 20L))]
        previous <- supported(ref)
        ref <- ref[profile %in% previous[supported == TRUE, profile]]
        n <- supported(d)
        lev <- c(reference, setdiff(sort(n[supported == TRUE, profile]), reference))
        stopifnot(reference %in% n[supported == TRUE, profile])
        d <- d[profile %in% lev]
        d[, `:=`(group, factor(profile, levels = lev))]
        for (v in c("sex", "center", "gen_batch", "smoke", "alcohol")) set(d, j = v, value = droplevels(d[[v]]))
        f <- reformulate(c("group", if (ses) cov0 else cov1), if (ses) 
            "ses_z"
        else "Surv(time,event)")
        warnings <- character()
        m <- withCallingHandlers(if (ses) 
            lm(f, d)
        else coxph(f, d, ties = "efron", nocenter = NULL, control = coxph.control(iter.max = 60)), warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
        stopifnot(!length(warnings), all(is.finite(coef(m))), if (ses) 
            TRUE
        else m$iter < 60L)
        sm <- coef(summary(m))
        sm <- sm[paste0("group", lev[-1]), , drop = FALSE]
        a <- rbind(data.table(profile = lev[-1], beta = sm[, 1], se = sm[, if (ses) 
            2
        else 3], contrast_p = sm[, grep("^Pr", colnames(sm))[1]]), data.table(profile = reference, beta = 0, se = NA_real_, 
            contrast_p = NA_real_))
        a <- merge(n, a, by = "profile", all.x = TRUE)
        a[, `:=`(outcome = outcome, analysis_n = nrow(d), analysis_events = if (ses) 
            NA_integer_
        else sum(d$event), relative_estimate = if (ses) 
            beta
        else exp(beta), relative_lower = if (ses) 
            beta - 1.96 * se
        else exp(beta - 1.96 * se), relative_upper = if (ses) 
            beta + 1.96 * se
        else exp(beta + 1.96 * se))]
        ref[, `:=`(group, factor(reference, levels = lev))]
        b <- setNames(c(0, coef(m)[paste0("group", lev[-1])]), lev)
        if (ses) 
            pred <- mean(predict(m, newdata = ref)) + b
        else {
            h <- basehaz(m, centered = TRUE)
            H10 <- tail(h$hazard[h$time <= 10], 1)
            lp <- predict(m, newdata = ref, type = "lp", reference = "sample")
            pred <- vapply(b, function(v) 100 * mean(-expm1(-H10 * exp(lp + v))), 0)
        }
        a[, `:=`(estimate = unname(pred[profile]), lower = NA_real_, upper = NA_real_, p = NA_real_, measure = if (ses) 
            "Adjusted SES mean (SD)"
        else "10-year net risk (%)", lowest = FALSE)]
        best <- a[supported == TRUE][if (ses) 
            which.max(estimate)
        else which.min(estimate), profile]
        a[profile == best, `:=`(lowest, TRUE)]
        a
    }
    results <- rbindlist(lapply(c("SES", forest_out), fit))
    results[, `:=`(q, p.adjust(contrast_p, "BH"))]
    description <- x[, .(N = .N, SES_n = sum(is.finite(ses_z)), MD = median(md), MD_low = quantile(md, 0.25), MD_high = quantile(md, 
        0.75), Delta = median(dtdi), Delta_low = quantile(dtdi, 0.25), Delta_high = quantile(dtdi, 0.75), Density = median(density_shift), 
        Density_low = quantile(density_shift, 0.25), Density_high = quantile(density_shift, 0.75), Cross = unique(cross_ttwa), 
        Route = as.character(unique(urban_route))), by = profile]
    setorder(description, profile)
    boxes <- rbindlist(lapply(c("md", "dtdi", "density_shift"), function(v) x[, {
        a <- get(v)
        q <- as.numeric(quantile(a, c(0.25, 0.5, 0.75)))
        iqr <- q[3] - q[1]
        list(N = .N, minimum = min(a), lower = min(a[a >= q[1] - 1.5 * iqr]), q1 = q[1], median = q[2], q3 = q[3], 
            upper = max(a[a <= q[3] + 1.5 * iqr]), maximum = max(a), outliers = sum(a < q[1] - 1.5 * iqr | a > 
                q[3] + 1.5 * iqr), variable = v)
    }, by = profile]))
    r <- list(description = description, boxes = boxes, results = results, best = results[lowest == TRUE], targets = c("SES", 
        forest_out), reference = reference, N = nrow(x), SES_n = sum(is.finite(x$ses_z)), clusters = current$clusters)
    saveRDS(r, file.path(out, "profiles.rds"), compress = "gzip")
    save_xlsx(list(Profiles = description, Distributions = boxes, Outcomes = results), file.path(out, "profiles.xlsx"))
    r
}

if (sys.nframe() == 0L) profile_fit()

