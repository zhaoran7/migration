# 1. Regional mixed-effects Cox models
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

use(c("assoc"))

suppressPackageStartupMessages({
    library(coxme)
    library(survival)
})

region_root <- function() root

region_dir <- function() Sys.getenv("MOVE_REGION", file.path(root, "data/region"))

region_path <- function(part, name) file.path(region_dir(), part, name)

region_require <- function(packages) {
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) 
        stop("Missing R packages: ", paste(missing, collapse = ", "))
    invisible(TRUE)
}

region_save <- function(x, file) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    tmp <- paste0(file, ".tmp")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(x, tmp, compress = "gzip")
    if (!file.rename(tmp, file)) 
        stop("Cannot finalize: ", file)
    invisible(file)
}


# 2. Prepare regional model inputs
region_data <- function(run = FALSE) {
    if (!isTRUE(run)) 
        stop("Preparation is explicit: region_data(run=TRUE).")
    region_require(c("data.table", "jsonlite"))
    phe <- file.path(region_root(), "data/phe")
    for (p in c("data", "fit", "tests", "res")) dir.create(file.path(region_dir(), p), recursive = TRUE, showWarnings = FALSE)
    x <- fig2_data()
    x <- x[complete.cases(x[, c(region_cov, "origin", "destination"), with = FALSE])]
    od <- defs(make_icd(x$eid))
    nm <- jsonlite::fromJSON(file.path(region_root(), "data/ref/ttwa.geojson"))$features$properties
    nm <- data.table::as.data.table(nm)[, .(code = TTWA11CD, area = TTWA11NM)]
    coverage <- data.table::rbindlist(lapply(c("origin", "destination"), function(v) x[, .(n = .N), by = .(code = get(v))][, 
        `:=`(context, v)]))
    coverage <- merge(coverage, nm, by = "code", all.x = TRUE)
    keepo <- coverage[context == "origin" & n >= 100, code]
    keepd <- coverage[context == "destination" & n >= 100, code]
    config <- list(endpoints = forest_out, covariates = region_cov, formula = paste("Surv(time,event) ~", paste(region_cov, 
        collapse = " + ")), model = "MD + Delta TDI; no Birthplace/Recruitment TDI; no Cross TTWA", censor_date = as.character(end), 
        min_n = 100L, min_events = 10L, source = tools::md5sum(ukb_inputs()), created = as.character(Sys.time()), 
        R = R.version.string)
    samples <- list()
    for (i in seq_along(forest_out)) {
        o <- forest_out[i]
        d <- ukb_surv(x, od, o)
        oldn <- nrow(d)
        olde <- sum(d$event)
        cap <- pmin(as.numeric(end - d$date_attend)/365.25, as.numeric(d$date_lost - d$date_attend)/365.25, na.rm = TRUE)
        if (!o %in% c("All-cause mortality", "CVD mortality")) 
            cap <- pmin(cap, as.numeric(d$date_death - d$date_attend)/365.25, na.rm = TRUE)
        if (o == "CVD mortality") {
            other <- !is.na(d$date_death) & (is.na(d$date_cvd) | d$date_cvd != d$date_death)
            cap[other] <- pmin(cap[other], as.numeric(d$date_death[other] - d$date_attend[other])/365.25)
        }
        d[, `:=`(cap, cap)]
        d[, `:=`(event, as.integer(event & time <= cap))]
        d[, `:=`(time, pmin(time, cap))]
        d <- d[is.finite(cap) & cap > 0 & is.finite(time) & time > 0]
        d[, `:=`(supported, origin %in% keepo & destination %in% keepd)]
        vars <- c("time", "event", "cap", "origin", "destination", "supported", config$covariates)
        d <- d[, ..vars]
        d[, `:=`(c("origin", "destination"), lapply(.SD, factor)), .SDcols = c("origin", "destination")]
        for (v in names(d)) if (is.factor(d[[v]])) 
            data.table::set(d, j = v, value = droplevels(d[[v]]))
        stopifnot(!anyNA(d), all(d$time <= d$cap), all(d$event %in% 0:1))
        qc <- data.table::rbindlist(lapply(c("origin", "destination"), function(v) d[, .(n = .N, events = sum(event), 
            person_years = sum(time)), by = .(code = as.character(get(v)))][, `:=`(context, v)]))
        qc <- merge(qc, nm, by = "code", all.x = TRUE)
        stopifnot(!anyNA(qc$area))
        qc[, `:=`(reliable, n >= config$min_n & events >= config$min_events)]
        item <- list(id = i, outcome = o, data = d, coverage = qc, config = config)
        region_save(item, region_path("data", paste0(i, ".rds")))
        samples[[i]] <- data.table::data.table(id = i, outcome = o, n = nrow(d), events = sum(d$event), original_n = oldn, 
            original_events = olde, supported_n = sum(d$supported), followup_median = median(d$time), followup_min = min(d$time), 
            followup_max = max(d$time))
        message(o, ": ", nrow(d), " participants; ", sum(d$event), " events")
    }
    region_save(list(config = config, coverage = coverage, samples = data.table::rbindlist(samples)), region_path("data", 
        "design.rds"))
    invisible(data.table::rbindlist(samples))
}


# 3. Fit crossed regional random effects
region_fitmodel <- function(d, f, initial = NULL, start = NULL) {
    region_require(c("survival", "coxme"))
    Surv <- survival::Surv
    environment(f) <- environment()
    ff <- update(f, . ~ . + (1 | origin) + (1 | destination))
    warnings <- character()
    args <- list(formula = ff, data = d, ties = "efron", control = coxme::coxme.control(iter.max = 40))
    if (!is.null(initial)) 
        args$vinit <- as.list(unname(pmax(initial, 1e-05)))
    if (!is.null(start)) 
        args$init <- unname(start)
    m <- withCallingHandlers(do.call(coxme::coxme, args), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    if (length(m$fail)) 
        stop(paste(m$fail, collapse = "; "))
    if (any(!is.finite(c(m$coefficients, m$loglik, unlist(m$vcoef))))) 
        stop("Nonfinite Cox fit")
    if (any(grepl("converg|iteration|infinite|singular", warnings, ignore.case = TRUE))) 
        stop(paste(warnings, collapse = "; "))
    attr(m, "warnings") <- warnings
    m
}

region_effect <- function(m, coverage) {
    region_require("data.table")
    v <- diag(as.matrix(m$var))
    k <- 0L
    z <- data.table::rbindlist(lapply(names(m$frail), function(context) {
        b <- m$frail[[context]]
        ix <- k + seq_along(b)
        k <<- k + length(b)
        data.table::data.table(context, code = names(b), beta = as.numeric(b), se = sqrt(v[ix]))
    }))
    z <- merge(z, coverage, by = c("context", "code"), sort = FALSE)
    z[, `:=`(HR = exp(beta), conditional_lower = exp(beta - 1.96 * se), conditional_upper = exp(beta + 1.96 * se), 
        p_approx = 2 * pnorm(-abs(beta/se)))]
    z
}

region_compact_fit <- function(x) {
    fields <- c("id", "outcome", "version", "module_version", "input_hash", "n", "events", "area", "variance", 
        "config", "packages", "seconds", "lrt")
    z <- x[intersect(fields, names(x))]
    z$model <- x$model[c("coefficients", "loglik", "vcoef")]
    z
}

region_fit <- function(id, run = FALSE, check_version = FALSE, fresh = FALSE) {
    region_require(c("data.table", "survival", "coxme"))
    stopifnot(length(id) == 1L, id %in% 1:10)
    input <- region_path("data", paste0(id, ".rds"))
    hash <- unname(tools::md5sum(input))
    if (is.na(hash)) 
        stop("Missing endpoint data: ", input)
    file <- region_path("fit", paste0(id, ".rds"))
    if (!fresh && file.exists(file)) {
        previous <- readRDS(file)
        valid <- identical(previous$input_hash, hash) && (!isTRUE(check_version) || identical(previous$version, 
            region_version))
        stopifnot(all(c("id", "outcome", "n", "events", "area", "variance", "model") %in% names(previous)), identical(names(previous$variance), 
            c("origin", "destination")), all(is.finite(c(previous$model$loglik, previous$model$coefficients, previous$variance))))
        if (valid) 
            return(invisible(region_compact_fit(previous)))
    }
    if (!isTRUE(run)) 
        stop("Missing fit; explicitly use region_fit(id, run=TRUE).")
    a <- readRDS(input)
    d <- a$data
    Surv <- survival::Surv
    f <- as.formula(a$config$formula, env = environment())
    message(a$outcome, ": full crossed Cox fit (all areas)")
    timer <- proc.time()[3]
    m <- region_fitmodel(d, f)
    area <- region_effect(m, a$coverage)
    variance <- unlist(m$vcoef, use.names = FALSE)
    stopifnot(identical(names(m$vcoef), c("origin", "destination")))
    names(variance) <- c("origin", "destination")
    result <- region_compact_fit(list(id = id, outcome = a$outcome, version = region_version, module_version = region_module_version, 
        input_hash = hash, n = nrow(d), events = sum(d$event), model = m, area = area, variance = variance, seconds = proc.time()[3] - 
            timer, config = a$config, packages = c(coxme = as.character(utils::packageVersion("coxme")), survival = as.character(utils::packageVersion("survival")))))
    region_save(result, file)
    invisible(result)
}


# 4. Test regional variance
region_test_one <- function(a, full, context) {
    region_require(c("data.table", "survival", "coxme"))
    stopifnot(context %in% c("origin", "destination"))
    d <- a$data
    stopifnot(nrow(d) == full$n, sum(d$event) == full$events, !anyNA(d))
    Surv <- survival::Surv
    f <- as.formula(a$config$formula, env = environment())
    keep <- setdiff(c("origin", "destination"), context)
    f <- update(f, paste(". ~ . + (1 |", keep, ")"))
    ll <- full$model$loglik[2]
    mh <- exp(qnorm(0.75) * sqrt(2 * full$variance[[context]]))
    timer <- proc.time()[3]
    warnings <- character()
    m <- withCallingHandlers(coxme::coxme(f, data = d, ties = "efron", vinit = list(max(full$variance[[keep]], 
        1e-05)), control = coxme::coxme.control(iter.max = 40)), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    stopifnot(!length(m$fail), all(is.finite(c(m$loglik, m$coefficients, unlist(m$vcoef)))), isTRUE(all.equal(as.numeric(m$y[, 
        1]), d$time)), isTRUE(all.equal(as.numeric(m$y[, 2]), as.numeric(d$event))), identical(names(m$coefficients), 
        names(full$model$coefficients)))
    if (any(grepl("converg|iteration|infinite|singular", warnings, ignore.case = TRUE))) 
        stop(paste(warnings, collapse = "; "))
    lr <- 2 * (ll - m$loglik[2])
    if (!is.finite(lr) || lr <= -0.05) 
        stop("Materially negative LR; inspect optimizer before publishing.")
    list(result = data.table::data.table(id = a$id, outcome = a$outcome, context, n = full$n, events = full$events, 
        MHR = mh, variance = full$variance[[context]], full_ll = ll, reduced_ll = m$loglik[2], LR_raw = lr, LR = max(lr, 
            0), df = 1, p = pchisq(max(lr, 0), 1, lower.tail = FALSE), seconds = proc.time()[3] - timer, warnings = paste(warnings, 
            collapse = "; ")), kept = keep, reduced_variance = m$vcoef, reduced_fixed = m$coefficients, input_hash = full$input_hash, 
        formula = deparse(f), packages = c(coxme = as.character(utils::packageVersion("coxme"))))
}

region_method <- function() list(test = "Conditional integrated-likelihood ratio test; conservative chi-square reference, 1 df", 
    correction = "Benjamini-Hochberg across all 20 outcome-by-context tests", selection = "Plot each outcome-by-context pair with q < 0.05", 
    boundary = "Small negative LR from numerical approximation truncated to zero; LR_raw must exceed -0.05", interval = "Approximate pointwise conditional EB interval; variance-component uncertainty not included", 
    area_test = "Exploratory normal-tail test of conditional random effect, BH across all estimated areas and outcomes", 
    mhr_interval = "Not estimated", maps = "All estimated TTWAs, including unreliable ones; no local significance mask; unestimated TTWAs grey", 
    labels = "Extreme EB effects with spatial separation, not significance-selected city tests", joint = "Legacy two-df joint LRT retained for audit only; NOT used for panels or their q values")

region_read <- function(cache = region_path("res", "tests.rds"), verify_hash = FALSE) {
    region_require("data.table")
    if (!file.exists(cache)) 
        stop("Missing numerical cache: ", cache, "; no automatic fitting.")
    d <- readRDS(cache)
    design <- region_path("data", "design.rds")
    if (file.exists(design)) 
        stopifnot(setequal(readRDS(design)$config$covariates, region_cov))
    stopifnot(all(c("tests", "models", "area", "joint") %in% names(d)))
    for (nm in c("tests", "area", "joint")) d[[nm]] <- data.table::as.data.table(d[[nm]])
    t <- d$tests
    a <- d$area
    r <- d$joint
    stopifnot(nrow(t) == 20L, nrow(r) == 10L, setequal(r$id, 1:10), !anyDuplicated(t[, .(id, context)]), !anyDuplicated(r$id), 
        setequal(t$context, c("origin", "destination")), all(t[, .N, by = id]$N == 2L), setequal(t$id, r$id), all(is.finite(t$q)), 
        all(t$warnings == ""), all(t$LR_raw > -0.05), all(t$df == 1), all(t$variance >= 0), isTRUE(all.equal(t$LR, 
            pmax(t$LR_raw, 0))), isTRUE(all.equal(t$LR_raw, 2 * (t$full_ll - t$reduced_ll))), isTRUE(all.equal(t$p, 
            pchisq(t$LR, 1, lower.tail = FALSE))), isTRUE(all.equal(t$q, p.adjust(t$p, "BH"))), isTRUE(all.equal(t$MHR, 
            exp(qnorm(0.75) * sqrt(2 * t$variance)))), !anyDuplicated(a[, .(id, context, code)]), isTRUE(all.equal(a$q_approx, 
            p.adjust(a$p_approx, "BH"))), isTRUE(all.equal(a$HR, exp(a$beta))), length(d$models) == 20L)
    for (j in seq_len(nrow(t))) {
        row <- t[j]
        s <- r[id == row$id]
        stopifnot(identical(row$outcome, s$outcome), row$n == s$n, row$events == s$events, isTRUE(all.equal(row$MHR, 
            s[[row$context]])))
        found <- Filter(function(m) m$result$id == row$id && m$result$context == row$context, d$models)
        stopifnot(length(found) == 1L, isTRUE(all.equal(found[[1]]$result$LR_raw, row$LR_raw)), isTRUE(all.equal(found[[1]]$result$MHR, 
            row$MHR)))
        if (isTRUE(verify_hash)) {
            hash <- unname(tools::md5sum(region_path("data", paste0(row$id, ".rds"))))
            stopifnot(!is.na(hash), identical(found[[1]]$input_hash, hash))
        }
    }
    d$cache_provenance <- list(path = normalizePath(cache), md5 = unname(tools::md5sum(cache)), input_hashes_verified = isTRUE(verify_hash), 
        module_version = region_module_version)
    d
}


# 5. Collect estimates and corrected tests
region_tests <- function(run = FALSE, cache = region_path("res", "tests.rds"), fresh = FALSE) {
    if (!isTRUE(run)) 
        return(region_read(cache))
    region_require(c("data.table", "survival", "coxme"))
    previous <- if (file.exists(cache)) 
        region_read(cache)
    else NULL
    models <- areas <- joint <- vector("list", 10L)
    for (id in 1:10) {
        a <- readRDS(region_path("data", paste0(id, ".rds")))
        full <- region_fit(id, run = TRUE, fresh = fresh)
        per_id <- lapply(c("destination", "origin"), function(context) {
            file <- region_path("tests", paste0(id, "_", context, ".rds"))
            old <- if (fresh) 
                list()
            else if (file.exists(file)) 
                list(readRDS(file))
            else Filter(function(z) z$result$id == id && z$result$context == context, previous$models)
            if (length(old) == 1L && identical(old[[1]]$input_hash, full$input_hash) && isTRUE(all.equal(old[[1]]$result$full_ll, 
                full$model$loglik[2])) && isTRUE(all.equal(old[[1]]$result$variance, full$variance[[context]])) && 
                old[[1]]$result$n == full$n && old[[1]]$result$events == full$events) 
                return(old[[1]])
            message(a$outcome, ": conditional test of ", context)
            z <- region_test_one(a, full, context)
            region_save(z, file)
            z
        })
        models[[id]] <- per_id
        ar <- data.table::copy(full$area)
        ar[, `:=`(id = a$id, outcome = a$outcome)]
        areas[[id]] <- ar
        joint[[id]] <- data.table::data.table(id = a$id, outcome = a$outcome, n = full$n, events = full$events, 
            origin = exp(qnorm(0.75) * sqrt(2 * full$variance[["origin"]])), destination = exp(qnorm(0.75) * sqrt(2 * 
                full$variance[["destination"]])), LR = if (is.null(full$lrt)) 
                NA_real_
            else full$lrt)
        rm(full, a)
        invisible(gc(FALSE))
    }
    models <- unlist(models, recursive = FALSE)
    t <- data.table::rbindlist(lapply(models, `[[`, "result"))
    t[, `:=`(q, p.adjust(p, "BH"))]
    a <- data.table::rbindlist(areas)
    a[, `:=`(q_approx, p.adjust(p_approx, "BH"))]
    r <- data.table::rbindlist(joint)
    r[, `:=`(p = NA_real_, q = NA_real_)]
    same <- c("id", "outcome", "n", "events", "origin", "destination")
    if (!is.null(previous) && isTRUE(all.equal(r[, ..same], previous$joint[, ..same]))) 
        r <- data.table::copy(previous$joint)
    result <- list(tests = t, models = models, area = a, joint = r, method = region_method(), module_version = region_module_version)
    region_save(result, cache)
    region_read(cache)
}


# 6. Link aggregate effects to boundaries
region_geometry <- function(file = file.path(region_root(), "data/ref/ttwa.geojson")) {
    region_require(c("data.table", "jsonlite"))
    areas <- jsonlite::fromJSON(file, simplifyVector = FALSE)
    stopifnot(identical(areas$crs$properties$name, "EPSG:27700"))
    poly <- data.table::rbindlist(lapply(areas$features, function(f) {
        g <- f$geometry
        stopifnot(g$type %in% c("Polygon", "MultiPolygon"))
        p <- if (g$type == "Polygon") 
            list(g$coordinates)
        else g$coordinates
        data.table::rbindlist(lapply(seq_along(p), function(j) data.table::rbindlist(lapply(seq_along(p[[j]]), 
            function(k) {
                z <- do.call(rbind, lapply(p[[j]][[k]], unlist))
                data.table::data.table(code = f$properties$TTWA11CD, x = z[, 1], y = z[, 2], piece = paste(f$properties$TTWA11CD, 
                  j, sep = "_"), ring = paste(f$properties$TTWA11CD, j, k, sep = "_"))
            }))))
    }))
    cen <- data.table::rbindlist(lapply(areas$features, function(f) {
        ps <- if (f$geometry$type == "Polygon") 
            list(f$geometry$coordinates)
        else f$geometry$coordinates
        z <- data.table::rbindlist(lapply(ps, function(p) {
            xy <- do.call(rbind, lapply(p[[1]], unlist))
            storage.mode(xy) <- "double"
            if (any(xy[1, ] != xy[nrow(xy), ])) 
                xy <- rbind(xy, xy[1, ])
            x <- xy[, 1]
            y <- xy[, 2]
            k <- length(x)
            cross <- x[-k] * y[-1] - x[-1] * y[-k]
            ar <- sum(cross)/2
            data.table::data.table(area = abs(ar), x = sum((x[-k] + x[-1]) * cross)/(6 * ar), y = sum((y[-k] + 
                y[-1]) * cross)/(6 * ar))
        }))
        data.table::data.table(code = f$properties$TTWA11CD, name = f$properties$TTWA11NM, x = weighted.mean(z$x, 
            z$area), y = weighted.mean(z$y, z$area))
    }))
    stopifnot(!anyDuplicated(cen$code), all(is.finite(cen$x)), all(is.finite(cen$y)))
    list(poly = poly, centroids = cen, crs = areas$crs$properties$name, file = normalizePath(file), md5 = unname(tools::md5sum(file)))
}

run_region <- function(force = FALSE, cache = region_path("res", "tests.rds"), verify_hash = FALSE) {
    if (!force && !stale(cache, ukb_inputs())) 
        return(region_read(cache, verify_hash = verify_hash))
    if (!all(file.exists(region_path("data", paste0(1:10, ".rds")))) || stale(region_path("data", "design.rds"), 
        ukb_inputs())) 
        region_data(run = TRUE)
    region_tests(run = TRUE, cache = cache, fresh = force)
}

if (sys.nframe() == 0L) {
    a <- commandArgs(TRUE)
    stage <- if (length(a)) 
        a[1]
    else "all"
    if (stage == "data") 
        region_data(run = TRUE)
    else if (stage == "fit") 
        region_fit(as.integer(a[2]), run = TRUE, fresh = TRUE)
    else if (stage == "tests") 
        region_tests(run = TRUE, fresh = FALSE)
    else if (stage == "all") {
        region_data(run = TRUE)
        run_region(force = TRUE)
    }
    else stop("Use data, fit <1-10>, or tests")
}

