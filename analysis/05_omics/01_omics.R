# 1. Protein and MRI associations and mediation
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

use("assoc")

# 2. Define imaging timing and covariates
med_visit <- function() {
    file <- file.path(phe, "mri_visit.rds")
    if (!file.exists(file)) {
        raw <- Sys.getenv("MOVE_UKB_RAW", file.path(ukb, "pheno.tsv"))
        x <- fread(raw, select = c("eid", "p53_i2", "p54_i2", "p21003_i2", "p25000_i2"))
        setnames(x, c("eid", "image_date", "image_center", "image_age", "head_scale"))
        x[, `:=`(image_date, as.IDate(image_date))]
        saveRDS(x, file, compress = "gzip")
    }
    x <- as.data.table(readRDS(file))
    x[, `:=`(eid = as.character(eid), image_date = as.IDate(image_date))]
    x
}

med_rhs <- function(layer, head = FALSE) {
    scan <- if (layer == "MRI") 
        c("image_age_z", "image_center", if (head) "head_scale_z")
    list(md = c("md_z", "dtdi_z", fig2_cov1, scan), dtdi = c("md_z", "dtdi_z", tdi_cov1, scan), ttwa = c("cross_ttwa", 
        "dtdi_z", cov1, scan))
}

med_design <- function(x, rhs) {
    X <- model.matrix(reformulate(rhs), droplevels(as.data.frame(x[, ..rhs])))
    qx <- qr(X)
    if (qx$rank != ncol(X)) 
        stop("Omics model is rank deficient")
    V <- chol2inv(qr.R(qx))
    dimnames(V) <- list(colnames(X), colnames(X))
    list(X = X, qr = qx, V = V, df = nrow(X) - ncol(X))
}


# 3. Estimate feature associations
med_lm <- function(om, x, rhs, info, layer, chunk = 160L) {
    f <- intersect(setdiff(names(om), "eid"), info[["Feature ID"]])
    head <- if (layer == "MRI") 
        info$Measure[match(f, info[["Feature ID"]])] %in% c("area", "volume")
    else rep(FALSE, length(f))
    groups <- split(f, head)
    ans <- rbindlist(lapply(names(groups), function(group) {
        design <- lapply(med_rhs(layer, group == "TRUE"), function(v) med_design(x, v))
        rbindlist(lapply(split(groups[[group]], ceiling(seq_along(groups[[group]])/chunk)), function(g) {
            Y <- apply(as.matrix(om[, ..g]), 2, z)
            if (is.null(dim(Y))) 
                Y <- matrix(Y, ncol = 1, dimnames = list(NULL, g))
            rbindlist(lapply(names(design), function(id) {
                d <- design[[id]]
                B <- qr.coef(d$qr, Y)
                if (is.null(dim(B))) 
                  B <- matrix(B, ncol = 1, dimnames = list(colnames(d$X), g))
                E <- Y - d$X %*% B
                vars <- switch(id, md = "md_z", dtdi = "dtdi_z", ttwa = "cross_ttwa")
                rbindlist(lapply(vars, function(v) {
                  se <- sqrt(d$V[v, v] * colSums(E^2)/d$df)
                  data.table(feature = g, layer, var = v, exposure = unname(med_exp[v]), n = nrow(Y), beta = as.numeric(B[v, 
                    ]), se = as.numeric(se), p = 2 * stats::pt(-abs(B[v, ]/se), d$df))
                }))
            }))
        }))
    }))
    ans[, `:=`(p_bonf, pmin(1, p * length(f) * length(med_exp)))]
    merge(ans, info, by.x = "feature", by.y = "Feature ID", all.x = TRUE, sort = FALSE)
}

med_data <- function(layer, followup = TRUE) {
    x <- fig2_data()
    x[, `:=`(eid, as.character(eid))]
    if (layer == "Protein") {
        om <- as.data.table(readRDS(file.path(phe, "prot.rds")))
        info <- as.data.table(readRDS(file.path(phe, "omics_info.rds"))$Protein)
    }
    else {
        x <- merge(x, med_visit()[, .(eid, image_date, image_center, image_age, head_scale)], by = "eid", all = FALSE, 
            sort = FALSE)
        if (followup) 
            x <- x[image_date < end]
        x[, `:=`(image_age_z = z(image_age), image_center = factor(image_center), head_scale_z = z(head_scale))]
        om <- as.data.table(readRDS(file.path(phe, "mri.rds")))
        info <- as.data.table(readRDS(file.path(phe, "omics_info.rds"))$MRI)
    }
    info <- unique(info, by = "Feature ID")
    om[, `:=`(eid, as.character(eid))]
    i <- match(x$eid, om$eid)
    keep <- !is.na(i)
    x <- x[keep]
    om <- om[i[keep]]
    f <- intersect(setdiff(names(om), "eid"), info[["Feature ID"]])
    rhs <- med_rhs(layer)
    required <- unique(unlist(med_rhs(layer, head = layer == "MRI")))
    keep <- complete.cases(x[, ..required]) & complete.cases(om[, ..f])
    x <- x[keep]
    om <- om[keep]
    factors <- intersect(c("sex", "center", "gen_batch", "smoke", "alcohol", "image_center"), names(x))
    x[, `:=`((factors), lapply(.SD, droplevels)), .SDcols = factors]
    x[, `:=`(rowid, .I)]
    list(ph = x, om = om, info = info, rhs = rhs, features = f)
}


# 4. Define post-measurement follow-up
med_surv <- function(x, ed = NULL, mort = FALSE, cvd = FALSE) {
    d <- copy(x)
    if (mort) {
        d[, `:=`(censor, fifelse(!is.na(date_death) & date_death <= end, date_death, end))]
        d[, `:=`(event, if (cvd) 
            !is.na(date_cvd) & date_cvd > image_date & date_cvd <= censor
        else !is.na(date_death) & date_death > image_date & date_death <= end)]
        d[, `:=`(prev = FALSE, end0 = fifelse(event, if (cvd) 
            date_cvd
        else date_death, censor))]
    }
    else {
        d <- merge(d, ed, by = "eid", all.x = TRUE, sort = FALSE)
        d[, `:=`(censor, fifelse(!is.na(date_death) & date_death <= end, date_death, end))]
        d[, `:=`(prev = !is.na(event_date) & event_date <= image_date, event = !is.na(event_date) & event_date > 
            image_date & event_date <= censor, end0 = fifelse(!is.na(event_date) & event_date > image_date & event_date <= 
            censor, event_date, censor))]
    }
    d[, `:=`(time, as.numeric(as.IDate(end0) - image_date)/365.25)]
    d[!prev & is.finite(time) & time > 0]
}

med_cox <- function(X, y) {
    focus <- intersect(c("md_z", "dtdi_z", "cross_ttwa", "mediator0"), colnames(X))
    one <- function(eps, init = NULL) {
        warn <- character()
        fit <- withCallingHandlers(survival::coxph.fit(X, y, strata = NULL, offset = rep(0, nrow(X)), init = init, 
            control = coxph.control(iter.max = 100, eps = eps), weights = rep(1, nrow(X)), method = "efron", rownames = NULL, 
            resid = TRUE), warning = function(w) {
            warn <<- c(warn, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
        list(fit = fit, warnings = unique(warn))
    }
    a <- one(1e-09)
    fit <- a$fit
    warn <- a$warnings
    target <- function(f) c(f$coefficients[focus], sqrt(diag(f$var))[match(focus, colnames(X))])
    valid <- all(is.finite(target(fit))) && all(diag(fit$var)[match(focus, colnames(X))] > 0)
    change <- 0
    boundary <- character()
    if (length(warn) || fit$iter >= 100L) {
        b <- one(1e-11, if (all(is.finite(fit$coefficients))) 
            fit$coefficients
        else NULL)
        change <- max(abs(target(b$fit) - target(fit))/pmax(1, abs(target(fit))))
        valid <- valid && all(is.finite(target(b$fit))) && is.finite(change) && change < 1e-05 && abs(diff(c(tail(fit$loglik, 
            1), tail(b$fit$loglik, 1)))) < 1e-05 && b$fit$iter < 100L
        warn <- unique(c(warn, b$warnings))
        fit <- b$fit
        for (w in warn) {
            if (!grepl("coefficient may be infinite", w, fixed = TRUE)) {
                valid <- FALSE
                next
            }
            idx <- as.integer(regmatches(strsplit(w, ";", fixed = TRUE)[[1]][1], gregexpr("[0-9]+", strsplit(w, 
                ";", fixed = TRUE)[[1]][1]))[[1]])
            boundary <- union(boundary, colnames(X)[idx])
        }
        if (any(boundary %in% focus)) 
            valid <- FALSE
    }
    attr(fit, "qc") <- list(valid = valid, status = if (!valid) "Unstable fit" else if (length(warn)) "Nuisance boundary; target stable" else "OK", 
        warnings = paste(warn, collapse = "; "), boundary = paste(boundary, collapse = ", "), iterations = fit$iter, 
        change = change)
    fit
}

med_influence <- function(fit, X, y, vars) {
    fit$x <- X
    fit$y <- y
    fit$method <- "efron"
    fit$terms <- terms(~1)
    class(fit) <- "coxph"
    score <- residuals(fit, type = "score")
    if (is.null(dim(score))) 
        score <- matrix(score, ncol = ncol(X))
    ans <- score %*% fit$var[, match(vars, colnames(X)), drop = FALSE]
    colnames(ans) <- vars
    ans
}

med_interval <- function(a, b, c, ia, ib, ic, seed) {
    S <- crossprod(cbind(ia, ib, ic))
    S <- (S + t(S))/2
    eig <- eigen(S, symmetric = TRUE)
    set.seed(seed)
    Z <- matrix(rnorm(30000), 10000, 3)
    sim <- sweep(Z %*% diag(sqrt(pmax(eig$values, 0)), 3) %*% t(eig$vectors), 2, c(a, b, c), "+")
    ind <- sim[, 1] * sim[, 2]
    ratio <- ind/sim[, 3]
    iq <- quantile(ind, c(0.025, 0.975))
    pq <- quantile(ratio, c(0.025, 0.975))
    list(indirect_lower = unname(iq[1]), indirect_upper = unname(iq[2]), proportion_lower = unname(pq[1]), proportion_upper = unname(pq[2]), 
        denominator_crossing = mean(sign(sim[, 3]) != sign(c)), covariance_ab = S[1, 2], covariance_ac = S[1, 3], 
        covariance_bc = S[2, 3])
}


# 5. Estimate individual-feature mediation
med_fit <- function(layer, d, association_data = d, variables = names(med_exp)) {
    assoc <- med_lm(association_data$om, association_data$ph, association_data$rhs, association_data$info, layer)
    cand <- assoc[p_bonf < 0.05 & var %in% variables]
    f <- unique(cand$feature)
    if (!length(f)) 
        return(list(association = assoc, mediation = data.table()))
    Y <- apply(as.matrix(d$om[, ..f]), 2, z)
    if (is.null(dim(Y))) 
        Y <- matrix(Y, ncol = 1, dimnames = list(NULL, f))
    ph <- d$ph
    od <- defs(make_icd(ph$eid))
    one <- function(outcome) {
        message(layer, " ", outcome)
        s <- if (layer == "Protein") 
            ukb_surv(ph, od, outcome)
        else if (outcome == "All-cause mortality") 
            med_surv(ph, mort = TRUE)
        else if (outcome == "CVD mortality") 
            med_surv(ph, mort = TRUE, cvd = TRUE)
        else med_surv(ph, od[[outcome]])
        y <- Surv(s$time, s$event)
        yy <- Y[s$rowid, , drop = FALSE]
        keys <- c(md_z = "md", dtdi_z = "dtdi", cross_ttwa = "ttwa")
        ids <- unique(unname(keys[cand$var]))
        des <- lapply(if (layer == "MRI") 
            c(FALSE, TRUE)
        else FALSE, function(head) lapply(med_rhs(layer, head)[ids], function(rhs) {
            lm <- med_design(s, rhs)
            X <- lm$X[, -1, drop = FALSE]
            fit <- med_cox(X, y)
            B <- qr.coef(lm$qr, yy)
            if (is.null(dim(B))) 
                B <- matrix(B, ncol = 1, dimnames = list(colnames(lm$X), f))
            E <- yy - lm$X %*% B
            SE <- sqrt(outer(diag(lm$V), colSums(E^2)/lm$df))
            dimnames(SE) <- list(colnames(lm$X), colnames(yy))
            targets <- intersect(c("md_z", "dtdi_z", "cross_ttwa"), colnames(X))
            list(X = X, fit = fit, B = B, se = SE, lmX = lm$X, V = lm$V, E = E, df = lm$df, influence = if (isTRUE(attr(fit, 
                "qc")$valid)) med_influence(fit, X, y, targets) else NULL)
        }))
        names(des) <- as.character(if (layer == "MRI") c(FALSE, TRUE) else FALSE)
        rbindlist(lapply(f, function(g) {
            z0 <- cand[feature == g]
            keys <- c(md_z = "md", dtdi_z = "dtdi", cross_ttwa = "ttwa")
            rbindlist(lapply(unique(unname(keys[z0$var])), function(id) {
                head <- layer == "MRI" && z0$Measure[1] %in% c("area", "volume")
                a <- des[[as.character(head)]][[id]]
                vars <- intersect(z0$var, names(keys)[keys == id])
                if (!length(vars)) 
                  return(NULL)
                fit <- med_cox(cbind(a$X, mediator0 = yy[, g]), y)
                b <- unname(fit$coefficients["mediator0"])
                sb <- sqrt(fit$var[ncol(a$X) + 1, ncol(a$X) + 1])
                h <- copy(z0[var %in% vars])
                h[, `:=`(outcome = outcome, n = nrow(s), event = sum(s$event), a = as.numeric(a$B[var, g]), se_a = as.numeric(a$se[var, 
                  g]), b = b, se_b = sb, p_b = 2 * pnorm(-abs(b/sb)), total = as.numeric(a$fit$coefficients[var]), 
                  se_total = sqrt(diag(a$fit$var)[match(var, colnames(a$X))]), direct = as.numeric(fit$coefficients[var]))]
                qa <- attr(a$fit, "qc")
                qb <- attr(fit, "qc")
                df <- a$df
                h[, `:=`(model_valid = isTRUE(qa$valid) && isTRUE(qb$valid), total_status = qa$status, mediator_status = qb$status, 
                  total_warning = qa$warnings, mediator_warning = qb$warnings, boundary_terms = paste(unique(c(qa$boundary, 
                    qb$boundary)), collapse = "; "), target_change = max(qa$change, qb$change), p_a = 2 * stats::pt(-abs(a/se_a), 
                    df))]
                mult <- length(d$features) * length(med_exp)
                h[, `:=`(p_a_bonf = pmin(1, p_a * mult), p_b_bonf = pmin(1, p_b * mult * length(forest_out)), indirect_lower = NA_real_, 
                  indirect_upper = NA_real_, proportion_lower = NA_real_, proportion_upper = NA_real_, denominator_crossing = NA_real_, 
                  covariance_ab = NA_real_, covariance_ac = NA_real_, covariance_bc = NA_real_)]
                eligible <- which(h$model_valid & h$p_bonf < 0.05 & h$p_a_bonf < 0.05 & h$p_b_bonf < 0.05 & 2 * 
                  pnorm(-abs(h$a * h$b/sqrt((h$b * h$se_a)^2 + (h$a * h$se_b)^2))) * mult * length(forest_out) < 
                  0.05)
                if (length(eligible)) {
                  ib <- med_influence(fit, cbind(a$X, mediator0 = yy[, g]), y, "mediator0")[, 1]
                  for (j in eligible) {
                    v <- h$var[j]
                    ia <- drop(a$lmX %*% a$V[, v]) * a$E[, g]
                    interval <- med_interval(h$a[j], h$b[j], h$total[j], ia, ib, a$influence[, v], 20261019L + 
                      match(layer, c("Protein", "MRI")) * 100000L + match(outcome, forest_out) * 10000L + match(g, 
                      d$features) * 3L + match(v, names(med_exp)))
                    for (nm in names(interval)) set(h, i = j, j = nm, value = interval[[nm]])
                  }
                }
                h
            }), fill = TRUE)
        }), fill = TRUE)
    }
    ans <- rbindlist(mclapply(forest_out, one, mc.cores = max(1, min(core, 6))), fill = TRUE)
    stopifnot(nrow(ans) == nrow(cand) * length(forest_out))
    ans[, `:=`(indirect = a * b, se_indirect = sqrt((b * se_a)^2 + (a * se_b)^2))]
    ans[, `:=`(p_indirect = 2 * pnorm(-abs(indirect/se_indirect)), p_total = 2 * pnorm(-abs(total/se_total)), indirect_hr = exp(indirect), 
        total_hr = exp(total), direct_hr = exp(direct), proportion = indirect/total)]
    ans[, `:=`(p_indirect_bonf = pmin(1, p_indirect * length(d$features) * length(med_exp) * length(forest_out)), 
        p_total_bonf = pmin(1, p_total * length(med_exp) * length(forest_out)))]
    for (v in c("p_b", "p_b_bonf", "p_indirect", "p_indirect_bonf", "p_total", "p_total_bonf")) set(ans, i = which(!ans$model_valid), 
        j = v, value = NA_real_)
    ans <- med_select(ans)
    list(association = assoc, mediation = ans)
}


# 6. Apply numerical and statistical criteria
med_select <- function(x) {
    x <- copy(x)
    x[, `:=`(path_significant, model_valid & p_bonf < 0.05 & p_a_bonf < 0.05 & p_b_bonf < 0.05 & p_indirect_bonf < 
        0.05)]
    x[, `:=`(proportion_supported, is.finite(proportion) & proportion > 0 & proportion < 1)]
    x[, `:=`(proportion_precise, is.finite(proportion_lower) & proportion_lower > 0 & is.finite(proportion_upper) & 
        proportion_upper < 1)]
    x[, `:=`(selected, path_significant & proportion_supported & proportion > 0.05 & sign(indirect) == sign(total))]
    x[, `:=`(decision, fcase(!model_valid, "Unstable Cox target", p_bonf >= 0.05, "Exposure association not significant", 
        p_a_bonf >= 0.05, "Path a not significant", p_b_bonf >= 0.05, "Path b not significant", p_indirect_bonf >= 
            0.05, "Indirect association not significant", !is.finite(proportion) | proportion <= 0 | proportion >= 
            1, "Proportion outside 0-100%", proportion <= 0.05, "Proportion no greater than 5%", default = "Qualifying pathway"))]
    x
}

med_pick <- function(x) if ("selected" %in% names(x)) x[selected == TRUE] else x[0]

med_show <- function(x, max_nodes = 40L) {
    x <- med_pick(copy(x))[exposure %in% med_view]
    if (!nrow(x)) 
        return(x)
    rank <- x[, .(prop = max(proportion), p = min(p_indirect_bonf)), by = feature]
    setorder(rank, -prop, p, feature)
    md <- rank[feature %in% x[exposure == "MD", feature], feature]
    ids <- head(c(md, setdiff(rank$feature, md)), max_nodes)
    x[feature %in% ids]
}


# 7. Save the primary results
run_fig6_fit <- function(fresh = FALSE, layers = c("Protein", "MRI")) {
    file <- file.path(out, "Fig6.rds")
    r <- if (!fresh && file.exists(file)) 
        readRDS(file)
    else list()
    if (!identical(r$model, med_model)) 
        r <- list(model = med_model)
    for (layer in layers) if (fresh || is.null(r[[layer]])) {
        message("Omics: ", layer)
        d <- med_data(layer)
        r[[layer]] <- if (layer == "MRI") 
            med_fit(layer, d, med_data(layer, FALSE))
        else med_fit(layer, d)
        saveRDS(r, file, compress = "gzip")
    }
    for (layer in intersect(c("Protein", "MRI"), names(r))) if (nrow(r[[layer]]$mediation)) 
        r[[layer]]$mediation <- med_select(r[[layer]]$mediation)
    for (layer in intersect(c("Protein", "MRI"), names(r))) r[[layer]]$association <- unique(r[[layer]]$association, 
        by = c("exposure", "feature"))
    saveRDS(r, file, compress = "gzip")
    r
}

supp_pairs <- function() data.table(comparison = c("Low to high vs High to low", "High to high vs Low to low", 
    "High to low vs High to high", "Low to high vs Low to low"), a = c("Low to high", "High to high", "High to low", 
    "Low to high"), b = c("High to low", "Low to low", "High to high", "Low to low"))

supp_ols <- function(X, Y, extra = 0L) {
    q <- qr(X, tol = 1e-08)
    keep <- sort(q$pivot[seq_len(q$rank)])
    X <- X[, keep, drop = FALSE]
    q <- qr(X)
    b <- qr.coef(q, Y)
    if (is.null(dim(b))) 
        b <- matrix(b, ncol = 1, dimnames = list(colnames(X), colnames(Y)))
    V <- chol2inv(qr.R(q))
    dimnames(V) <- list(colnames(X), colnames(X))
    df <- nrow(X) - ncol(X) - extra
    stopifnot(df > 0)
    e <- Y - X %*% b
    list(b = b, V = V, s2 = colSums(e * e)/df, df = df, rank = ncol(X))
}

supp_contrast <- function(m, C, features, layer, n) {
    C <- C[, rownames(m$b), drop = FALSE]
    b <- C %*% m$b
    se <- sqrt(outer(diag(C %*% m$V %*% t(C)), m$s2))
    rbindlist(lapply(seq_len(nrow(C)), function(i) data.table(layer, feature = features, comparison = rownames(C)[i], 
        n, beta = as.numeric(b[i, ]), se = as.numeric(se[i, ]), lower = as.numeric(b[i, ] - qt(0.975, m$df) * se[i, 
            ]), upper = as.numeric(b[i, ] + qt(0.975, m$df) * se[i, ]), p = as.numeric(2 * stats::pt(-abs(b[i, 
            ]/se[i, ]), m$df)))))
}


# 8. Estimate TDI-group contrasts
supp_omics <- function() {
    layer <- "Protein"
    message("Supplement: four TDI-group ", layer, " comparisons")
    d <- med_data(layer, followup = FALSE)
    x <- d$ph
    om <- d$om
    info <- d$info
    f <- d$features
    x[, `:=`(tdi4, factor(as.character(tdi4), levels = tdi_lev))]
    pairs <- supp_pairs()
    results <- profiles <- list()
    cv <- setdiff(med_rhs(layer, FALSE)$dtdi, c("md_z", "dtdi_z"))
    X <- model.matrix(reformulate(c("tdi4", "md_z", cv)), droplevels(as.data.frame(x)))
    C <- matrix(0, nrow(pairs), ncol(X), dimnames = list(pairs$comparison, colnames(X)))
    for (i in seq_len(nrow(pairs))) {
        if (pairs$a[i] != tdi_lev[1]) 
            C[i, paste0("tdi4", pairs$a[i])] <- 1
        if (pairs$b[i] != tdi_lev[1]) 
            C[i, paste0("tdi4", pairs$b[i])] <- -1
    }
    for (g in split(f, ceiling(seq_along(f)/120))) {
        Y <- apply(as.matrix(om[, ..g]), 2, z)
        if (is.null(dim(Y))) 
            Y <- matrix(Y, ncol = 1, dimnames = list(NULL, g))
        results[[length(results) + 1]] <- supp_contrast(supp_ols(X, Y), C, g, layer, nrow(x))
        profiles[[length(profiles) + 1]] <- rbindlist(lapply(tdi_lev, function(v) {
            yy <- Y[x$tdi4 == v, , drop = FALSE]
            data.table(layer = layer, feature = g, group = v, n = nrow(yy), mean = colMeans(yy), sd = apply(yy, 
                2, sd), se = apply(yy, 2, sd)/sqrt(nrow(yy)))
        }))
        message("  ", layer, " ", match(tail(g, 1), f), "/", length(f))
    }
    a <- rbindlist(results)
    a[, `:=`(bonferroni, pmin(1, p * 4 * length(f)))]
    a <- merge(a, info, by.x = "feature", by.y = "Feature ID", all.x = TRUE, sort = FALSE)
    means <- merge(rbindlist(profiles), info, by.x = "feature", by.y = "Feature ID", all.x = TRUE, sort = FALSE)
    list(contrasts = a, profiles = means, counts = x[, .(N = .N), by = tdi4], methods = list(groups = "Four prespecified linear contrasts", 
        correction = paste("Bonferroni across four contrasts and", length(f), "features"), profiles = "Standardized group means"))
}

supp_pattern <- function(r) {
    a <- r$omics$Protein
    features <- unique(a$contrasts[bonferroni < 0.05, feature])
    d <- dcast(a$profiles[feature %in% features], feature ~ group, value.var = "mean")
    M <- as.matrix(d[, ..tdi_lev])
    rownames(M) <- d$feature
    Z <- t(scale(t(M)))
    set.seed(7)
    fits <- lapply(2:6, function(k) kmeans(Z, k, nstart = 100, iter.max = 100))
    ds <- dist(Z)
    score <- vapply(fits, function(x) mean(cluster::silhouette(x$cluster, ds)[, "sil_width"]), 0)
    fit <- fits[[which.max(score)]]
    order <- as.integer(names(sort(table(fit$cluster), decreasing = TRUE)))
    map <- data.table(feature = rownames(Z), pattern = match(fit$cluster, order))
    profiles <- as.data.table(Z, keep.rownames = "feature")
    profiles <- merge(melt(profiles, id.vars = "feature", variable.name = "group", value.name = "shape"), map, 
        by = "feature", sort = FALSE)
    list(map = map, profiles = profiles, selection = data.table(k = 2:6, silhouette = score))
}

if (sys.nframe() == 0L) {
    a <- commandArgs(TRUE)
    run_fig6_fit(fresh = TRUE, layers = if (length(a)) 
        a
    else c("Protein", "MRI"))
}

