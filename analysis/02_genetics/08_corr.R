# 1. Phenotypic and genetic correlations
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

use("assoc")

.fig2a_require <- function(packages) {
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) 
        stop("Fig2a requires: ", paste(missing, collapse = ", "))
}

fig2a_traits <- function() {
    data.frame(order = seq_len(13), trait = c("migration_distance", "migration_direction", "tdi", "income", "college", 
        "edu_years", "camsis", "isei", "siops", "work_hours", "work_stand", "work_heavy", "work_shift"), variable = c("migration_distance_rint", 
        "delta_tdi_rint", "tdi_rint", "income_rint", "college_rint", "edu.years_rint", "camsis_rint", "isei_rint", 
        "siops_rint", "work_hours_rint", "work_stand_rint", "work_heavy_rint", "work_shift_rint"), label = c("Migration distance", 
        "ΔTDI", "TDI", "Income", "College", "Education years", "CAMSIS", "ISEI", "SIOPS", "Work hours", "Work standing", 
        "Heavy work", "Shift work"), stringsAsFactors = FALSE)
}

fig2a_data <- function(project = get("root", envir = environment(fig2a_data), inherits = TRUE), phe_file = file.path(project, 
    "data/phe/move.rds"), h2_file = file.path(project, "data/ldsc/latest/all.h2.res"), rg_file = file.path(project, 
    "data/ldsc/latest/all.rg.res"), phe_data = NULL, h2_data = NULL, rg_data = NULL) {
    .fig2a_require("data.table")
    tr <- fig2a_traits()
    inputs <- c(phenotype = phe_file, heritability = h2_file, genetic_correlation = rg_file)
    supplied <- c(!is.null(phe_data), !is.null(h2_data), !is.null(rg_data))
    if (any(!file.exists(inputs[!supplied]))) 
        stop("Missing Fig2a inputs: ", paste(inputs[!supplied][!file.exists(inputs[!supplied])], collapse = ", "))
    if (is.null(phe_data)) 
        phe_data <- readRDS(phe_file)
    if (is.null(h2_data)) 
        h2_data <- data.table::fread(h2_file)
    if (is.null(rg_data)) 
        rg_data <- data.table::fread(rg_file)
    need <- function(x, cols, label) {
        missing <- setdiff(cols, names(x))
        if (length(missing)) 
            stop(label, " missing columns: ", paste(missing, collapse = ", "))
    }
    candidates <- intersect(c("rtdi_rint", "tdi_rint"), names(phe_data))
    if (!length(candidates)) 
        stop("Missing recruitment TDI RINT: rtdi_rint or tdi_rint")
    selected_tdi <- candidates[1]
    if (length(candidates) == 2L && !isTRUE(all.equal(phe_data[[candidates[1]]], phe_data[[candidates[2]]], tolerance = 1e-12, 
        check.attributes = FALSE))) 
        stop("rtdi_rint and tdi_rint differ; cannot silently change the fixed Fig2a")
    tr$variable[tr$trait == "tdi"] <- selected_tdi
    raw_tdi <- intersect(c("rtdi", "tdi_recruitment", "tdi"), names(phe_data))
    tdi_check <- data.frame(selected = selected_tdi, raw = if (length(raw_tdi)) 
        raw_tdi[1]
    else NA_character_, n_compared = NA_integer_, max_abs_difference = NA_real_, missing_mismatches = NA_integer_, 
        status = "Raw unavailable; existing transform retained")
    if (length(raw_tdi)) {
        raw <- suppressWarnings(as.numeric(phe_data[[raw_tdi[1]]]))
        expected <- rep(NA_real_, length(raw))
        ok <- is.finite(raw)
        if (sum(ok) >= 3L) 
            expected[ok] <- stats::qnorm((rank(raw[ok], ties.method = "average") - 0.5)/sum(ok))
        actual <- phe_data[[selected_tdi]]
        if (!is.numeric(actual)) 
            stop("Recruitment TDI RINT must be numeric")
        both <- is.finite(expected) & is.finite(actual)
        tdi_check$n_compared <- sum(both)
        tdi_check$max_abs_difference <- if (any(both)) 
            max(abs(actual[both] - expected[both]))
        else 0
        tdi_check$missing_mismatches <- sum(is.finite(actual) != is.finite(expected))
        if (tdi_check$max_abs_difference > 1e-12 || tdi_check$missing_mismatches > 0) 
            stop("Existing recruitment TDI RINT differs from make_data() raw transform")
        tdi_check$status <- "Validated against make_data() RINT; original values retained"
    }
    need(phe_data, tr$variable, "Phenotype data")
    need(h2_data, c("trait", "h2", "h2_se"), "Heritability data")
    need(rg_data, c("p1", "p2", "rg", "p"), "LDSC correlations")
    x <- as.data.frame(phe_data)[, tr$variable, drop = FALSE]
    if (!all(vapply(x, is.numeric, logical(1)))) 
        stop("RINT columns must be numeric")
    if (any(vapply(x, function(v) any(is.infinite(v)), logical(1)))) 
        stop("Infinite RINT values: resolve these before matching the reference plot")
    ntrait <- nrow(tr)
    r <- suppressWarnings(stats::cor(x, use = "pairwise.complete.obs", method = "pearson"))
    pn <- pp <- matrix(NA_real_, ntrait, ntrait)
    for (i in seq_len(ntrait)) for (j in seq_len(i - 1L)) {
        ok <- is.finite(x[[i]]) & is.finite(x[[j]])
        pn[i, j] <- sum(ok)
        if (sum(ok) >= 3L && stats::sd(x[[i]][ok]) > 0 && stats::sd(x[[j]][ok]) > 0) 
            pp[i, j] <- suppressWarnings(stats::cor.test(x[[i]][ok], x[[j]][ok], method = "pearson")$p.value)
    }
    tr$n_nonmissing <- vapply(x, function(v) sum(is.finite(v)), integer(1))
    n_rows <- nrow(x)
    rm(x, phe_data)
    rg <- as.data.frame(rg_data)
    gr <- gp <- gs <- gz <- matrix(NA_real_, ntrait, ntrait)
    for (k in seq_len(nrow(rg))) {
        i <- match(sub(".*/", "", rg$p1[k]), tr$trait)
        j <- match(sub(".*/", "", rg$p2[k]), tr$trait)
        if (is.na(i) || is.na(j) || i == j) 
            next
        gr[i, j] <- gr[j, i] <- rg$rg[k]
        gp[i, j] <- gp[j, i] <- rg$p[k]
        if ("se" %in% names(rg)) 
            gs[i, j] <- gs[j, i] <- rg$se[k]
        if ("z" %in% names(rg)) 
            gz[i, j] <- gz[j, i] <- rg$z[k]
    }
    d <- do.call(rbind, lapply(seq_len(ntrait), function(i) {
        do.call(rbind, lapply(seq_len(ntrait), function(j) {
            lower <- i > j
            val <- if (lower) 
                r[i, j]
            else if (i < j) 
                gr[i, j]
            else NA_real_
            pv <- if (lower) 
                pp[i, j]
            else if (i < j) 
                gp[i, j]
            else NA_real_
            data.frame(row = i, col = j, row_trait = tr$trait[i], col_trait = tr$trait[j], x = tr$label[j], y = tr$label[i], 
                val = val, p = pv, half = if (lower) 
                  "Phenotypic"
                else if (i < j) 
                  "LDSC"
                else "Blank", n = if (lower) 
                  pn[i, j]
                else NA_real_, se = if (i < j) 
                  gs[i, j]
                else NA_real_, z = if (i < j) 
                  gz[i, j]
                else NA_real_, label = if (is.finite(val)) 
                  paste0(sprintf("%.2f", val), if (!is.na(pv) && pv < 0.05) 
                    "*"
                  else "")
                else "", stringsAsFactors = FALSE)
        }))
    }))
    d$colour_value <- pmax(-1, pmin(1, d$val))
    d$magnitude <- pmin(abs(d$val), 1)
    d$circle_radius_cm <- ifelse(d$half == "LDSC", 0.6 * 0.45 * sqrt(d$magnitude), NA_real_)
    h <- as.data.frame(h2_data)
    selected <- h$trait %in% tr$trait
    if (anyDuplicated(h$trait[selected])) 
        stop("Duplicate heritability traits")
    if (!all(tr$trait %in% h$trait)) 
        stop("Missing heritability traits")
    h <- h[match(tr$trait, h$trait), c("trait", "h2", "h2_se"), drop = FALSE]
    rownames(h) <- NULL
    if (!is.numeric(h$h2) || !is.numeric(h$h2_se) || any(!is.finite(h$h2)) || any(!is.finite(h$h2_se)) || any(h$h2 < 
        0) || any(h$h2_se < 0)) 
        stop("Invalid heritability estimates or SEs")
    h$label <- tr$label
    h$text <- sprintf("%.3f", h$h2)
    h$lo_unclipped <- h$h2 - 1.96 * h$h2_se
    h$lo <- pmax(h$lo_unclipped, 0)
    h$hi <- h$h2 + 1.96 * h$h2_se
    h$row <- tr$order
    if (max(h$hi) <= 0) 
        stop("Heritability axis has zero range")
    inputs[supplied] <- "<supplied in memory>"
    hashes <- rep(NA_character_, length(inputs))
    hashes[!supplied] <- unname(tools::md5sum(inputs[!supplied]))
    structure(list(traits = tr, labs = tr$label, d = d, h = h, n_rows = n_rows, tdi_check = tdi_check, inputs = data.frame(role = names(inputs), 
        file = unname(inputs), md5 = hashes)), class = "move_fig2a_data")
}

if (sys.nframe() == 0L) saveRDS(fig2a_data(), file.path(out, "correlation.rds"))

