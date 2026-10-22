# 1. Individual PGS associations
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

use("assoc")

prs_fit <- function() {
    x <- direction_data()
    gridfile <- file.path(root, "data/ref/population_2011.asc")
    header <- strsplit(trimws(readLines(gridfile, n = 6)), "[[:space:]]+")
    h <- setNames(vapply(header, function(v) as.numeric(v[2]), 0), tolower(vapply(header, "[", "", 1)))
    grid <- matrix(scan(gridfile, skip = 6, quiet = TRUE), nrow = h["nrows"], byrow = TRUE)
    density <- function(e, n) {
        j <- floor((e - h["xllcorner"])/h["cellsize"]) + 1L
        i <- h["nrows"] - floor((n - h["yllcorner"])/h["cellsize"])
        ok <- is.finite(i) & is.finite(j) & i >= 1 & i <= nrow(grid) & j >= 1 & j <= ncol(grid)
        v <- rep(NA_real_, length(e))
        v[ok] <- grid[cbind(i[ok], j[ok])]
        v[v < 0 | v == h["nodata_value"]] <- NA_real_
        v/(h["cellsize"]/1000)^2
    }
    x[, `:=`(density_birth = density(birth_east, birth_north), density_home = density(home_east, home_north))]
    pc <- as.data.table(readRDS(file.path(phe, "all.rds")))[, c("eid", paste0("PC", 6:40)), with = FALSE]
    pc[, `:=`(eid, as.character(eid))]
    p <- as.data.table(readRDS(file.path(phe, "prs.rds")))
    p[, `:=`(eid, as.character(eid))]
    x <- merge(merge(x, pc, by = "eid"), p, by = "eid")
    meta <- fread(file.path(code, "reference/pgs.tsv"))
    cv <- c("age", "age2", "sex", "center", "gen_batch", paste0("PC", 1:40))
    x[, `:=`(age2, as.numeric(age)^2)]
    req <- c(cv, meta$score, "md", "dtdi", "density_birth", "density_home", "density_shift", "tdi4", "urban_route", 
        "cross_ttwa")
    stopifnot(all(req %in% names(x)))
    x <- x[complete.cases(x[, ..req])]
    quart <- function(v) factor(cut(v, quantile(v, seq(0, 1, 0.25)), include.lowest = TRUE, labels = paste0("Q", 
        1:4)), levels = paste0("Q", 1:4))
    x[, `:=`(md_score = z(rint(md)), dtdi_score = z(rint(dtdi)), birth_score = z(log1p(density_birth)), home_score = z(log1p(density_home)), 
        change_score = z(density_shift), md_q = quart(md), dtdi_q = quart(dtdi), density_birth_q = quart(density_birth), 
        density_home_q = quart(density_home), density_shift_q = quart(density_shift), urban4 = factor(tolower(as.character(urban_route)), 
            levels = tolower(urban_levels)), ttwa2 = factor(cross_ttwa, levels = 0:1, labels = c("Within TTWA", 
            "Cross TTWA")))]
    x[, `:=`(urban4, factor(c("Urban to urban", "Urban to rural", "Rural to urban", "Rural to rural")[as.integer(urban4)], 
        levels = urban_levels))]
    groups <- c("md_q", "dtdi_q", "tdi4", "urban4", "density_birth_q", "density_home_q", "density_shift_q", "ttwa2")
    continuous <- c(md_q = "md_score", dtdi_q = "dtdi_score", density_birth_q = "birth_score", density_home_q = "home_score", 
        density_shift_q = "change_score")
    bars <- list()
    associations <- list()
    for (g in groups) for (s in meta$score) {
        a <- x[, .(value = get(s), group = factor(get(g)))]
        m <- lm(value ~ group, a)
        pg <- anova(m)[1, "Pr(>F)"]
        bars[[length(bars) + 1L]] <- a[, .(n = .N, estimate = mean(value), se = sd(value)/sqrt(.N)), by = group][, 
            `:=`(grouping = g, prs = s, p_group = pg)]
        if (g %in% names(continuous)) {
            d <- copy(x)
            d[, `:=`(score, get(s))]
            fit <- lm(reformulate(c("score", cv), continuous[[g]]), d)
            sm <- coef(summary(fit))["score", ]
            associations[[length(associations) + 1L]] <- data.table(grouping = g, prs = s, group = "Per SD PRS", 
                reference = NA_character_, n = nobs(fit), beta = sm[1], se = sm[2], p = sm[4], model = "Standardized migration phenotype ~ standardized PRS + covariates")
        }
        else {
            d <- copy(x)
            d[, `:=`(value = get(s), group = factor(get(g)))]
            fit <- lm(reformulate(c("group", cv), "value"), d)
            sm <- coef(summary(fit))
            lev <- levels(d$group)
            sm <- sm[paste0("group", lev[-1]), , drop = FALSE]
            associations[[length(associations) + 1L]] <- data.table(grouping = g, prs = s, group = lev[-1], reference = lev[1], 
                n = nobs(fit), beta = sm[, 1], se = sm[, 2], p = sm[, 4], model = "Standardized PRS ~ migration group + covariates")
        }
    }
    b <- rbindlist(bars)
    a <- rbindlist(associations)
    a[, `:=`(lower = beta - 1.96 * se, upper = beta + 1.96 * se, q = p.adjust(p, "BH"))]
    tests <- unique(b[, .(grouping, prs, p_group)])
    tests[, `:=`(q_group, p.adjust(p_group, "BH"))]
    b[tests, on = .(grouping, prs), `:=`(q_group, i.q_group)]
    labels <- c("Education", "Income", "Mortality proxy", "CVD proxy", "CAD", "Stroke", "HF", "AF", "T2D", "Depression", 
        "Anxiety", "Dementia")
    b[, `:=`(prs_label, factor(labels[match(prs, meta$score)], levels = labels))]
    r <- list(bars = b, associations = a, definitions = meta, display = c("md_q", "dtdi_q", "tdi4", "urban4", "density_shift_q", 
        "ttwa2"))
    saveRDS(r, file.path(out, "prs.rds"), compress = "gzip")
    save_xlsx(list(Groups = b, Associations = a, Sources = meta), file.path(out, "prs.xlsx"))
    r
}

if (sys.nframe() == 0L) prs_fit()

