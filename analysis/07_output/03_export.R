# 1. Figure-specific source workbooks
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

use(c("assoc", "interaction", "corr", "region", "external", "omics", "describe", "genetic", "enrich"))


# 2. Remove participant identifiers
clean_source <- function(x) {
    if (is.data.frame(x)) {
        x <- copy(x)
        remove <- intersect(c("eid", "IID", "FID", "Source", "source_file"), names(x))
        if (length(remove)) 
            x[, `:=`((remove), NULL)]
        return(x)
    }
    if (is.list(x)) {
        x[] <- lapply(x, clean_source)
        return(x)
    }
    x
}


# 3. Collect the data displayed in each figure
figure_source <- function(id) {
    get <- function(n) readRDS(file.path(out, paste0(n, ".rds")))
    geography <- function() region_geometry(file.path(root, "data/ref/ttwa.geojson"))
    if (id == "Fig1") 
        return(list(description = fig1_description_data(), map = list(plotting = fig1_map_data()$plotting)))
    if (id == "Fig2") 
        return(list(result = get("Fig2"), a = fig2a_data(), bounds = int_bounds(), selected = int_main()))
    if (id == "Fig3") {
        r <- region_read()
        return(list(result = r[c("tests", "area", "joint", "method", "cache_provenance")], geography = geography()))
    }
    if (id == "Fig4") 
        return(list(result = get("Fig4")))
    if (id == "Fig5") 
        return(list(result = get("Fig5")))
    if (id == "Fig6") 
        return(list(result = enrich_data()))
    if (id == "ED1") 
        return(list(result = ed1_data()))
    if (id == "ED2") 
        return(list(result = ed2_data()))
    if (id == "ED3") 
        return(list(manhattan = gwas_points(), qc = ed3_qc(), comparison = ed3_comparison()))
    if (id %in% c("ED4", "ED5", "ED6")) 
        return(list(prs = get("prs")))
    if (id == "ED7") 
        return(list(result = get("ED7"), main_result = get("Fig2")))
    if (id == "ED8") {
        if (file.exists(file.path(out, "profiles.rds"))) 
            return(list(result = get("profiles")))
        r <- get("ED8")
        r$plot$state$geometry <- list(text = list(), squares = list(), regional = list())
        r$plot$state$r <- r$plot$state$r[c("description", "results", "best", "targets", "reference", "N", "SES_n")]
        r$clusters <- NULL
        return(list(result = r))
    }
    if (id == "ED9") 
        return(list(result = get("ED9"), bounds = int_bounds(), selected = int_main()))
    if (id == "ED10") 
        return(list(result = get("ED10")[c("Effects", "Curves", "Within", "Methods")]))
    s <- get("Supp")
    if (id == "S1") 
        return(list(result = s$flow))
    if (id == "S2") 
        return(list(supp = list(base = list(ses = s$base$ses)), result = get("Fig2")))
    if (id == "S3") {
        file <- file.path(out, "supp/S3_source_data.xlsx")
        r <- s$base$geo
        if ("km" %in% names(r$nearest)) {
            a <- r$nearest[, .(N = .N), by = .(X = log1p(km))][order(X)]
            a[, `:=`(Y, cumsum(N)/sum(N))]
            r$nearest <- a[, .(X, Y)]
            width <- diff(range(c(r$tdi$btdi, r$tdi$rtdi), finite = TRUE))/65
            a <- r$tdi[, .(MT = .N), by = .(ix = ceiling(btdi/width) - 1L, iy = ceiling(rtdi/width) - 1L)]
            a[, `:=`(Xmin = ix * width, Xmax = (ix + 1) * width, Ymin = iy * width, Ymax = (iy + 1) * width)]
            r$tdi <- a[, .(Xmin, Xmax, Ymin, Ymax, MT)]
        }
        else {
            r$nearest <- as.data.table(openxlsx::read.xlsx(file, "a", check.names = FALSE))
            r$tdi <- as.data.table(openxlsx::read.xlsx(file, "b", check.names = FALSE))
        }
        return(list(supp = list(base = list(geo = r), matching = s$matching)))
    }
    if (id == "S4") 
        return(list(supp = list(adjustment = s$adjustment)))
    if (id == "S5") 
        return(list(supp = list(base = list(regional = s$base$regional))))
    if (id == "S6") 
        return(list(result = readRDS(file.path(out, "supp/S6.rds")), geography = geography()))
    if (id == "S7") {
        f <- file.path(out, "supp/S7_source_data.xlsx")
        if (file.exists(f)) 
            return(setNames(lapply(letters[1:4], function(n) as.data.table(openxlsx::read.xlsx(f, n, check.names = FALSE))), 
                letters[1:4]))
        e <- s$base$external
        a <- rbindlist(lapply(names(e), function(k) copy(e[[k]]$intervals)[, .(Cohort = k, Observed.Intervals = observed, 
            N = n, Fraction = n/sum(n))]))
        b <- rbindlist(lapply(names(e), function(k) copy(e[[k]]$coverage)[order(fraction), .(Cohort = k, Recorded.Fraction = fraction, 
            Cumulative.Share = cumsum(n)/sum(n))]))
        curves <- function(field, adjust) {
            values <- unlist(lapply(e, function(x) x[[field]]$years))
            range <- range(values, finite = TRUE)
            rbindlist(lapply(names(e), function(k) {
                d <- stats::density(e[[k]][[field]]$years, adjust = adjust, from = range[1], to = range[2], n = 512, 
                  na.rm = TRUE)
                data.table(Cohort = k, Years = d$x, Density = d$y)
            }))
        }
        c <- curves("follow", 1.15)
        d <- curves("interval_length", 1.1)
        limit <- quantile(unlist(lapply(e, function(x) x$interval_length$years)), 0.99, na.rm = TRUE)
        return(list(a = a, b = b, c = c, d = d[Years <= limit]))
    }
    if (id == "S8") 
        return(list(supp = s[c("omics", "pattern", "display")]))
    if (id == "S9") {
        folder <- file.path(root, "analysis/enrich/extra")
        return(list(cells = readRDS(file.path(folder, "cells.rds")), umap = setNames(lapply(c("heart", "liver"), 
            function(a) {
                x <- fread(file.path(folder, a, "umap.tsv.gz"))
                x[, intersect(c("x", "y", "label"), names(x)), with = FALSE]
            }), c("Heart", "Liver"))))
    }
    if (id == "S10") {
        folder <- file.path(root, "analysis/enrich/spatial")
        tables <- function(part) setNames(lapply(c("cs8", "dlpfc"), function(a) setNames(lapply(c("MD", "DS"), 
            function(k) fread(file.path(folder, a, part, if (part == "cauchy_combination") 
                paste0(a, "_", k, ".Cauchy.csv.gz")
            else paste0(a, "_", k, ".csv.gz")))), c("MD", "DS"))), c("cs8", "dlpfc"))
        return(list(coordinates = setNames(lapply(c("cs8", "dlpfc"), function(a) fread(file.path(folder, paste0(a, 
            "_spots.csv")))), c("cs8", "dlpfc")), tissues = tables("cauchy_combination"), spatial = tables("spatial_ldsc")))
    }
    stop("Unknown figure: ", id)
}

if (sys.nframe() == 0L) {
    args <- commandArgs(TRUE)
    ids <- if (length(args)) 
        args
    else c(paste0("Fig", 1:6), paste0("ED", 1:10), paste0("S", 1:10))
    folder <- Sys.getenv("MOVE_SOURCE", file.path(root, "source"))
    for (id in ids) {
        message("Source: ", id)
        group <- if(startsWith(id,"Fig")) "01_Main" else if(startsWith(id,"ED")) "02_Extended_Data" else "03_Supplementary"
        destination <- file.path(folder,group,"source")
        dir.create(destination,recursive=TRUE,showWarnings=FALSE)
        if (id == "ED2") {
            use("map")
            ed2_export(destination)
            next
        }
        x <- clean_source(figure_source(id))
        source_write(x, file.path(destination, paste0(id, "_source_data.xlsx")))
        rm(x)
        invisible(gc())
    }
}

