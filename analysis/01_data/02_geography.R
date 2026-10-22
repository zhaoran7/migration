# 1. TTWA linkage and geographic summaries
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

use(c("ukb"))

map_require_columns <- function(x, columns) {
    missing <- setdiff(columns, names(x))
    if (length(missing)) 
        stop("Missing columns: ", paste(missing, collapse = ", "))
    invisible(x)
}

fig1_map_labels <- function() {
    data.table(name = c("London", "Birmingham", "Manchester", "Liverpool", "Leeds", "Sheffield", "Newcastle", "Glasgow", 
        "Edinburgh", "Bristol", "Cardiff", "Nottingham", "Leicester", "Oxford", "Cambridge", "Reading", "Southampton", 
        "Brighton", "Norwich", "Plymouth", "Aberdeen", "Swansea", "Belfast"), dx = c(28000, -38000, -40000, -42000, 
        26000, 30000, 28000, -42000, 30000, -38000, -42000, 30000, 30000, -32000, 33000, -34000, 26000, 30000, 
        36000, -42000, 32000, -42000, -42000), dy = c(-28000, -12000, 24000, -4000, 22000, -18000, 26000, 20000, 
        -14000, -22000, -18000, -4000, -17000, -24000, 17000, -28000, -28000, -24000, 18000, -18000, 16000, 16000, 
        0))
}

map_ring <- function(ring) {
    x <- do.call(rbind, lapply(ring, function(z) as.numeric(unlist(z))[1:2]))
    colnames(x) <- c("x", "y")
    x
}

map_polys <- function(g) {
    if (g$type == "Polygon") 
        return(list(lapply(g$coordinates, map_ring)))
    if (g$type != "MultiPolygon") 
        stop("Expected Polygon or MultiPolygon geometry")
    lapply(g$coordinates, function(x) lapply(x, map_ring))
}

map_centroid <- function(ps) {
    z <- rbindlist(lapply(ps, function(p) {
        m <- p[[1]]
        x <- m[, 1]
        y <- m[, 2]
        if (x[1] != x[length(x)] || y[1] != y[length(y)]) {
            x <- c(x, x[1])
            y <- c(y, y[1])
        }
        cr <- x[-length(x)] * y[-1] - x[-1] * y[-length(y)]
        a <- sum(cr)/2
        data.table(a = abs(a), cx = sum((x[-length(x)] + x[-1]) * cr)/(6 * a), cy = sum((y[-length(y)] + y[-1]) * 
            cr)/(6 * a))
    }))
    c(cx = weighted.mean(z$cx, z$a), cy = weighted.mean(z$cy, z$a), area = sum(z$a))
}


# 2. Read TTWA boundaries
read_ttwa <- function(geography = file.path(root, "data/ref/ttwa.geojson")) {
    geo <- if (is.list(geography)) 
        geography
    else jsonlite::fromJSON(geography, simplifyVector = FALSE)
    if (is.null(geo$features) || !length(geo$features)) 
        stop("No TTWA features found")
    if (isTRUE(geo$exceededTransferLimit)) 
        stop("Incomplete TTWA service response")
    crs <- geo$crs$properties$name
    if (!is.null(crs) && !grepl("27700", crs, fixed = TRUE)) {
        stop("TTWA geometry must be EPSG:27700; found ", crs)
    }
    ft <- lapply(geo$features, function(f) {
        if (is.null(f$properties$TTWA11CD) || is.null(f$properties$TTWA11NM)) {
            stop("TTWA11CD and TTWA11NM properties are required")
        }
        ps <- map_polys(f$geometry)
        xy <- do.call(rbind, lapply(ps, function(p) do.call(rbind, p)))
        cc <- map_centroid(ps)
        list(code = f$properties$TTWA11CD, name = f$properties$TTWA11NM, polys = ps, bbox = c(range(xy[, 1]), range(xy[, 
            2])), cx = cc["cx"], cy = cc["cy"], area = cc["area"])
    })
    if (anyDuplicated(vapply(ft, `[[`, character(1), "code"))) 
        stop("Duplicate TTWA codes")
    ft
}

map_pip <- function(x, y, ring) {
    rx <- ring[, 1]
    ry <- ring[, 2]
    inside <- rep(FALSE, length(x))
    j <- length(rx)
    for (i in seq_along(rx)) {
        hit <- ((ry[i] > y) != (ry[j] > y)) & (x < (rx[j] - rx[i]) * (y - ry[i])/(ry[j] - ry[i] + 1e-12) + rx[i])
        hit[is.na(hit)] <- FALSE
        inside <- xor(inside, hit)
        j <- i
    }
    inside
}


# 3. Assign endpoint areas
assign_ttwa <- function(x, y, ft) {
    if (length(x) != length(y)) 
        stop("Coordinate vectors must have equal lengths")
    code <- name <- rep(NA_character_, length(x))
    ok <- is.finite(x) & is.finite(y)
    for (f in ft) {
        cand0 <- ok & is.na(code) & x >= f$bbox[1] & x <= f$bbox[2] & y >= f$bbox[3] & y <= f$bbox[4]
        if (!any(cand0)) 
            next
        hit <- rep(FALSE, length(x))
        for (p in f$polys) {
            o <- p[[1]]
            bb <- c(range(o[, 1]), range(o[, 2]))
            cand <- cand0 & x >= bb[1] & x <= bb[2] & y >= bb[3] & y <= bb[4]
            if (!any(cand)) 
                next
            ins <- map_pip(x[cand], y[cand], o)
            if (length(p) > 1 && any(ins)) 
                for (h in p[-1]) ins <- ins & !map_pip(x[cand], y[cand], h)
            hit[which(cand)[ins]] <- TRUE
        }
        code[hit] <- f$code
        name[hit] <- f$name
    }
    data.table(code, name)
}

map_assign_unique <- function(x, y, ft) {
    xy <- data.table(x = as.numeric(x), y = as.numeric(y))
    keys <- unique(xy[is.finite(x) & is.finite(y)])
    a <- assign_ttwa(keys$x, keys$y, ft)
    keys[, `:=`(code = a$code, name = a$name)]
    keys[xy, on = .(x, y), .(code, name)]
}

map_coverage <- function(x, origin, destination) {
    data.table(metric = c("input_rows", "finite_birth_coordinates", "finite_home_coordinates", "origin_assigned", 
        "destination_assigned", "both_assigned", "same_ttwa", "cross_ttwa", "unassigned_either"), n = c(nrow(x), 
        sum(is.finite(as.numeric(x$birth_east)) & is.finite(as.numeric(x$birth_north))), sum(is.finite(as.numeric(x$home_east)) & 
            is.finite(as.numeric(x$home_north))), sum(!is.na(origin)), sum(!is.na(destination)), sum(!is.na(origin) & 
            !is.na(destination)), sum(origin == destination, na.rm = TRUE), sum(origin != destination, na.rm = TRUE), 
        sum(is.na(origin) | is.na(destination))))
}

fig1_map_data <- function(all_file = file.path(root, "data/phe/all.rds"), ids = read_ids("migration_distance"), 
    geography = file.path(root, "data/ref/ttwa.geojson"), data = NULL) {
    default_ids <- missing(ids)
    ft <- read_ttwa(geography)
    labels <- fig1_map_labels()
    poly <- rbindlist(lapply(ft, function(f) {
        rbindlist(lapply(seq_along(f$polys), function(i) {
            rbindlist(lapply(seq_along(f$polys[[i]]), function(j) {
                m <- f$polys[[i]][[j]]
                data.table(code = f$code, name = f$name, group = paste(f$code, i, j, sep = "_"), x = m[, 1], y = m[, 
                  2])
            }))
        }))
    }))
    cen <- rbindlist(lapply(ft, function(f) {
        data.table(code = f$code, name = f$name, cx = f$cx, cy = f$cy, area = f$area, is_core = f$name %in% labels$name)
    }))
    x <- if (is.null(data)) 
        as.data.table(readRDS(all_file))
    else as.data.table(data)
    map_require_columns(x, c("eid", "birth_east", "birth_north", "home_east", "home_north", "move", "tdi"))
    ids <- unique(as.character(ids))
    d <- x[as.character(eid) %chin% ids, .(birth_east = as.numeric(birth_east), birth_north = as.numeric(birth_north), 
        home_east = as.numeric(home_east), home_north = as.numeric(home_north), move = as.numeric(move), tdi = as.numeric(tdi))]
    a <- map_assign_unique(d$birth_east, d$birth_north, ft)
    b <- map_assign_unique(d$home_east, d$home_north, ft)
    d[, `:=`(from = a$code, from_name = a$name, to = b$code, to_name = b$name)]
    tdi_area <- d[!is.na(to) & is.finite(tdi), .(tdi_med = median(tdi), tdi_n = .N), by = .(code = to)]
    poly <- merge(poly, tdi_area[, .(code, tdi_med)], by = "code", all.x = TRUE, sort = FALSE)
    od <- d[!is.na(move) & !is.na(from) & !is.na(to) & from != to, .(N = .N), by = .(from, to)]
    od <- merge(od, cen[, .(from = code, from_name = name, x0 = cx, y0 = cy)], by = "from")
    od <- merge(od, cen[, .(to = code, to_name = name, x1 = cx, y1 = cy)], by = "to")
    nodes <- merge(cen, od[, .(in_n = sum(N)), by = .(code = to)], by = "code", all.x = TRUE)
    nodes <- merge(nodes, od[, .(out_n = sum(N)), by = .(code = from)], by = "code", all.x = TRUE)
    nodes[is.na(in_n), `:=`(in_n, 0)]
    nodes[is.na(out_n), `:=`(out_n, 0)]
    nodes[, `:=`(total = in_n + out_n, net = in_n - out_n)]
    edges <- od[from %in% nodes[is_core == TRUE, code] & to %in% nodes[is_core == TRUE, code] & N >= 100]
    edges[, `:=`(id, .I)]
    areas <- merge(nodes, tdi_area, by = "code", all.x = TRUE, sort = FALSE)
    areas[is.na(tdi_n), `:=`(tdi_n, 0L)]
    areas[, `:=`(core_pie = is_core & total > 0, other_point = !is_core & total > 0)]
    flows <- copy(od)
    flows[, `:=`(plotted, from %in% nodes[is_core == TRUE, code] & to %in% nodes[is_core == TRUE, code] & N >= 
        100)]
    coverage <- rbindlist(list(map_coverage(d, d$from, d$to), data.table(metric = c("saved_gwas_ids", "gwas_ids_found", 
        "gwas_ids_absent", "nonmissing_move", "inter_ttwa_migrants_with_move", "all_flow_pairs", "plotted_core_flow_pairs", 
        "plotted_core_migrants", "areas_in_reference", "core_names_in_reference", "core_labels_with_flow"), n = c(length(ids), 
        sum(ids %chin% as.character(x$eid)), sum(!ids %chin% as.character(x$eid)), sum(!is.na(d$move)), sum(od$N), 
        nrow(od), nrow(edges), sum(edges$N), nrow(cen), sum(labels$name %chin% cen$name), nodes[is_core == TRUE & 
            total > 0, .N]))))
    label_coverage <- merge(labels, nodes[, .(name, code, total)], by = "name", all.x = TRUE, sort = FALSE)
    label_coverage[, `:=`(in_reference = !is.na(code), drawn = !is.na(total) & total > 0)]
    provenance <- data.table(field = c("all_file", "ids", "geography", "crs", "flow_filter", "background", "labels"), 
        value = c(if (is.null(data)) all_file else "supplied data", if (default_ids) file.path(root, "data/phe/migration_distance.ids") else "explicit ids argument", 
            if (is.character(geography)) geography else "supplied GeoJSON object", "EPSG:27700", "!is.na(move); both assigned; from != to; displayed core pairs N >= 100", 
            "median finite recruitment tdi by destination; all selected GWAS rows", "exact TTWA11NM matches; original 23 core labels and offsets"))
    message("map: ", nrow(d), " distance-GWAS samples; ", sum(edges$N), " plotted core-to-core migrants")
    source <- list(nodes = areas, flows = flows, coverage = coverage, label_coverage = label_coverage, provenance = provenance)
    plotting <- list(poly = poly, nodes = nodes, edges = edges, labels = labels)
    structure(list(source = source, plotting = plotting, poly = poly, nodes = nodes, edges = edges, areas = areas, 
        flows = flows, labels = labels, label_coverage = label_coverage, coverage = coverage, provenance = provenance), 
        class = "fig1_map_data")
}


# 4. Save participant linkage
make_ttwa <- function(all_file = file.path(root, "data/phe/all.rds"), geography = file.path(root, "data/ref/ttwa.geojson"), 
    data = NULL) {
    x <- if (is.null(data)) 
        as.data.table(readRDS(all_file))
    else copy(as.data.table(data))
    map_require_columns(x, c("eid", "birth_east", "birth_north", "home_east", "home_north", "ethnic.c", "age", 
        "sex", "center", "gen_batch", "tdi.birth", paste0("PC", 1:5)))
    ft <- read_ttwa(geography)
    a <- map_assign_unique(x$birth_east, x$birth_north, ft)
    b <- map_assign_unique(x$home_east, x$home_north, ft)
    x[, `:=`(origin_ttwa = a$code, destination_ttwa = b$code)]
    x[, `:=`(cross_ttwa, as.integer(origin_ttwa != destination_ttwa))]
    assignments <- x[, .(FID = eid, IID = eid, origin_ttwa, destination_ttwa, cross_ttwa)]
    phenotype <- x[as.character(ethnic.c) == "White", .(FID = eid, IID = eid, cross_ttwa, origin_ttwa, destination_ttwa, 
        age = as.numeric(age), age2 = as.numeric(age)^2, sex = as.numeric(sex), center, gen_batch = fifelse(is.na(gen_batch), 
            NA_character_, paste0("b", gen_batch)), tdi_birth = as.numeric(tdi.birth), birth_east = as.numeric(birth_east), 
        birth_north = as.numeric(birth_north), PC1, PC2, PC3, PC4, PC5)]
    covs <- c("age", "age2", "sex", "center", "gen_batch", "tdi_birth", "birth_east", "birth_north", paste0("PC", 
        1:5))
    eligible <- !is.na(phenotype$cross_ttwa) & complete.cases(phenotype[, ..covs])
    coverage <- rbindlist(list(map_coverage(x, x$origin_ttwa, x$destination_ttwa), data.table(metric = c("white_rows", 
        "white_assigned_both", "white_complete_covariates", "phenotype_rows", "phenotype_controls", "phenotype_cases", 
        "phenotype_origin_ttwas"), n = c(nrow(phenotype), sum(!is.na(phenotype$cross_ttwa)), sum(complete.cases(phenotype[, 
        ..covs])), sum(eligible), sum(phenotype$cross_ttwa[eligible] == 0L), sum(phenotype$cross_ttwa[eligible] == 
        1L), uniqueN(phenotype$origin_ttwa[eligible])))))
    phenotype <- phenotype[eligible]
    phenotype[, `:=`(birth_east = as.numeric(scale(birth_east)), birth_north = as.numeric(scale(birth_north)))]
    phenotype[, `:=`(birth_east2 = birth_east^2, birth_north2 = birth_north^2, birth_cross = birth_east * birth_north)]
    message("TTWA phenotype: ", nrow(phenotype), " rows; ", uniqueN(phenotype$origin_ttwa), " origin TTWAs")
    list(phenotype = phenotype, assignments = assignments, coverage = coverage, counts = phenotype[, .N, by = cross_ttwa], 
        geography = geography, covariates = covs)
}

ttwa_table <- function(fresh = FALSE) {
    file <- file.path(phe, "ttwa.rds")
    if (!fresh && file.exists(file)) 
        return(readRDS(file))
    legacy <- file.path(phe, "ukb.pheno")
    if (!fresh && file.exists(legacy)) {
        x <- fread(legacy, select = c("IID", "cross_ttwa", "origin_ttwa", "destination_ttwa"))
    }
    else x <- make_ttwa()$phenotype[, .(IID, cross_ttwa, origin_ttwa, destination_ttwa)]
    setnames(x, c("IID", "origin_ttwa", "destination_ttwa"), c("eid", "origin", "destination"))
    x[, `:=`(eid, as.character(eid))]
    if (fresh || !file.exists(legacy)) 
        saveRDS(x, file, compress = "gzip")
    x
}


# 5. Summarize exposure distributions
fig1_description_data <- function() {
    fields <- c("eid", "home_east", "home_north", "birth_east", "birth_north", "tdi.birth", "tdi")
    x <- as.data.table(readRDS(file.path(phe, "all.rds")))[, ..fields]
    x[, `:=`(eid, as.character(eid))]
    md_ids <- read_ids("migration_distance")
    stopifnot(!anyDuplicated(x$eid), all(md_ids %in% x$eid))
    x <- x[eid %in% md_ids]
    x[, `:=`(md, sqrt((home_east - birth_east)^2 + (home_north - birth_north)^2)/1000)]
    stopifnot(all(is.finite(x$md)), all(x$md > 0))
    base <- as.data.table(readRDS(file.path(phe, "move.rds")))
    cutoff <- median(base$rtdi, na.rm = TRUE)
    rm(base)
    x[, `:=`(tdi4, fcase(tdi.birth <= cutoff & tdi <= cutoff, "Low to low", tdi.birth > cutoff & tdi <= cutoff, 
        "High to low", tdi.birth <= cutoff & tdi > cutoff, "Low to high", tdi.birth > cutoff & tdi > cutoff, "High to high"))]
    u <- as.data.table(readRDS(file.path(root, "data/ref/urban/ukb.rds")))
    u[, `:=`(eid, as.character(eid))]
    t <- fread(file.path(phe, "ukb.pheno"), select = c("IID", "cross_ttwa"))
    setnames(t, "IID", "eid")
    t[, `:=`(eid, as.character(eid))]
    stopifnot(!anyDuplicated(u$eid), !anyDuplicated(t$eid))
    x <- merge(x, u[, .(eid, urban_route)], by = "eid", all.x = TRUE, sort = FALSE)
    x <- merge(x, t, by = "eid", all.x = TRUE, sort = FALSE)
    x[, `:=`(ttwa, fifelse(cross_ttwa == 1, "Cross TTWA", fifelse(cross_ttwa == 0, "Within TTWA", NA_character_)))]
    file <- file.path(root, "data/ref/population_2011.asc")
    header <- strsplit(trimws(readLines(file, n = 6)), "[[:space:]]+")
    h <- setNames(vapply(header, function(v) as.numeric(v[2]), 0), tolower(vapply(header, "[", "", 1)))
    g <- matrix(scan(file, skip = 6, quiet = TRUE), nrow = h["nrows"], byrow = TRUE)
    population <- function(east, north) {
        col <- floor((east - h["xllcorner"])/h["cellsize"]) + 1L
        row <- h["nrows"] - floor((north - h["yllcorner"])/h["cellsize"])
        ok <- is.finite(row) & is.finite(col) & row >= 1 & row <= nrow(g) & col >= 1 & col <= ncol(g)
        v <- rep(NA_real_, length(east))
        v[ok] <- g[cbind(row[ok], col[ok])]
        v[v < 0 | v == h["nodata_value"]] <- NA_real_
        v/(h["cellsize"]/1000)^2
    }
    x[, `:=`(density_shift, log1p(population(home_east, home_north)) - log1p(population(birth_east, birth_north)))]
    composition <- rbindlist(lapply(c("tdi4", "urban_route", "ttwa"), function(v) {
        d <- x[!is.na(get(v)), .(N = .N), by = .(Route = get(v))]
        d[, `:=`(c("Variable", "Share"), list(v, N/sum(N)))]
        d
    }))
    density_rows <- function(v, section, transformation) {
        v <- v[is.finite(v)]
        d <- stats::density(v, from = min(v), to = max(v), n = 1024)
        data.table(section = section, transformation = transformation, n = length(v), bandwidth = d$bw, x = d$x, 
            density = d$y)
    }
    density <- rbindlist(list(density_rows(rint(x$md), "MD RINT", "RINT within distance GWAS IDs"), density_rows(x$density_shift, 
        "Density change", "log1p recruitment density minus log1p birthplace density")))
    coverage <- data.table(metric = c("MD descriptive N", "TDI descriptive N", "TTWA observed N", "Urban route observed N", 
        "Population density observed N"), n = c(nrow(x), sum(!is.na(x$tdi4)), sum(!is.na(x$ttwa)), sum(!is.na(x$urban_route)), 
        sum(is.finite(x$density_shift))))
    list(density = density, composition = composition, coverage = coverage, definition = data.table(parameter = "Shared TDI threshold", 
        value = cutoff))
}

ed2_data <- function() {
    proj <- file.path(dirname(dirname(R.home())), "share/proj")
    if (!nzchar(Sys.getenv("PROJ_DATA")) && file.exists(file.path(proj, "proj.db"))) 
        Sys.setenv(PROJ_DATA = proj)
    if (!requireNamespace("sf", quietly = TRUE)) 
        stop("ED2 requires the sf package")
    ref <- file.path(root, "data/ref")
    areas <- sf::st_read(file.path(ref, "ttwa.geojson"), quiet = TRUE)
    stopifnot(sf::st_crs(areas)$epsg == 27700, !anyDuplicated(areas$TTWA11CD))
    areas <- sf::st_make_valid(areas)
    fields <- c("eid", "birth_east", "birth_north", "home_east", "home_north")
    x <- as.data.table(readRDS(file.path(phe, "all.rds")))[, ..fields]
    ids <- read_ids("migration_distance")
    x <- x[as.character(eid) %chin% ids]
    stopifnot(nrow(x) == length(ids), !anyDuplicated(x$eid))
    message("ED2: individual birthplace and residence points")
    a <- x[is.finite(birth_east) & is.finite(birth_north), .(easting = birth_east, northing = birth_north)]
    b <- x[is.finite(home_east) & is.finite(home_north), .(easting = home_east, northing = home_north)]
    old <- fig1_map_data()$areas
    c <- old[, .(code, area = name, n = as.integer(tdi_n), tdi = tdi_med)]
    stopifnot(nrow(c) == nrow(areas), setequal(c$code, areas$TTWA11CD))
    coverage <- data.table(metric = c("Selected MD GWAS sample", "Birthplace points", "Recruitment residence points", 
        "TDI areas copied from Fig1c"), n = c(nrow(x), nrow(a), nrow(b), sum(is.finite(c$tdi))))
    rm(x)
    invisible(gc())
    message("ED2: rural and urban reference geography")
    ew <- sf::st_read(file.path(ref, "urban/ew_lsoa2011.geojson"), quiet = TRUE)
    look <- fread(file.path(ref, "urban/ew_ruc2011_lsoa.tsv"))
    ew$class <- ifelse(startsWith(look$RUC11[match(ew$LSOA11CD, look$LSOA11CD)], "Urban"), "Urban", "Rural")
    sc <- sf::st_read(file.path(ref, "urban/scotland_datazone2011.geojson"), quiet = TRUE)
    look <- fread(file.path(ref, "urban/scotland_ruc2020_datazone.csv"))
    sc$class <- ifelse(look$UrbanRural2fold2020[match(sc$CODE, look$DataZone)] == 1, "Urban", "Rural")
    ni <- sf::st_read(paste0("/vsizip/", file.path(ref, "urban/ni_2015.zip")), quiet = TRUE)
    ni <- sf::st_transform(ni, 27700)
    ni$class <- ifelse(ni$Band %in% LETTERS[1:5], "Urban", "Rural")
    stopifnot(!anyNA(ew$class), !anyNA(sc$class), !anyNA(ni$class))
    ni_land <- areas[startsWith(areas$TTWA11CD, "N"), ]
    ni_land$class <- "Rural"
    urban_source <- rbindlist(list(data.table(code = ew$LSOA11CD, classification = ew$class, country = "England / Wales", 
        year = 2011L), data.table(code = sc$CODE, classification = sc$class, country = "Scotland", year = 2020L), 
        data.table(code = as.character(seq_len(nrow(ni))), classification = ni$class, country = "Northern Ireland", 
            year = 2015L)))
    urban <- rbind(ew[, "class"], sc[, "class"], ni_land[, "class"], ni[ni$class == "Urban", "class"])
    urban <- sf::st_simplify(sf::st_make_valid(urban), dTolerance = 100, preserveTopology = TRUE)
    urban <- urban[!sf::st_is_empty(urban), ]
    message("ED2: 2011 population grid")
    file <- file.path(ref, "population_2011.asc")
    header <- strsplit(trimws(readLines(file, n = 6)), "[[:space:]]+")
    h <- setNames(vapply(header, function(v) as.numeric(v[2]), 0), tolower(vapply(header, "[", "", 1)))
    g <- matrix(scan(file, skip = 6, quiet = TRUE), nrow = h["nrows"], byrow = TRUE)
    ix <- which(g >= 0 & g != h["nodata_value"], arr.ind = TRUE)
    pop <- data.table(easting = h["xllcorner"] + (ix[, 2] - 0.5) * h["cellsize"], northing = h["yllcorner"] + (h["nrows"] - 
        ix[, 1] + 0.5) * h["cellsize"], population = g[ix]/(h["cellsize"]/1000)^2)
    area_source <- data.table(code = areas$TTWA11CD, area = areas$TTWA11NM, area_km2 = as.numeric(sf::st_area(areas))/1e+06, 
        year = 2011L)
    list(areas = areas, urban = urban, pop = pop, cellsize = h["cellsize"], birth = a, home = b, tdi = c, coverage = coverage, 
        urban_source = urban_source, area_source = area_source)
}

if (sys.nframe() == 0L) {
    ttwa_table(fresh = TRUE)
}

