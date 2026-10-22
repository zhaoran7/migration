# 1. Official rural-urban assignment
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))


# 2. Read official classifications
urban_polygons <- function() {
    ref <- file.path(root, "data/ref/urban")
    ew <- sf::st_read(file.path(ref, "ew_lsoa2011.geojson"), quiet = TRUE)
    key <- fread(file.path(ref, "ew_ruc2011_lsoa.tsv"))
    ew$class <- ifelse(startsWith(key$RUC11[match(ew$LSOA11CD, key$LSOA11CD)], "Urban"), "Urban", "Rural")
    sc <- sf::st_read(file.path(ref, "scotland_datazone2011.geojson"), quiet = TRUE)
    key <- fread(file.path(ref, "scotland_ruc2020_datazone.csv"))
    sc$class <- ifelse(key$UrbanRural2fold2020[match(sc$CODE, key$DataZone)] == 1, "Urban", "Rural")
    ni <- sf::st_transform(sf::st_read(paste0("/vsizip/", file.path(ref, "ni_2015.zip")), quiet = TRUE), 27700)
    ni$class <- ifelse(ni$Band %in% LETTERS[1:5], "Urban", "Rural")
    land <- sf::st_read(file.path(root, "data/ref/ttwa.geojson"), quiet = TRUE)
    land <- land[startsWith(land$TTWA11CD, "N"), ]
    land <- sf::st_difference(sf::st_union(sf::st_make_valid(land)), sf::st_union(sf::st_make_valid(ni[ni$class == 
        "Urban", ])))
    rural <- sf::st_sf(class = "Rural", geometry = land)
    sf::st_make_valid(rbind(ew[, "class"], sc[, "class"], rural, ni[ni$class == "Urban", "class"]))
}


# 3. Assign endpoints and recruitment grid cells
urban_assign <- function(east, north, areas, grid = FALSE) {
    xy <- data.table(east, north)
    points <- unique(xy[is.finite(east) & is.finite(north)])
    if (!grid) {
        geom <- sf::st_as_sf(points, coords = c("east", "north"), crs = 27700, remove = FALSE)
        hits <- sf::st_intersects(geom, areas)
        values <- lapply(hits, function(i) unique(areas$class[i]))
        points[, `:=`(class, vapply(values, function(v) if (length(v) == 1L) 
            v
        else NA_character_, ""))]
    }
    else {
        polys <- lapply(seq_len(nrow(points)), function(i) {
            e <- points$east[i]
            n <- points$north[i]
            sf::st_polygon(list(matrix(c(e, n, e + 1000, n, e + 1000, n + 1000, e, n + 1000, e, n), ncol = 2, byrow = TRUE)))
        })
        geom <- sf::st_sf(cell = seq_len(nrow(points)), geometry = sf::st_sfc(polys, crs = 27700))
        hits <- sf::st_intersects(geom, areas)
        values <- lapply(hits, function(i) unique(areas$class[i]))
        fraction <- rep(NA_real_, nrow(points))
        pure <- which(lengths(values) == 1L)
        fraction[pure] <- vapply(values[pure], function(v) as.numeric(v == "Urban"), 0)
        mixed <- which(lengths(values) > 1L)
        if (length(mixed)) {
            pieces <- suppressWarnings(sf::st_intersection(geom[mixed, ], areas[unique(unlist(hits[mixed])), ]))
            pieces$area <- as.numeric(sf::st_area(pieces))
            share <- as.data.table(sf::st_drop_geometry(pieces))[, .(urban = sum(area[class == "Urban"])/sum(area)), 
                by = cell]
            fraction[share$cell] <- share$urban
        }
        points[, `:=`(class, fifelse(is.finite(fraction), ifelse(fraction > 0.5, "Urban", "Rural"), NA_character_))]
    }
    points[xy, on = .(east, north), class]
}

