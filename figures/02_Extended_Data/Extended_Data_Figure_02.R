# 1. Extended Data Figure 02
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

ed2_plot <- function(r) {
    gray <- "#DEDEDE"
    point_col <- rep("#01A6AC", 2)
    ttwa_fill <- "#EEE2EF"
    ttwa_line <- "#BDA0BF"
    tdi_palette <- scales::col_numeric(c("#2166AC", "#F7F7F7", "#B2182B"), c(-5, 5), na.color = "#F2F2F2")
    tdi_col <- tdi_palette(seq(-5, 5, length.out = 256))
    pop_colors <- c("#F5F9FC", "#C5DDEC", "#85B6D5", "#4C8FB9", "#235B8C", "#123252")
    pop_lim <- range(log1p(r$pop$population))
    pop_palette <- scales::col_numeric(pop_colors, pop_lim, na.color = "#F3F3F3")
    pop_col <- pop_palette(seq(pop_lim[1], pop_lim[2], length.out = 256))
    paths <- function(g) {
        g <- sf::st_cast(sf::st_geometry(g), "MULTIPOLYGON", warn = FALSE)
        lapply(g, function(p) unlist(p, recursive = FALSE))
    }
    area_paths <- paths(r$areas)
    urban_paths <- paths(r$urban)
    values <- function(d, v) d[[v]][match(r$areas$TTWA11CD, d$code)]
    tdi_fill <- scales::alpha(tdi_palette(pmax(-5, pmin(5, values(r$tdi, "tdi")))), 0.58)
    pp <- copy(r$pop)
    pp[, `:=`(colour, pop_palette(log1p(population)))]
    text <- function(s, x, y, size = 6, face = "plain", just = "left") grid::grid.text(s, grid::unit(x, "cm"), 
        grid::unit(y, "cm"), just = just, gp = grid::gpar(fontfamily = font, fontsize = size, fontface = face, 
            col = "#222222"))
    box <- function(x, y, w, h, fill) grid::grid.rect(grid::unit(x, "cm"), grid::unit(y, "cm"), grid::unit(w, "cm"), 
        grid::unit(h, "cm"), gp = grid::gpar(fill = fill, col = NA))
    shapes <- function(p, col, border = NA) {
        for (fill in unique(col)) {
            rings <- unlist(p[col == fill], recursive = FALSE)
            if (!length(rings)) 
                next
            xy <- do.call(rbind, rings)
            grid::grid.path(grid::unit(xy[, 1], "native"), grid::unit(xy[, 2], "native"), id.lengths = vapply(rings, 
                nrow, 1L), rule = "evenodd", gp = grid::gpar(fill = fill, col = border, lwd = 0.22, linejoin = "round"))
        }
    }
    bar <- function(cols, limits, breaks, labels, title) {
        x <- 0.78
        w <- 4.78
        y <- 0.8
        for (i in seq_along(cols)) box(x + (i - 0.5) * w/length(cols), y, w/length(cols) + 0.001, 0.15, cols[i])
        for (i in seq_along(breaks)) text(labels[i], x + w * (breaks[i] - limits[1])/diff(limits), 0.54, just = "centre")
        text(title, x + w/2, 1.11, just = "centre")
    }
    key <- function(label, col, x, y) {
        box(x + 0.06, y, 0.12, 0.12, col)
        text(label, x + 0.18, y)
    }
    titles <- c("Birthplace distribution", "Recruitment residence distribution", "Recruitment TDI", "Urban and rural areas", 
        "Travel-to-work areas", "Population density")
    ntext <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)
    notes <- c(paste0("N = ", ntext(r$coverage$n[2])), paste0("N = ", ntext(r$coverage$n[3])), "Median among participants in each TTWA", 
        "Country-specific classification", paste0(nrow(r$areas), " TTWAs (2011)"), "2011 residential population; 1 km grid")
    function() {
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(width = grid::unit(21.5, "cm"), height = grid::unit(24.8, "cm")))
        box(10.75, 12.4, 21.5, 24.8, "white")
        for (i in 1:6) {
            col <- (i - 1)%%3
            row <- (i - 1)%/%3
            grid::pushViewport(grid::viewport(x = grid::unit(0.45 + col * 6.93, "cm"), y = grid::unit(24.42 - row * 
                12.05, "cm"), width = grid::unit(6.55, "cm"), height = grid::unit(11.7, "cm"), just = c("left", 
                "top")))
            text(letters[i], 0, 11.5, 10, "bold")
            text(titles[i], 0.42, 11.48, 7)
            text(notes[i], 0.42, 11.09, 6)
            grid::pushViewport(grid::viewport(x = grid::unit(3.275, "cm"), y = grid::unit(6.2, "cm"), width = grid::unit(5.69, 
                "cm"), height = grid::unit(9.35, "cm"), xscale = c(-45000, 710000), yscale = c(-20000, 1220000), 
                clip = "on"))
            if (i <= 2) {
                shapes(area_paths, rep("#F2F2F2", length(area_paths)))
                image <- png::readPNG(file.path(dirname(input), r$layers[i]))
                grid::grid.raster(image, width = 1, height = 1, interpolate = FALSE)
            }
            if (i == 3) 
                shapes(area_paths, tdi_fill)
            if (i == 4) {
                shapes(area_paths, rep(gray, length(area_paths)))
                order <- order(r$urban$class == "Urban")
                shapes(urban_paths[order], ifelse(r$urban$class[order] == "Urban", "#C58D55", "#F0DFC9"))
            }
            if (i == 5) 
                shapes(area_paths, rep(ttwa_fill, length(area_paths)), border = ttwa_line)
            if (i == 6) {
                shapes(area_paths, rep("#F3F3F3", length(area_paths)))
                for (fill in unique(pp$colour)) {
                  d <- pp[colour == fill]
                  s <- r$cellsize/2
                  xx <- as.vector(t(cbind(d$easting - s, d$easting + s, d$easting + s, d$easting - s)))
                  yy <- as.vector(t(cbind(d$northing - s, d$northing - s, d$northing + s, d$northing + s)))
                  grid::grid.path(grid::unit(xx, "native"), grid::unit(yy, "native"), id.lengths = rep(4L, nrow(d)), 
                    gp = grid::gpar(fill = fill, col = NA))
                }
            }
            grid::popViewport()
            if (i <= 2) {
                grid::grid.circle(grid::unit(2.1, "cm"), grid::unit(0.8, "cm"), r = grid::unit(0.015, "cm"), gp = grid::gpar(fill = point_col[i], 
                  col = NA))
                text("One dot per participant", 2.24, 0.8)
            }
            if (i == 3) {
                bar(tdi_col, c(-5, 5), c(-5, 0, 5), c("-5", "0", "5"), "TDI")
                key("No data", "#F2F2F2", 2.65, 0.19)
            }
            if (i == 4) {
                key("Urban", "#C58D55", 1.65, 0.82)
                key("Rural", "#F0DFC9", 3.55, 0.82)
                text("England / Wales 2011; Scotland 2020", 3.275, 0.49, just = "centre")
                text("Northern Ireland 2015", 3.275, 0.2, just = "centre")
            }
            if (i == 5) {
                key("TTWA boundary", ttwa_line, 2.23, 0.82)
                sx <- 2.52
                sw <- 2e+05 * 5.69/755000
                grid::grid.lines(grid::unit(c(sx, sx + sw), "cm"), grid::unit(c(0.36, 0.36), "cm"), gp = grid::gpar(col = "#555555", 
                  lwd = 0.55))
                text("0", sx, 0.12, just = "centre")
                text("200 km", sx + sw, 0.12, just = "centre")
            }
            if (i == 6) {
                ticks <- c(0, 10, 100, 1000, 10000)
                bar(pop_col, pop_lim, log1p(ticks), ntext(ticks), expression("People / km"^2))
                key("No grid data", "#F3F3F3", 2.4, 0.19)
            }
            grid::popViewport()
        }
        grid::popViewport()
    }
}

run_ed2 <- function(r = ed2_data()) {
    draw <- ed2_plot(r)
    message("ED2: vector PowerPoint and source data")
    draw_ppt(file.path(out, "ED2.pptx"), draw, 21.5 * cm, 24.8 * cm)
    note <- function(...) data.table(definition = c(...))
    invisible(NULL)
    print(r$coverage)
    invisible(r)
}

run_ed2(plot_data$result)

finish("ED2")

