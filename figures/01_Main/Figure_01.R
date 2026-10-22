# 1. Figure 1
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

map_curve <- function(x0, y0, x1, y1, n = 41) {
    dx <- x1 - x0
    dy <- y1 - y0
    d <- sqrt(dx^2 + dy^2)
    nx <- -dy/d
    ny <- dx/d
    a <- pmin(70000, pmax(18000, 0.18 * d))
    cx <- (x0 + x1)/2 + nx * a
    cy <- (y0 + y1)/2 + ny * a
    t <- seq(0, 1, length.out = n)
    data.table(t = t, x = (1 - t)^2 * x0 + 2 * (1 - t) * t * cx + t^2 * x1, y = (1 - t)^2 * y0 + 2 * (1 - t) * 
        t * cy + t^2 * y1)
}

fig1_map_size <- function() {
    map_h <- 15 * cm
    list(width = map_h * 760000/1200000, height = map_h, units = "in")
}

fig1_map_plot <- function(r = fig1_map_data()) {
    stopifnot(is.list(r), is.list(r$plotting))
    x <- r$plotting
    map_pal <- list(bg_low = "#2166AC", bg_mid = "#F7F7F7", bg_high = "#B2182B", in_col = "#2AA99A", out_col = "#C18543")
    fs <- 6/ggplot2::.pt
    fs_title <- 7/ggplot2::.pt
    poly <- x$poly
    nodes <- copy(x$nodes)
    edges <- copy(x$edges)
    max_edge <- if (nrow(edges)) 
        max(edges$N)
    else 100
    lw <- function(n) 0.045 + 1.45 * sqrt(n)/sqrt(max_edge)
    la <- function(n) 0.24 + 0.5 * (n/max_edge)^0.7
    if (nrow(edges)) {
        cdf <- rbindlist(Map(function(x0, y0, x1, y1, n, id) {
            z <- map_curve(x0, y0, x1, y1)
            z[, `:=`(id = id, N = n, lw = lw(n), la = la(n))]
        }, edges$x0, edges$y0, edges$x1, edges$y1, edges$N, edges$id))
        cdf <- rbindlist(list(copy(cdf[t <= 0.5])[, `:=`(half = paste0(id, "_out"), col = map_pal$out_col)], copy(cdf[t >= 
            0.5])[, `:=`(half = paste0(id, "_in"), col = map_pal$in_col)]))
    }
    else {
        cdf <- data.table(x = numeric(), y = numeric(), half = character(), col = character(), lw = numeric(), 
            la = numeric())
        warning("No core-to-core flows meet N >= 100; no paths drawn")
    }
    core_nodes <- nodes[is_core == TRUE & total > 0]
    core_nodes[, `:=`(r, if (.N) 
        6900 + 9800 * sqrt(total)/sqrt(max(total))
    else numeric())]
    pie <- function(code, cx, cy, r, inn, outn) {
        vals <- c(inflow = inn, outflow = outn)
        cols <- c(inflow = map_pal$in_col, outflow = map_pal$out_col)
        st <- pi/2
        rbindlist(lapply(names(vals), function(nm) {
            if (vals[[nm]] <= 0) 
                return(NULL)
            en <- st + 2 * pi * vals[[nm]]/sum(vals)
            th <- seq(st, en, length.out = max(8, ceiling(45 * vals[[nm]]/sum(vals))))
            st <<- en
            data.table(piece = paste(code, nm), x = c(cx, cx + r * cos(th), cx), y = c(cy, cy + r * sin(th), cy), 
                col = cols[[nm]])
        }))
    }
    pies <- if (nrow(core_nodes)) 
        rbindlist(Map(pie, core_nodes$code, core_nodes$cx, core_nodes$cy, core_nodes$r, core_nodes$in_n, core_nodes$out_n))
    else data.table(piece = character(), x = numeric(), y = numeric(), col = character())
    other <- nodes[is_core == FALSE & total > 0]
    other[, `:=`(pt = if (.N) 
        0.48 + 2.05 * sqrt(total)/sqrt(max(total))
    else numeric(), col = fifelse(net >= 0, map_pal$in_col, map_pal$out_col))]
    lab <- merge(core_nodes, x$labels, by = "name", all.x = TRUE)
    lab[is.na(dx), `:=`(dx = 9000, dy = 9000)]
    lab[, `:=`(lx = cx + dx, ly = cy + dy)]
    leg_x <- 585000
    flow <- data.table(N = c(100, 500, 2000), x = leg_x, xend = leg_x + 52000, y = c(950000, 915000, 880000))
    flow[, `:=`(lw = lw(N), la = la(N))]
    flow2 <- rbindlist(list(flow[, .(N, x = x, xend = (x + xend)/2, y, lw, la, col = map_pal$out_col)], flow[, 
        .(N, x = (x + xend)/2, xend = xend, y, lw, la, col = map_pal$in_col)]))
    pts <- data.table(N = c(1000, 5000, 20000), x = leg_x + 26000, y = c(742000, 7e+05, 658000))
    pts[, `:=`(pt, 0.48 + 2.05 * sqrt(N)/sqrt(max(1, nodes$total)))]
    tbar <- data.table(v = seq(-5, 5, length.out = 80))
    tbar[, `:=`(xmin = leg_x + 20000, xmax = leg_x + 31000, ymin = 410000 + (v + 5)/10 * 140000, ymax = 410000 + 
        (shift(v, type = "lead", fill = 5) + 5)/10 * 140000)]
    ggplot() + geom_polygon(data = poly, aes(x, y, group = group, subgroup = group, fill = tdi_med), color = NA, 
        alpha = 0.58) + geom_point(data = other, aes(cx, cy, color = col, size = pt), alpha = 1) + geom_path(data = cdf, 
        aes(x, y, group = half, color = col, linewidth = lw, alpha = la), lineend = "butt") + geom_polygon(data = pies, 
        aes(x, y, group = piece), fill = pies$col, color = NA, alpha = 1) + geom_segment(data = lab, aes(cx, cy, 
        xend = lx, yend = ly), color = "#5A5A5A", linewidth = 0.1, alpha = 0.22) + geom_label(data = lab, aes(lx, 
        ly, label = name), size = fs_title, linewidth = 0.1, label.padding = grid::unit(0.055, "lines"), fill = scales::alpha("white", 
        0.76), color = "#252525") + annotate("text", x = leg_x + 26000, y = 1018000, label = "Migration direction", 
        hjust = 0.5, color = "#333333", size = fs_title) + annotate("text", x = leg_x + 13000, y = 985000, label = "Out", 
        hjust = 0.5, color = map_pal$out_col, size = fs) + annotate("text", x = leg_x + 39000, y = 985000, label = "In", 
        hjust = 0.5, color = map_pal$in_col, size = fs) + geom_segment(data = flow2, aes(x, y, xend = xend, yend = y, 
        color = col, linewidth = lw, alpha = la), lineend = "butt") + annotate("text", x = leg_x + 62500, y = flow$y, 
        label = scales::comma(flow$N), hjust = 0, color = "#333333", size = fs) + annotate("text", x = leg_x + 
        26000, y = 790000, label = "Total flow", hjust = 0.5, color = "#333333", size = fs_title) + geom_point(data = pts, 
        aes(x, y, size = pt), color = "#666666", alpha = 1) + annotate("text", x = leg_x + 62500, y = pts$y, label = scales::comma(pts$N), 
        hjust = 0, color = "#333333", size = fs) + annotate("text", x = leg_x + 25500, y = 585000, label = "TDI", 
        hjust = 0.5, color = "#333333", size = fs_title) + geom_rect(data = tbar, aes(xmin = xmin, xmax = xmax, 
        ymin = ymin, ymax = ymax, fill = v), color = NA) + annotate("text", x = leg_x + 39000, y = c(410000, 480000, 
        550000), label = c("-5", "0", "5"), hjust = 0, color = "#333333", size = fs) + scale_fill_gradient2(low = map_pal$bg_low, 
        mid = map_pal$bg_mid, high = map_pal$bg_high, midpoint = 0, limits = c(-5, 5), oob = scales::squish, na.value = "#F2F2F2") + 
        scale_color_identity() + scale_size_identity() + scale_linewidth_identity() + scale_alpha_identity() + 
        coord_fixed(xlim = c(-45000, 715000), ylim = c(15000, 1215000), expand = FALSE, clip = "off") + guides(fill = "none", 
        color = "none", size = "none", linewidth = "none", alpha = "none") + theme_void(base_size = 6) + theme(plot.margin = margin(0, 
        0, 0, 0), panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", 
        color = NA))
}

fig1_description <- function(r) {
    ink <- "#252525"
    teal <- "#168D83"
    red <- "#B2182B"
    blue <- "#2667AD"
    purple <- "#77519D"
    urban <- "#D6A474"
    rural <- "#B59B7D"
    density_col <- pal$density
    label <- function(s, x, y, size = 6, col = ink, face = "plain", just = "centre", rot = 0) grid::grid.text(s, 
        grid::unit(x, "cm"), grid::unit(y, "cm"), just = just, rot = rot, gp = grid::gpar(fontfamily = font, fontsize = size, 
            col = col, fontface = face))
    line <- function(x, y, col = ink, width = 0.35, tip = FALSE, dash = "solid") grid::grid.lines(grid::unit(x, 
        "cm"), grid::unit(y, "cm"), gp = grid::gpar(col = col, lwd = width, lineend = "round", lty = dash), arrow = if (tip) 
        grid::arrow(length = grid::unit(0.055, "cm"), type = "closed")
    else NULL)
    rect <- function(x, y, w, h, fill, col = NA, width = 0.35) grid::grid.rect(grid::unit(x, "cm"), grid::unit(y, 
        "cm"), grid::unit(w, "cm"), grid::unit(h, "cm"), gp = grid::gpar(fill = fill, col = col, lwd = width))
    poly <- function(x, y, fill, col = NA, width = 0.35) grid::grid.polygon(grid::unit(x, "cm"), grid::unit(y, 
        "cm"), gp = grid::gpar(fill = fill, col = col, lwd = width))
    count <- function(variable, route) r$composition[Variable == variable & Route == route, N]
    pop <- if (!is.null(r)) 
        r$density[section == "Density change"]
    else NULL
    pop_density <- if (!is.null(pop)) 
        list(x = pop$x, y = pop$density)
    else NULL
    descriptive <- function() {
        H <- 15
        txt <- function(s, x, y, ...) label(s, x, H - y, ...)
        seg <- function(x, y, ...) line(x, H - y, ...)
        box <- function(x, y, w, h, ...) rect(x, H - y, w, h, ...)
        area <- function(x, y, fill, ...) poly(x, H - y, fill, ...)
        fmt <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)
        density <- function(d, x, y, w, h, col, ticks, ticklabels = ticks, unit = "", zero = FALSE) {
            rx <- range(d$x)
            xx <- function(v) x + (v - rx[1])/diff(rx) * w
            yt <- pretty(c(0, max(d$y)), n = 2)
            ymax <- max(yt)
            yy <- function(v) y + h - v/ymax * h
            area(c(xx(d$x), rev(xx(d$x))), c(yy(d$y), rep(y + h, length(d$x))), adjustcolor(col, alpha.f = 0.21))
            seg(xx(d$x), yy(d$y), col, 0.7)
            if (zero) 
                seg(rep(xx(0), 2), c(y + h, yy(approx(d$x, d$y, xout = 0)$y)), col, 0.35, dash = "dotted")
            seg(c(x, x + w), rep(y + h, 2), "#717979", 0.35)
            seg(rep(x, 2), c(y, y + h), "#717979", 0.35)
            for (v in yt) {
                seg(x + c(-0.055, 0), rep(yy(v), 2), "#717979", 0.35)
                txt(format(v, trim = TRUE), x - 0.12, yy(v), 6, just = "right")
            }
            txt("Density", 0.13, y + h/2, 6, rot = 90)
            for (i in seq_along(ticks)) {
                seg(rep(xx(ticks[i]), 2), y + h + c(0, 0.055), "#717979", 0.35)
                txt(as.character(ticklabels[i]), xx(ticks[i]), y + h + 0.18, 6)
            }
            txt(unit, x + w/2, y + h + 0.47, 6)
        }
        sankey <- function(variable, routes, groups, cols, y, height, flowcols) {
            n <- vapply(routes, function(z) count(variable, z), numeric(1))
            p <- n/sum(n)
            gap <- 0.22
            scale <- height - gap
            left <- c(y, y + (p[1] + p[2]) * scale + gap)
            right <- c(y, y + (p[1] + p[3]) * scale + gap)
            starts <- c(left[1], left[1] + p[1] * scale, left[2], left[2] + p[3] * scale)
            ends <- c(right[1], right[2], right[1] + p[1] * scale, right[2] + p[2] * scale)
            x1 <- 0.77
            x2 <- 4.55
            u <- seq(0, 1, length.out = 120)
            ease <- 3 * u^2 - 2 * u^3
            for (i in seq_along(n)) {
                top <- starts[i] + (ends[i] - starts[i]) * ease
                area(c(x1 + (x2 - x1) * u, rev(x1 + (x2 - x1) * u)), c(top, rev(top + p[i] * scale)), adjustcolor(flowcols[i], 
                  alpha.f = 0.36))
            }
            for (k in 1:2) {
                hs <- c(sum(p[1:2]), sum(p[3:4]))[k] * scale
                hd <- c(sum(p[c(1, 3)]), sum(p[c(2, 4)]))[k] * scale
                box(x1, left[k] + hs/2, 0.095, hs, cols[k])
                box(x2, right[k] + hd/2, 0.095, hd, cols[k])
                txt(groups[k], x1 - 0.13, left[k] + hs/2, 6, cols[k], just = "right")
                txt(groups[k], x2 + 0.13, right[k] + hd/2, 6, cols[k], just = "left")
            }
            txt("Birthplace", 1.06, y - 0.24, 6)
            txt("Recruitment", 4.2, y - 0.24, 6)
            pos <- c(0.5, 0.22, 0.79, 0.5)
            for (i in seq_along(n)) {
                e <- 3 * pos[i]^2 - 2 * pos[i]^3
                cy <- starts[i] + (ends[i] - starts[i]) * e + p[i] * scale/2
                cx <- x1 + (x2 - x1) * pos[i]
                if (p[i] * scale < 0.12) {
                  seg(rep(cx, 2), c(cy, y + height + 0.13), flowcols[i], 0.35)
                  cy <- y + height + 0.25
                }
                txt(fmt(n[i]), cx, cy, 6, ink)
            }
        }
        txt("MD", 2.7, 0.19, 7, teal)
        d <- r$density[section == "MD RINT"]
        density(list(x = d$x, y = d$density), 0.84, 0.6, 4.15, 1.58, teal, c(-4, 0, 4), unit = "RINT")
        txt("ΔTDI", 2.7, 3.02, 7, blue)
        sankey("tdi4", c("High to high", "High to low", "Low to high", "Low to low"), c("High", "Low"), c(red, 
            blue), 3.63, 1.91, c(red, "#6497C4", "#D78691", blue))
        txt("Rural–urban", 2.7, 6.14, 7, "#A58A6D")
        sankey("urban_route", c("Urban to Urban", "Urban to Rural", "Rural to Urban", "Rural to Rural"), c("Urban", 
            "Rural"), c(urban, rural), 6.75, 1.91, c(urban, rural, "#E2B993", "#9E876D"))
        txt("Cross TTWA", 2.7, 9.35, 7, purple)
        nn <- c(count("ttwa", "Within TTWA"), count("ttwa", "Cross TTWA"))
        x0 <- 1.1
        w <- 3.04
        xmax <- 225000
        for (i in 1:2) {
            y <- 9.91 + (i - 1) * 0.52
            bw <- nn[i]/xmax * w
            box(x0 + bw/2, y, bw, 0.28, c("#C8B8D9", purple)[i])
            txt(c("Within", "Cross")[i], x0 - 0.13, y, 6, just = "right")
            txt(fmt(nn[i]), x0 + bw + 0.1, y, 6, just = "left")
        }
        seg(c(x0, x0 + w), rep(10.77, 2), "#717979", 0.35)
        for (v in c(0, 1e+05, 2e+05)) {
            xx <- x0 + v/xmax * w
            seg(rep(xx, 2), 10.77 + c(0, 0.055), "#717979", 0.35)
            txt(fmt(v), xx, 10.97, 6)
        }
        txt("Participants", x0 + w/2, 11.27, 6)
        txt("Density change", 2.7, 11.92, 7, density_col)
        density(pop_density, 0.84, 12.31, 4.15, 1.81, density_col, c(-8, 0, 8), unit = "Log density change", zero = TRUE)
    }
    descriptive
}

fig1_assemble <- function(file, map_plot, description_plot, map_file = NULL) {
    tmp <- tempfile("fig1_")
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    ns <- c(p = "http://schemas.openxmlformats.org/presentationml/2006/main", a = "http://schemas.openxmlformats.org/drawingml/2006/main", 
        r = "http://schemas.openxmlformats.org/officeDocument/2006/relationships")
    sizes <- list(concept = c(14.39, 4.3), map = c(8.64, 15), description = c(5.4, 15))
    stopifnot(file.copy(file.path(dirname(input), "Fig1a.pptx"), file.path(tmp, "concept.pptx")))
    plots <- list(map = map_plot, description = description_plot)
    for (nm in names(plots)) {
        dest <- file.path(tmp, paste0(nm, ".pptx"))
        if (is.function(plots[[nm]])) 
            draw_ppt(dest, plots[[nm]], sizes[[nm]][1] * cm, sizes[[nm]][2] * cm)
        else write_ppt(dest, list(plots[[nm]]), sizes[[nm]][1] * cm, sizes[[nm]][2] * cm)
    }
    art <- lapply(names(sizes), function(nm) {
        pd <- file.path(tmp, nm)
        utils::unzip(file.path(tmp, paste0(nm, ".pptx")), exdir = pd)
        doc <- xml2::read_xml(file.path(pd, "ppt/slides/slide1.xml"))
        if (length(xml2::xml_find_all(doc, "//p:pic | //a:blip | //*[@r:embed or @r:link]", ns))) 
            stop("Fig1 panel must contain only editable native shapes: ", nm)
        for (f in xml2::xml_find_all(doc, "//a:latin | //a:ea | //a:cs", ns)) xml2::xml_set_attr(f, "typeface", 
            font)
        as.list(xml2::xml_children(xml2::xml_find_first(doc, "//p:spTree", ns))[-(1:2)])
    })
    names(art) <- names(sizes)
    base <- file.path(tmp, "concept")
    slide_file <- file.path(base, "ppt/slides/slide1.xml")
    doc <- xml2::read_xml(slide_file)
    tree <- xml2::xml_find_first(doc, "//p:spTree", ns)
    xml2::xml_remove(xml2::xml_children(tree)[-(1:2)])
    emu <- function(x) sprintf("%.0f", x * 360000)
    place <- function(name, x, y, w, h, rotate = FALSE) {
        wh <- sizes[[name]]
        expected <- if (rotate) 
            rev(wh)
        else wh
        stopifnot(isTRUE(all.equal(expected, c(w, h))))
        xy <- c(x + w/2, y + h/2) - wh/2
        g <- xml2::read_xml(sprintf(paste0("<p:grpSp xmlns:p=\"%s\" xmlns:a=\"%s\"><p:nvGrpSpPr>", "<p:cNvPr id=\"1\" name=\"fig1.%s\"/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr>", 
            "<p:grpSpPr><a:xfrm rot=\"%s\"><a:off x=\"%s\" y=\"%s\"/><a:ext cx=\"%s\" cy=\"%s\"/>", "<a:chOff x=\"0\" y=\"0\"/><a:chExt cx=\"%s\" cy=\"%s\"/></a:xfrm></p:grpSpPr></p:grpSp>"), 
            ns["p"], ns["a"], name, if (rotate) 
                "16200000"
            else "0", emu(xy[1]), emu(xy[2]), emu(wh[1]), emu(wh[2]), emu(wh[1]), emu(wh[2])))
        for (child in art[[name]]) xml2::xml_add_child(g, child)
        xml2::xml_add_child(tree, g)
    }
    place("concept", 0.35, 0.45, 14.39, 4.3)
    place("description", 0.35, 5.1, 5.4, 15)
    place("map", 6.1, 5.1, 8.64, 15)
    if (!is.null(map_file)) {
        old <- xml2::read_xml(unzip(map_file, files = "ppt/slides/slide1.xml", exdir = file.path(tmp, "map_saved")))
        path <- "//p:grpSp[p:nvGrpSpPr/p:cNvPr[@name=\"fig1.map\"]]"
        node <- xml2::xml_find_first(old, path, ns)
        stopifnot(!inherits(node, "xml_missing"))
        off <- xml2::xml_find_first(node, "./p:grpSpPr/a:xfrm/a:off", ns)
        xml2::xml_set_attr(off, "x", emu(6.1))
        xml2::xml_replace(xml2::xml_find_first(doc, path, ns), node)
    }
    tag <- function(label, x, y) {
        shape <- xml2::read_xml(sprintf(paste0("<p:sp xmlns:p=\"%s\" xmlns:a=\"%s\"><p:nvSpPr>", "<p:cNvPr id=\"1\" name=\"panel %s\"/><p:cNvSpPr txBox=\"1\"/><p:nvPr/></p:nvSpPr>", 
            "<p:spPr><a:xfrm><a:off x=\"%s\" y=\"%s\"/><a:ext cx=\"144000\" cy=\"144000\"/>", "</a:xfrm><a:prstGeom prst=\"rect\"><a:avLst/></a:prstGeom><a:noFill/></p:spPr>", 
            "<p:txBody><a:bodyPr lIns=\"0\" rIns=\"0\" tIns=\"0\" bIns=\"0\"/><a:lstStyle/><a:p><a:r>", "<a:rPr sz=\"1000\" b=\"1\"><a:solidFill><a:srgbClr val=\"000000\"/></a:solidFill>", 
            "<a:latin typeface=\"%s\"/></a:rPr><a:t>%s</a:t></a:r></a:p></p:txBody></p:sp>"), ns["p"], ns["a"], 
            label, emu(x), emu(y), font, label))
        xml2::xml_add_child(tree, shape)
    }
    tag("a", 0.1, 0.15)
    tag("b", 0.1, 4.85)
    tag("c", 6.1, 4.85)
    ids <- xml2::xml_find_all(doc, "//p:cNvPr", ns)
    for (i in seq_along(ids)) xml2::xml_set_attr(ids[[i]], "id", as.character(i))
    xml2::write_xml(doc, slide_file)
    pres_file <- file.path(base, "ppt/presentation.xml")
    pres <- xml2::read_xml(pres_file)
    xml2::xml_set_attrs(xml2::xml_find_first(pres, "//p:sldSz", ns), c(cx = emu(15.04), cy = emu(20.25)))
    xml2::write_xml(pres, pres_file)
    target <- file.path(tmp, "Fig1.pptx")
    zip::zipr(target, list.files(base, recursive = TRUE, all.files = TRUE, no.. = TRUE), root = base, mode = 'mirror')
    if (!file.copy(target, file, overwrite = TRUE)) 
        stop("Cannot replace Fig1 presentation: ", file)
    invisible(file)
}

fig1_assemble(file.path(out, "Fig1.pptx"), fig1_map_plot(plot_data$map), fig1_description(plot_data$description))

finish("Fig1")

