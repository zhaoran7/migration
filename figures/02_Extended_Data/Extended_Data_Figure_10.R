# 1. Extended Data Figure 10
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

enrich_cloud <- function(d, colour, size = 0.33) {
    xr <- range(d$x)
    yr <- range(d$y)
    xr <- xr + c(-1, 1) * diff(xr) * 0.01
    yr <- yr + c(-1, 1) * diff(yr) * 0.01
    x <- (d$x - xr[1])/diff(xr)
    y <- (d$y - yr[1])/diff(yr)
    radius <- if (size > 0.2) 
        0.0026
    else 0.0012
    angle <- seq(0, 2 * pi, length.out = 13)[-13]
    rx <- radius * cos(angle)
    ry <- radius * diff(xr)/diff(yr) * sin(angle)
    pieces <- unlist(lapply(unique(colour), function(col) {
        i <- which(colour == col)
        lapply(split(i, ceiling(seq_along(i)/300)), function(k) grid::pathGrob(as.vector(t(outer(x[k], rx, "+"))), 
            as.vector(t(outer(y[k], ry, "+"))), id = rep(seq_along(k), each = 12), rule = "winding", gp = grid::gpar(fill = col, 
                col = NA)))
    }), recursive = FALSE)
    grob <- do.call(grid::grobTree, unname(pieces))
    p <- ggplot() + annotation_custom(grob, xmin = xr[1], xmax = xr[2], ymin = yr[1], ymax = yr[2]) + scale_x_continuous(limits = xr, 
        expand = c(0, 0)) + scale_y_continuous(limits = yr, expand = c(0, 0)) + coord_fixed() + theme_void(base_family = font, 
        base_size = 6) + theme(plot.margin = margin(0, 0, 0, 0))
    p
}

enrich_umap <- function(r, atlas) {
    d <- r$umap[Atlas == atlas]
    keys <- sort(unique(d$Cell))
    cols <- setNames(grDevices::hcl.colors(length(keys), "Dark 3"), keys)
    selected <- r$cells[Atlas == atlas & Display == TRUE, unique(Cell)]
    lab <- r$counts[Atlas == atlas]
    lab[, `:=`(priority, Cell %in% selected)]
    lab <- lab[order(-priority, -Cells), head(.SD, 8)]
    short <- c(`near-projecting glutamatergic cortical neuron` = "Near-projecting neuron", `intratelencephalic-projecting glutamatergic cortical neuron` = "Intratelencephalic neuron", 
        `caudal ganglionic eminence derived interneuron` = "CGE interneuron", `medial ganglionic eminence derived interneuron` = "MGE interneuron", 
        `cerebral cortex neuron` = "Cortical neuron", `oligodendrocyte precursor cell` = "OPC", `capillary endothelial cell` = "Capillary endothelium", 
        `endothelial cell of artery` = "Arterial endothelium", `hematopoietic precursor cell` = "Hematopoietic precursor", 
        `CD4-positive alpha-beta T cell` = "CD4 T cell", `CD4-positive, alpha-beta T cell` = "CD4 T cell")
    lab[, `:=`(name, gsub("_", " ", Cell))]
    lab[, `:=`(label, fifelse(name %in% names(short), unname(short[name]), name))]
    lab[, `:=`(label, vapply(label, function(s) paste(strwrap(s, width = 26), collapse = "\n"), ""))]
    lab[, `:=`(dx = 0, dy = 0)]
    if (atlas == "Immune") {
        lab[Cell == "Cycling T&NK", `:=`(dx = 2.4, dy = 1.6)]
        lab[Cell == "Erythroid", `:=`(dx = -0.7, dy = -1.2)]
    }
    enrich_cloud(d, unname(cols[d$Cell]), 0.1) + geom_text_repel(data = lab, aes(x, y, label = label), inherit.aes = FALSE, 
        size = pt(6), family = font, color = "#222222", seed = 13, box.padding = 0.25, point.padding = 0.06, min.segment.length = 0, 
        segment.size = 0.1, nudge_x = lab$dx, nudge_y = lab$dy, max.overlaps = Inf, force = 4, force_pull = 0.05, 
        max.iter = 1e+05, max.time = Inf) + labs(title = atlas) + theme(plot.title = element_text(size = 7, hjust = 0.5, 
        margin = margin(b = 3)))
}

enrich_heat <- function(d, group, x, y, h, tile = 0.27, pitch = tile, category = NULL, limits = c(0, 15), colours = enrich_pal) {
    rows <- d[, .(best = min(p)), by = c(if (!is.null(category)) "Atlas", group)]
    if (is.null(category)) 
        rows <- rows[order(best)]
    else rows <- rows[order(match(Atlas, names(category)), best)]
    rows[, `:=`(yy, y + (seq_len(.N) - 0.5) * pitch)]
    text <- function(s, xx, yy, size = 6, just = "left", angle = 0) grid::grid.text(s, unit(xx, "cm"), unit(h - 
        yy, "cm"), just = just, rot = angle, gp = grid::gpar(fontfamily = font, fontsize = size, col = "#222222"))
    if (!is.null(category)) {
        blocks <- rows[, .(top = min(yy) - pitch/2, bottom = max(yy) + pitch/2), by = Atlas]
        for (i in seq_len(nrow(blocks))) {
            b <- blocks[i]
            grid::grid.rect(unit(x - 0.5, "cm"), unit(h - (b$top + b$bottom)/2, "cm"), unit(0.7, "cm"), unit(b$bottom - 
                b$top, "cm"), gp = grid::gpar(fill = category[b$Atlas], col = "#C8C8C8", lwd = 0.35))
            text(b$Atlas, x - 0.5, (b$top + b$bottom)/2, 6, "centre", 90)
        }
    }
    for (j in 1:2) text(c("MD", dtdi_lab)[j], x + (j - 0.5) * tile, y - 0.35, 6, "centre", 90)
    ramp <- colorRampPalette(colours)(256)
    for (i in seq_len(nrow(rows))) {
        yy <- rows$yy[i]
        for (j in 1:2) {
            take <- d[[group]] == rows[[group]][i] & d$Trait == c("MD", dtdi_lab)[j]
            if (!is.null(category)) 
                take <- take & d$Atlas == rows$Atlas[i]
            a <- d[which(take)]
            stopifnot(nrow(a) == 1)
            score <- pmax(limits[1], pmin(a$score, limits[2]))
            fill <- if (is.finite(score)) 
                ramp[1 + round((score - limits[1])/diff(limits) * 255)]
            else "white"
            grid::grid.rect(unit(x + (j - 0.5) * tile, "cm"), unit(h - yy, "cm"), unit(tile, "cm"), unit(tile, 
                "cm"), gp = grid::gpar(fill = fill, col = "#C8C8C8", lwd = 0.35))
            if (a$sig) 
                text("*", x + (j - 0.5) * tile, yy, 7, "centre")
        }
        label <- gsub("_", " ", rows[[group]][i])
        if (group != "Tissue") 
            label <- paste(strwrap(label, width = 40), collapse = "\n")
        text(label, x + 2 * tile + 0.13, yy, 6)
    }
    invisible(y + nrow(rows) * pitch)
}

enrich_bar <- function(x, y, w, cap, h, limits = c(0, cap), colours = enrich_pal, title = expression(paste("-", 
    log[10], "(", italic(p), ")"))) {
    colour <- colorRampPalette(colours)(256)
    grid::grid.rect(unit(x + ((1:256) - 0.5) * w/256, "cm"), unit(h - y, "cm"), width = unit(w/256 * 1.015, "cm"), 
        height = unit(0.16, "cm"), gp = grid::gpar(fill = colour, col = NA))
    ticks <- sort(unique(c(limits, pretty(limits, n = 3))))
    ticks <- ticks[ticks >= limits[1] & ticks <= limits[2]]
    grid::grid.text(format(ticks, trim = TRUE), unit(x + (ticks - limits[1])/diff(limits) * w, "cm"), unit(h - 
        y - 0.24, "cm"), gp = grid::gpar(fontfamily = font, fontsize = 6))
    grid::grid.text(title, unit(x + w/2, "cm"), unit(h - y + 0.25, "cm"), gp = grid::gpar(fontfamily = font, fontsize = 7))
}

enrich_ppt <- function(file, plots, positions, draw, w, h) {
    vector <- identical(Sys.getenv("MOVE_VECTOR", "0"), "1")
    stopifnot(requireNamespace("ragg", quietly = TRUE))
    if (vector) 
        stopifnot(requireNamespace("svglite", quietly = TRUE))
    stage <- tempfile(fileext = ".pptx")
    svg <- vapply(plots, function(p) tempfile(fileext = ".svg"), "")
    png <- vapply(plots, function(p) tempfile(fileext = ".png"), "")
    package <- tempfile()
    dir.create(package)
    on.exit({
        unlink(c(stage, svg, png))
        unlink(package, recursive = TRUE)
    }, add = TRUE)
    draw_ppt(stage, function() draw(FALSE), w * cm, h * cm)
    deck <- read_pptx(stage)
    for (i in seq_along(plots)) {
        plate <- plots[[i]]
        plate$layers <- plate$layers[1]
        plate <- plate + theme(plot.title = element_text(color = "transparent"))
        if (vector) {
            svglite::svglite(svg[i], width = positions$w[i] * cm, height = positions$h[i] * cm)
            print(plate)
            dev.off()
        }
        ragg::agg_png(png[i], width = positions$w[i], height = positions$h[i], units = "cm", res = 600)
        print(plate)
        dev.off()
        deck <- ph_with(deck, external_img(png[i], alt = "Enrichment point cloud"), location = ph_location(left = positions$x[i] * 
            cm, top = positions$y[i] * cm, width = positions$w[i] * cm, height = positions$h[i] * cm))
    }
    for (i in seq_along(plots)) {
        label <- plots[[i]]
        label$layers <- label$layers[-1L]
        label <- label + theme(plot.background = element_rect(fill = "transparent", color = NA), panel.background = element_rect(fill = "transparent", 
            color = NA))
        deck <- ph_with(deck, rvg::dml(ggobj = label, bg = "transparent"), location = ph_location(left = positions$x[i] * 
            cm, top = positions$y[i] * cm, width = positions$w[i] * cm, height = positions$h[i] * cm))
    }
    print(deck, target = stage)
    if (vector) {
        unzip(stage, exdir = package)
        slide <- file.path(package, "ppt/slides/slide1.xml")
        relfile <- file.path(package, "ppt/slides/_rels/slide1.xml.rels")
        doc <- xml2::read_xml(slide)
        rel <- xml2::read_xml(relfile)
        ns <- c(p = "http://schemas.openxmlformats.org/presentationml/2006/main", a = "http://schemas.openxmlformats.org/drawingml/2006/main")
        blips <- xml2::xml_find_all(doc, ".//p:pic/p:blipFill/a:blip", ns)
        stopifnot(length(blips) == length(svg))
        for (i in seq_along(svg)) {
            name <- paste0("enrich", i, ".svg")
            id <- paste0("rIdVector", i)
            file.copy(svg[i], file.path(package, "ppt/media", name))
            xml2::xml_add_child(rel, xml2::read_xml(paste0("<Relationship xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\"", 
                " Id=\"", id, "\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/image\"", 
                " Target=\"../media/", name, "\"/>")))
            xml2::xml_add_child(blips[[i]], xml2::read_xml(paste0("<a:extLst xmlns:a=\"", ns["a"], "\"><a:ext uri=\"{96DAC541-7B7A-43D3-8B79-37D633B846F1}\">", 
                "<asvg:svgBlip xmlns:asvg=\"http://schemas.microsoft.com/office/drawing/2016/SVG/main\"", " xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships\"", 
                " r:embed=\"", id, "\"/></a:ext></a:extLst>")))
        }
        types <- file.path(package, "[Content_Types].xml")
        content <- xml2::read_xml(types)
        if (!length(xml2::xml_find_all(content, ".//*[local-name()='Default'][@Extension='svg']"))) 
            xml2::xml_add_child(content, xml2::read_xml("<Default xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\" Extension=\"svg\" ContentType=\"image/svg+xml\"/>"))
        xml2::write_xml(doc, slide)
        xml2::write_xml(rel, relfile)
        xml2::write_xml(content, types)
        zip::zipr(stage, list.files(package, recursive = TRUE, all.files = TRUE, no.. = TRUE), root = package, 
            mode = "mirror")
    }
    if (!file.copy(stage, file, overwrite = TRUE)) 
        stop("Close the presentation before replacing: ", file)
    invisible(file)
}

run_ed10 <- function() {
    r <- plot_data$result
    zlim <- range(r$cells[Display == TRUE]$Z, 0, na.rm = TRUE)
    zlim <- c(floor(zlim[1]), ceiling(zlim[2]))
    zcol <- c("#F4FAFB", med_flow[["MD"]], med_col[["MD"]])
    maps <- lapply(c("MD", dtdi_lab), function(k) {
        d <- r$spatial[Trait == k]
        col <- colorRampPalette(enrich_pal)(256)[1 + round(pmin(d$score, 15)/15 * 255)]
        enrich_cloud(d, col) + labs(title = k) + theme(plot.title = element_text(size = 7, hjust = 0.5, margin = margin(b = 3)))
    })
    umap <- lapply(r$atlas$Atlas, function(k) enrich_umap(r, k))
    colours <- setNames(c("#CDE8DF", "#DCD5ED", "#F4E5BA", "#D2E5F1"), r$atlas$Atlas)
    positions <- list(a = data.table(x = c(3.8, 11.55), y = 0.6, w = 7.4, h = 7.9), b = data.table(x = rep(c(0.35, 
        6.75), 2), y = rep(c(0.75, 6.5), each = 2), w = 6, h = 5.3))
    plots <- list(a = maps, b = umap)
    sizes <- list(a = c(19.5, 9.45), b = c(19.5, 12.1))
    draws <- lapply(c("a", "b"), function(id) {
        force(id)
        h <- sizes[[id]][2]
        pos <- positions[[id]]
        function(points = TRUE) {
            grid::grid.text(id, unit(0.2, "cm"), unit(h - 0.18, "cm"), just = c("left", "top"), gp = grid::gpar(fontfamily = font, 
                fontsize = 10, fontface = "bold"))
            if (points) 
                for (i in seq_along(plots[[id]])) print(plots[[id]][[i]], vp = grid::viewport(x = unit(pos$x[i], 
                  "cm"), y = unit(h - pos$y[i], "cm"), width = unit(pos$w[i], "cm"), height = unit(pos$h[i], "cm"), 
                  just = c("left", "top")))
            if (id == "a") {
                enrich_heat(r$tissue, "Tissue", 0.5, 0.85, h, tile = 0.28)
                enrich_bar(6.3, 9, 9, 15, h)
            }
            else {
                bottom <- enrich_heat(r$cells[Display == TRUE], "Cell", 13.85, 1.35, h, tile = 0.34, pitch = 0.46, 
                  category = colours, limits = zlim, colours = zcol)
                enrich_bar(13.85, bottom + 0.7, 4.65, max(zlim), h, limits = zlim, colours = zcol, title = "Z score")
            }
        }
    })
    names(draws) <- names(plots)
    offset <- c(a = 0, b = sizes$a[2] + 0.3)
    size <- c(19.5, sum(vapply(sizes, `[`, numeric(1), 2)) + 0.3)
    draw <- function(points = TRUE) for (id in names(draws)) {
        grid::pushViewport(grid::viewport(x = unit(0, "cm"), y = unit(size[2] - offset[id], "cm"), width = unit(size[1], 
            "cm"), height = unit(sizes[[id]][2], "cm"), just = c("left", "top")))
        draws[[id]](points)
        grid::upViewport()
    }
    pos <- rbindlist(lapply(names(positions), function(id) {
        d <- copy(positions[[id]])
        d[, `:=`(y, y + offset[id])]
        d
    }))
    enrich_ppt(file.path(out, "ED10.pptx"), c(maps, umap), pos, draw, size[1], size[2])
    invisible(NULL)
    unlink(file.path(out, paste0("ED10", c("a", "b"), ".pptx")))
    invisible(list(data = r, draw = draw, size = size))
}

run_ed10()

finish("ED10")

