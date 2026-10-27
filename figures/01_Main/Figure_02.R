# 1. Figure 02
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

fig2b_plot <- function(r) {
    cols <- c(md = "#5AA245", dtdi = "#B2182B", ttwa = "#7A4EA3", density = pal$density)
    labs0 <- c(md = "MD", dtdi = dtdi_lab, ttwa = "Migration\ncross TTWA", density = "Density change")
    a <- rbindlist(list(r$main, r$density), fill = TRUE)
    a[, `:=`(filled, p.adjust(p, "BH") < 0.05)]
    a[, `:=`(exposure, factor(exposure, levels = rev(names(labs0)), labels = rev(labs0)))]
    acol <- setNames(cols, labs0)
    p1 <- ggplot(a, aes(estimate, exposure, color = exposure)) + geom_vline(xintercept = 0, color = "grey65", linewidth = 0.14) + 
        geom_segment(aes(x = lower, xend = upper, yend = exposure), linewidth = 0.28) + geom_point(aes(shape = filled), 
        size = 1.65, stroke = 0.45) + scale_color_manual(values = acol, guide = "none") + scale_shape_manual(values = c(`FALSE` = 1, 
        `TRUE` = 16), guide = "none") + scale_x_continuous(expression(italic(beta) ~ "(95% CI)"), breaks = breaks_pretty(4), 
        expand = expansion(mult = c(0.08, 0.08))) + labs(x = expression(italic(beta) ~ "(95% CI)"), y = NULL, title = "SES") + 
        theme0() + theme(axis.ticks.y = element_blank(), aspect.ratio = 1)
    p2 <- ggplot(r$md, aes(x, fit)) + geom_hline(yintercept = 0, color = "grey75", linewidth = 0.12) + geom_ribbon(aes(ymin = lo, 
        ymax = hi), fill = "#ACDFD2", alpha = 0.45, linewidth = 0) + geom_line(color = "#168D98", linewidth = 0.55) + 
        labs(x = "MD (km)", y = "SES (SD)", title = "MD") + theme0() + theme(aspect.ratio = 1)
    td <- copy(r$tdi)
    m <- max(abs(td$fit), na.rm = TRUE)
    td[, `:=`(light, abs(fit) > 0.62 * m)]
    p3 <- ggplot(td, aes(factor(destination, c("Low", "High")), factor(origin, c("Low", "High")), fill = fit)) + 
        geom_tile(color = "white", linewidth = 0.35) + geom_text(aes(label = sprintf("%.2f", fit), color = light), 
        family = font, size = pt(6)) + scale_fill_gradient(low = "#2667AD", high = "#F7FBFF", limits = c(-m, 0), 
        breaks = c(-m, -m/2, 0), labels = sprintf("%.1f", c(-m, -m/2, 0)), name = "SES (SD)") + scale_color_manual(values = c(`FALSE` = "#222222", 
        `TRUE` = "white"), guide = "none") + labs(x = "Recruitment TDI", y = "Birthplace TDI", title = "TDI groups") + 
        coord_equal() + theme0() + guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(1.75, 
        "cm"), barheight = unit(0.12, "cm"))) + theme(legend.position = "bottom", legend.key.width = unit(1.75, 
        "cm"), legend.key.height = unit(0.12, "cm"), legend.margin = margin(0, 0, 0, 0), aspect.ratio = 1)
    cp <- copy(r$urban)
    m <- max(abs(cp$fit), na.rm = TRUE)
    cp[, `:=`(light, abs(fit) > 0.62 * m)]
    p4 <- ggplot(cp, aes(factor(destination, c("Urban", "Rural")), factor(origin, c("Urban", "Rural")), fill = fit)) + 
        geom_tile(color = "white", linewidth = 0.35) + geom_text(aes(label = sprintf("%.2f", fit), color = light), 
        family = font, size = pt(6)) + scale_fill_gradient(low = "#A66A3F", high = "#FFF5E8", limits = c(-m, 0), 
        breaks = c(-m, -m/2, 0), labels = sprintf("%.1f", c(-m, -m/2, 0)), name = "SES (SD)") + scale_color_manual(values = c(`FALSE` = "#222222", 
        `TRUE` = "white"), guide = "none") + guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, 
        barwidth = unit(1.75, "cm"), barheight = unit(0.12, "cm"))) + labs(x = "Recruitment", y = "Birthplace", 
        title = "Rural–urban routes") + coord_equal() + theme0() + theme(axis.ticks = element_blank(), legend.position = "bottom", 
        legend.key.width = unit(1.75, "cm"), legend.key.height = unit(0.12, "cm"), legend.margin = margin(0, 0, 
            0, 0), aspect.ratio = 1)
    (p3 | p4)/(p1 | p2) + plot_layout(widths = c(1, 1), heights = c(1, 1)) & theme(plot.margin = margin(2, 2, 2, 
        2), plot.tag = element_text(family = font, size = 7, face = "bold"), plot.tag.position = c(0.01, 0.99))
}

fig2c_keys <- function(file) {
    tmp <- tempfile("fig2c_")
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    unzip(file, exdir = tmp)
    path <- file.path(tmp, "ppt/slides/slide1.xml")
    doc <- xml2::read_xml(path)
    ns <- c(p = "http://schemas.openxmlformats.org/presentationml/2006/main", a = "http://schemas.openxmlformats.org/drawingml/2006/main")
    if (length(xml2::xml_find_all(doc, ".//p:cNvPr[@name='route_keys']", ns))) 
        return(invisible(file))
    sp <- xml2::xml_find_all(doc, ".//p:sp", ns)
    tx <- vapply(sp, function(z) paste(xml2::xml_text(xml2::xml_find_all(z, ".//a:t", ns)), collapse = ""), "")
    off <- xml2::xml_find_first(sp, "./p:spPr/a:xfrm/a:off", ns)
    ext <- xml2::xml_find_first(sp, "./p:spPr/a:xfrm/a:ext", ns)
    xy <- data.table(x = as.numeric(xml2::xml_attr(off, "x")), y = as.numeric(xml2::xml_attr(off, "y")), w = as.numeric(xml2::xml_attr(ext, 
        "cx")), h = as.numeric(xml2::xml_attr(ext, "cy")))/360000
    axis <- which(tx == "HR (95% CI)")
    stopifnot(length(axis) == 4)
    drop <- integer()
    top <- numeric(1)
    for (j in 2L) {
        last <- axis[j]
        nextplot <- which(seq_along(sp) > last & xy$w > 15 & xy$h > 1)[1]
        drop <- c(drop, seq.int(last + 1L, nextplot - 1L))
        first <- max(which(seq_along(sp) < last & xy$w > 15 & xy$h > 1))
        a <- xy$y[first]
        b <- a + xy$h[first]
        top[j - 1L] <- a - 0.13
        shift <- 0.45
        scale <- (b - a - shift)/(b - a)
        for (i in seq.int(first, last)) {
            yy <- a + shift + (xy$y[i] - a) * scale
            geom <- xml2::xml_attr(xml2::xml_find_first(sp[[i]], "./p:spPr/a:prstGeom", ns), "prst")
            fixed <- nzchar(tx[i]) || identical(geom, "ellipse")
            hh <- if (fixed) 
                xy$h[i]
            else xy$h[i] * scale
            if (fixed) 
                yy <- yy + (xy$h[i] * scale - hh)/2
            xml2::xml_set_attr(off[[i]], "y", sprintf("%.0f", yy * 360000))
            xml2::xml_set_attr(ext[[i]], "cy", sprintf("%.0f", hh * 360000))
        }
    }
    i <- which(seq_along(sp) > axis[4] & tx == "Not significant")[1]
    boxes <- which(seq_along(sp) > axis[4] & seq_along(sp) < i & xy$w > 1 & xy$h < 0.3 & tx == "")
    drop <- c(drop, seq.int(max(boxes), length(sp)))
    xml2::xml_remove(sp[unique(drop)])
    xml2::write_xml(doc, path)
    basefile <- tempfile(fileext = ".pptx")
    on.exit(unlink(basefile), add = TRUE)
    zip::zipr(basefile, list.files(tmp, recursive = TRUE, all.files = TRUE, no.. = TRUE), root = tmp, mode = "mirror")
    deck <- read_pptx(basefile)
    draw <- function() {
        label <- function(s, x, y, size = 6, col = "#30363B", just = "centre") grid::grid.text(s, unit(x, "cm"), 
            unit(14 - y, "cm"), just = just, gp = grid::gpar(fontfamily = font, fontsize = size, col = col))
        dot <- function(x, y, col) grid::grid.points(unit(x, "cm"), unit(14 - y, "cm"), pch = 16, size = unit(1.6, 
            "mm"), gp = grid::gpar(col = col))
        routes <- list(HH = c(0.21, 0.21), HL = c(0.28, 0.44), LH = c(0.44, 0.28), LL = c(0.51, 0.51))
        ribbon <- function(x, y, k, w, col) {
            t <- seq(0, 1, length.out = 65)
            xx <- x + 1.6 * (3 * (1 - t)^2 * t * 0.45 + 3 * (1 - t) * t^2 * 0.55 + t^3)
            yy <- y + routes[[k]][1] + diff(routes[[k]]) * (3 * t^2 - 2 * t^3)
            grid::grid.polygon(unit(c(xx, rev(xx)), "cm"), unit(14 - c(yy - w/2, rev(yy + w/2)), "cm"), gp = grid::gpar(fill = col, 
                col = NA))
        }
        colours <- list(c("#8FB8D8", "#6497C4", "#3F78B0", "#1F5A91"), strat_spec()$tdi$colors)
        keys <- list(c("LH", "HH", "HL", "LH"), c("LL", "HL", "LH", "HH"))
        refs <- c("HL", "LL", "HH", "LL")
        for (j in 1L) {
            y <- top[j] + 0.28
            label("Birthplace TDI  →  Recruitment TDI", 11.45, y + 0.015)
            label("High", 6.93, y + 0.245, just = "right")
            label("Low", 6.93, y + 0.475, just = "right")
            label("High", 16.47, y + 0.245, just = "left")
            label("Low", 16.47, y + 0.475, just = "left")
            if (j == 1) {
                label("TDI group contrasts", 3.4, y + 0.12, just = "left")
                dot(3.4, y + 0.4, "#B2182B")
                label(paste("Per SD", dtdi_lab), 3.58, y + 0.4, just = "left")
            }
            else label(paste("MD by", dtdi_lab), 3.4, y + 0.36, just = "left")
            for (i in 1:4) {
                x <- 7.1 + (i - 1) * 2.4
                for (k in names(routes)) ribbon(x, y, k, 0.018, "#E5E7EB")
                if (j == 1) 
                  ribbon(x, y, refs[i], 0.06, "#989EA5")
                ribbon(x, y, keys[[j]][i], 0.06, colours[[j]][i])
                for (xx in c(x, x + 1.6)) for (yy in c(0.245, 0.475)) grid::grid.rect(unit(xx, "cm"), unit(14 - 
                  y - yy, "cm"), unit(0.025, "cm"), unit(0.13, "cm"), gp = grid::gpar(fill = "#727B83", col = NA))
                if (j == 1) {
                  yy <- y + routes[[refs[i]]][1]
                  grid::grid.rect(unit(x + 0.25, "cm"), unit(14 - yy, "cm"), unit(0.33, "cm"), unit(0.2, "cm"), 
                    gp = grid::gpar(fill = "white", col = NA))
                  label("Ref", x + 0.25, yy, col = "#5E646B")
                }
            }
        }
    }
    deck <- ph_with(deck, rvg::dml(code = draw(), bg = "transparent"), location = ph_location(left = 0, top = 0, 
        width = 21 * cm, height = 14 * cm))
    print(deck, target = basefile)
    unzip(basefile, exdir = tmp)
    doc <- xml2::read_xml(path)
    group <- tail(xml2::xml_find_all(doc, ".//p:grpSp/p:nvGrpSpPr/p:cNvPr", ns), 1)
    xml2::xml_set_attr(group, "name", "route_keys")
    for (s in xml2::xml_find_all(doc, ".//p:sp[.//a:t]", ns)) {
        text <- xml2::xml_find_all(s, ".//a:rPr[@sz]", ns)
        if (!length(text)) 
            next
        ratio <- 600/as.numeric(xml2::xml_attr(text[[1]], "sz"))
        off <- xml2::xml_find_first(s, "./p:spPr/a:xfrm/a:off", ns)
        ext <- xml2::xml_find_first(s, "./p:spPr/a:xfrm/a:ext", ns)
        for (v in c("x", "y")) {
            w <- as.numeric(xml2::xml_attr(ext, paste0("c", v)))
            a <- as.numeric(xml2::xml_attr(off, v))
            xml2::xml_set_attr(off, v, sprintf("%.0f", a + w * (1 - ratio)/2))
            xml2::xml_set_attr(ext, paste0("c", v), sprintf("%.0f", w * ratio))
        }
        xml2::xml_set_attr(text, "sz", "600")
    }
    xml2::write_xml(doc, path)
    zip::zipr(basefile, list.files(tmp, recursive = TRUE, all.files = TRUE, no.. = TRUE), root = tmp, mode = "mirror")
    if (!file.copy(basefile, file, overwrite = TRUE)) 
        stop("Close the PowerPoint file: ", file)
    invisible(file)
}

run_fig2b <- function() {
    p <- fig2b_plot(plot_data$result$ses) + plot_annotation(tag_levels = list("b"), theme = theme(plot.tag = element_text(family = font,
        size = 7, face = "bold"), plot.tag.position = c(0.01, 0.99)))
    write_ppt(file.path(out, "fig2b.pptx"), list(p), 9 * cm, 9 * cm)
}

run_fig2c <- function() {
    r <- plot_data$result$forest
    s <- plot_data$result$strata
    direction <- plot_data$result$direction
    c1 <- c(`MD Q1` = "#E9D9A4", `MD Q2` = "#C5C66A", `MD Q3` = "#83A94C", `MD Q4` = "#3F7F45", MD = "#C18543", 
        `Migration cross TTWA` = "#8F65B7")
    c2 <- c(`Per SD ΔTDI` = "#B2182B", `Low to high vs High to low` = "#8FB8D8", `High to high vs Low to low` = "#6497C4", 
        `High to low vs High to high` = "#3F78B0", `Low to high vs Low to low` = "#1F5A91")
    lim <- function(id) {
        v <- unlist(s[strat == id, .(lo, hi)])
        v <- c(v[is.finite(v) & v > 0 & v < 3], 1)
        z <- range(v)
        z + c(-1, 1) * max(0.02, diff(z) * 0.08)
    }
    p <- forest_panel(r$main, c1, FALSE)/forest_panel(r$tdi, c2, FALSE)/direction_plot(direction, FALSE)/fig2c_strat_panel(s, 
        "work", lim("work"), TRUE) + plot_layout(heights = rep(1, 4)) + plot_annotation(tag_levels = list("c"), 
        theme = theme(plot.tag = element_text(family = font, size = 7, face = "bold"), plot.tag.position = c(0.003, 
            0.998))) & theme(plot.margin = margin(0.5, 1, 0.5, 1))
    write_ppt(file.path(out, "fig2c.pptx"), list(p), 21 * cm, 14 * cm)
    fig2c_keys(file.path(out, "fig2c.pptx"))
    invisible(p)
}

direction_plot <- function(r, show_x = TRUE) {
    labels <- c("Density shift (per SD)", paste("Urban", intToUtf8(8594), "Urban (Ref.)"), paste("Urban", intToUtf8(8594), 
        "Rural"), paste("Rural", intToUtf8(8594), "Urban"), paste("Rural", intToUtf8(8594), "Rural"))
    key <- setNames(labels, c("density_z", "urban_reference", paste0("urban_route", c("Urban to Rural", "Rural to Urban", 
        "Rural to Rural"))))
    d <- direction_table(r)
    d[, `:=`(term, unname(key[term]))]
    cols <- setNames(c(pal$density, "#F0C09A", "#E99559", "#D4713B", "#AE4E25"), labels)
    forest_panel(d, cols, show_x) + guides(shape = "none", color = guide_legend(nrow = 1, override.aes = list(shape = c(16, 
        1, 16, 16, 16), linewidth = 0))) + theme(legend.background = element_blank(), legend.key.width = unit(0.2, 
        "cm"), legend.spacing.x = unit(0.06, "cm"))
}

fig2c_strat_panel <- function(d, id, lim, show_x = TRUE) {
    s <- strat_spec()[[id]]
    lab <- c("All-cause\nmortality", "CVD\nmortality", "CAD", "Stroke", "HF", "AF", "T2D", "Depression", "Anxiety", 
        "Dementia")
    x <- copy(d[strat == id])
    x[, `:=`(outcome = factor(outcome, levels = forest_out), group = factor(group, levels = s$levels), sig = is.finite(FDR) & 
        FDR < 0.05)]
    off <- setNames(seq(-0.24, 0.24, length.out = length(s$levels)), s$levels)
    x[, `:=`(xpos, as.numeric(outcome) + off[as.character(group)])]
    ggplot(x, aes(xpos, HR, color = group)) + geom_hline(yintercept = 1, linewidth = 0.14, color = "#666666") + 
        geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.24) + geom_point(aes(shape = sig), size = 1.6, 
        stroke = 0.45) + scale_color_manual(values = setNames(s$colors, s$levels), name = if (id == "tdi") 
        "MD by TDI group"
    else "MD by physical demand", drop = FALSE) + scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = NULL, 
        labels = c("Not significant", "FDR < 0.05")) + scale_x_continuous(NULL, breaks = seq_along(forest_out), 
        labels = lab) + scale_y_continuous("HR (95% CI)") + coord_cartesian(ylim = c(lim[1], lim[2] + 0.18 * diff(lim))) + 
        labs(title = NULL) + theme_classic(base_family = font, base_size = 5.8) + theme(axis.text.x = if (show_x) 
        element_text(size = 6, color = "#222222")
    else element_blank(), axis.ticks.x = if (show_x) 
        element_line(linewidth = 0.14)
    else element_blank(), axis.text.y = element_text(size = 5.2, color = "#222222"), axis.title.y = element_text(size = 5.8), 
        axis.line = element_line(linewidth = 0.14), axis.ticks = element_line(linewidth = 0.14), legend.position = "inside", 
        legend.position.inside = c(0.5, 0.99), legend.justification = c(0.5, 1), legend.direction = "horizontal", 
        legend.box = "horizontal", legend.background = element_rect(fill = "white", color = NA), legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), legend.key.width = unit(0.13, "in"), legend.key.height = unit(0.08, 
            "in"), legend.spacing.x = unit(0.08, "cm"), legend.box.spacing = unit(0.1, "cm"), legend.margin = margin(0, 
            0, 0, 0), plot.margin = margin(0.5, 2, 0.5, 2)) + guides(color = guide_legend(order = 1, nrow = 1, 
        byrow = TRUE, title.position = "left", override.aes = list(shape = 16, linewidth = 0)), shape = guide_legend(order = 2, 
        nrow = 1, override.aes = list(color = "#222222", size = 1.6)))
}

int_term_level <- function(x, term) {
    y <- gsub("md_z|dtdi_z|age_z|sex2|tdi4|work_g|urban_route|:", "", term)
    fifelse(x %chin% c("MD x TDI group", "MD x physical demand", "MD x sex", "dTDI x sex", "MD x rural-urban route"), 
        y, "")
}

int_heat_data <- function(z, tests = NULL) {
    rows <- data.table(interaction = c("MD x age", "dTDI x age", "MD x sex", "dTDI x sex", "MD x physical demand", 
        "MD x dTDI", "MD x density change", rep("MD x rural-urban route", 3)), level = c("", "", "Male", "Male", 
        "Usually", "", "", "Urban to Rural", "Rural to Urban", "Rural to Rural"), label = c("MD × age", "ΔTDI × age", 
        "MD × men", "ΔTDI × men", "MD × Occupational\nphysical demand: usually", "MD × ΔTDI", "MD × density change", 
        "MD × Urban → Rural", "MD × Rural → Urban", "MD × Rural → Rural"), block = rep(c("Health outcomes", 
        "SES"), each = 5))
    rows[, `:=`(key, paste(interaction, level, sep = "|"))]
    d0 <- copy(z[model %chin% c("Cox", "Linear")])
    if (!is.null(tests)) 
        d0 <- merge(d0, tests[, .(outcome, interaction, omnibus_q = FDR)], by = c("outcome", "interaction"), all.x = TRUE, 
            sort = FALSE)
    else d0[, `:=`(omnibus_q, 0)]
    d0[, `:=`(eligible, is.finite(FDR) & FDR < 0.05 & is.finite(omnibus_q) & omnibus_q < 0.05)]
    d0[, `:=`(key, paste(interaction, int_term_level(interaction, term), sep = "|"))]
    stopifnot(!anyDuplicated(d0[, .(outcome, key)]))
    h <- rows[block == "Health outcomes"]
    s <- rows[block == "SES"]
    health <- forest_out[forest_out %chin% d0[model == "Cox" & eligible == TRUE, unique(outcome)]]
    cells <- rbindlist(list(CJ(outcome = health, key = h$key), CJ(outcome = "SES", key = s$key)))
    d <- merge(cells, d0, by = c("outcome", "key"), all.x = TRUE, sort = FALSE)
    d <- merge(d, rows[, .(key, label, block)], by = "key", all.x = TRUE, sort = FALSE)
    d[, `:=`(x = fifelse(block == "SES", 3.18, 3 + (match(outcome, health) - 0.5) * 0.55), y = fifelse(block == 
        "SES", 1.88 - (match(key, s$key) - 1) * 0.38, 6.15 - 0.275 - (match(key, h$key) - 1) * 0.55 - fifelse(key == 
        h$key[5], 0.12, 0)), cell = fifelse(block == "SES", 0.36, 0.55), score = fifelse(eligible, pmin(-log10(FDR), 
        12), NA_real_))]
    d
}

int_heat <- function(z, tests = NULL) {
    d <- int_heat_data(z, tests)
    children <- list()
    add <- function(x) children[[length(children) + 1L]] <<- x
    txt <- function(label, x, y, size = 6, just = "centre", colour = "#252A30", rot = 0, ...) {
        grid::textGrob(label, x = unit(x, "native"), y = unit(y, "native"), just = just, rot = rot, gp = grid::gpar(fontfamily = font, 
            fontsize = size, col = colour, ...))
    }
    rect <- function(x, y, w, h, fill = "white", colour = "#E2E6EA", lwd = 0.28) {
        grid::rectGrob(x = unit(x, "native"), y = unit(y, "native"), width = unit(w, "native"), height = unit(h, 
            "native"), gp = grid::gpar(fill = fill, col = colour, lwd = lwd))
    }
    circle <- function(x, y, q, fill) {
        grid::circleGrob(x = unit(x, "native"), y = unit(y, "native"), r = unit(0.16 * sqrt(pmin(-log10(q), 12)/12), 
            "native"), gp = grid::gpar(fill = fill, col = "white", lwd = 0.35))
    }
    lim <- log(c(0.9, 1.18))
    ramp <- grDevices::colorRamp(c("#2E6791", "#F7F9FA", "#BE6550"), space = "Lab")
    shade <- function(b) {
        p <- ifelse(b < 0, 0.5 * (b - lim[1])/-lim[1], 0.5 + 0.5 * b/lim[2])
        grDevices::rgb(ramp(pmax(0, pmin(1, p))), maxColorValue = 255)
    }
    h <- unique(d[block == "Health outcomes", .(outcome, x)])[order(x)]
    add(txt("Health outcomes", mean(h$x), 6.44, 7))
    add(txt("SES", 3.18, 2.16, 7, colour = "#73518F"))
    add(txt(expression(italic(beta)[int]), 3.86, 2.16, 7, colour = "#73518F"))
    labels <- unique(d[, .(label, y, block)])
    for (i in seq_len(nrow(labels))) {
        a <- labels[i]
        add(txt(a$label, 2.85, a$y, 6.3, just = "right", lineheight = 0.95))
    }
    for (i in seq_len(nrow(d))) {
        a <- d[i]
        ses <- a$block == "SES"
        add(rect(a$x, a$y, a$cell, a$cell, fill = if (ses) 
            "#FCFAFE"
        else "white", colour = if (ses) 
            "#DAD3E1"
        else "#E2E6EA"))
        if (!isTRUE(a$eligible)) 
            next
        add(circle(a$x, a$y, a$FDR, if (ses) 
            "#82609F"
        else shade(log(a$ratio))))
        if (ses) {
            label <- paste0(if (a$beta < 0) 
                intToUtf8(8722)
            else "+", sprintf("%.3f", abs(a$beta)))
            add(txt(label, 3.86, a$y, 6.3, colour = "#654B7C"))
        }
    }
    short <- c(`All-cause mortality` = "All-cause", `CVD mortality` = "CVD death", `CAD/MI` = "CAD", `Heart failure` = "HF", 
        `Type 2 diabetes` = "T2D", Depression = "Depression", `Anxiety disorder` = "Anxiety", Dementia = "Dementia")
    for (i in seq_len(nrow(h))) add(txt(short[h$outcome[i]], h$x[i] + 0.04, 3.17, 6, just = "right", rot = 45))
    add(txt("RHR", 6.1, 2.15, 6.3))
    z0 <- seq(lim[1], lim[2], length.out = 257)
    xx <- 5 + (seq_len(256) - 0.5)/256 * 2.2
    for (i in seq_len(256)) add(rect(xx[i], 1.82, 2.2/256 * 1.02, 0.13, fill = shade(mean(z0[c(i, i + 1)])), colour = NA, 
        lwd = 0))
    for (v in c(0.9, 1, 1.18)) {
        at <- 5 + (log(v) - lim[1])/diff(lim) * 2.2
        add(txt(sprintf("%.2f", v), at, 1.55, 6))
    }
    add(txt(expression(-log[10](italic(q))), 6.1, 1.14, 6.3))
    for (i in 1:3) {
        q <- c(2, 6, 12)[i]
        at <- c(5.25, 6.1, 6.95)[i]
        add(circle(at, 0.75, 10^-q, "#ADB5BC"))
        add(txt(c("2", "6", "12+")[i], at, 0.42, 6))
    }
    grob <- grid::grobTree(children = do.call(grid::gList, children), vp = grid::viewport(xscale = c(0, 8), yscale = c(0, 
        6.65), clip = "off"))
    cowplot::ggdraw() + cowplot::draw_grob(grob)
}

plot_interaction <- function() {
    device <- cowplot::set_null_device("agg")
    on.exit(cowplot::set_null_device(device), add = TRUE)
    r <- plot_data$result$interaction
    t <- copy(r$Tests)
    z <- copy(r$Terms)
    e <- copy(r$Effects)
    bounds <- plot_data$bounds
    selected <- plot_data$selected
    g <- lapply(seq_len(nrow(selected)), function(i) {
        a <- selected[i]
        int_curve(e, t, a$interaction, a$outcome, ed_out(a$outcome), a$type, bounds)
    })
    p <- ggdraw() + draw_plot(int_heat(z, t), 0, 16.85/24, 1, 6.65/24)
    for (i in seq_along(g)) {
        col <- (i - 1)%%2
        row <- (i - 1)%/%2
        p <- p + draw_plot(g[[i]], col/2, (0.25 + (3 - row) * 4.15)/24, 1/2, 4/24)
    }
    p <- p + draw_label("d", x = 0.005, y = 0.995, hjust = 0, vjust = 1, fontfamily = font, fontface = "bold", 
        size = 10)
    write_ppt(file.path(out, "fig2d.pptx"), list(p), 8 * cm, 24 * cm)
    ragg::agg_png(file.path(out, "fig2d.png"), width = 8, height = 24, units = "cm", res = 600, background = "white")
    print(p)
    dev.off()
    invisible(p)
}

run_fig2 <- function() {
    src <- file.path(out, paste0("fig2", letters[1:4], ".pptx"))
    stopifnot(all(file.exists(src)))
    tmp <- tempfile("fig2_")
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    ns <- c(p = "http://schemas.openxmlformats.org/presentationml/2006/main", a = "http://schemas.openxmlformats.org/drawingml/2006/main", 
        r = "http://schemas.openxmlformats.org/officeDocument/2006/relationships")
    dirs <- file.path(tmp, letters[1:4])
    for (i in 1:4) unzip(src[i], exdir = dirs[i])
    size <- t(vapply(dirs, function(dir) {
        doc <- xml2::read_xml(file.path(dir, "ppt/presentation.xml"))
        z <- xml2::xml_find_first(doc, "//p:sldSz", ns)
        as.numeric(xml2::xml_attrs(z)[c("cx", "cy")])
    }, numeric(2)))
    slides <- lapply(dirs, function(dir) xml2::read_xml(file.path(dir, "ppt/slides/slide1.xml")))
    invisible(lapply(slides, fig2_text))
    shift <- matrix(0, 4, 2)
    bounds <- matrix(0, 4, 4)
    for (i in 1:4) {
        shapes <- xml2::xml_find_all(slides[[i]], "//p:sp", ns)
        tags <- vapply(shapes, function(z) paste(xml2::xml_text(xml2::xml_find_all(z, ".//a:t", ns)), collapse = ""), 
            "")
        xml2::xml_remove(shapes[tags == letters[i]])
        shapes <- xml2::xml_find_all(slides[[i]], "//p:sp | //p:pic", ns)
        extent <- lapply(shapes, function(z) {
            tr <- xml2::xml_find_first(z, "p:spPr/a:xfrm", ns)
            xy <- as.numeric(xml2::xml_attrs(xml2::xml_find_first(tr, "a:off", ns))[c("x", "y")])
            wh <- as.numeric(xml2::xml_attrs(xml2::xml_find_first(tr, "a:ext", ns))[c("cx", "cy")])
            fill <- xml2::xml_attr(xml2::xml_find_first(z, "p:spPr/a:solidFill/a:srgbClr", ns), "val")
            line <- xml2::xml_find_all(z, "p:spPr/a:ln/a:solidFill/*", ns)
            blank <- !length(xml2::xml_find_all(z, ".//a:t", ns)) && fill %in% c("FFFFFF", "ffffff") && (!length(line) || 
                all(xml2::xml_attr(line, "val") %in% c("FFFFFF", "ffffff"))) && all(wh > size[i, ] * 0.1)
            if (blank) {
                xml2::xml_remove(z)
                return(NULL)
            }
            angle <- as.numeric(xml2::xml_attr(tr, "rot"))/60000 * pi/180
            if (!is.finite(angle)) 
                angle <- 0
            half <- c(abs(wh[1] * cos(angle)) + abs(wh[2] * sin(angle)), abs(wh[2] * cos(angle)) + abs(wh[1] * 
                sin(angle)))/2
            c(xy + wh/2 - half, xy + wh/2 + half)
        })
        extent <- do.call(rbind, extent)
        bounds[i, ] <- c(apply(extent[, 1:2], 2, min, na.rm = TRUE), apply(extent[, 3:4], 2, max, na.rm = TRUE))
    }
    shift <- bounds[, 1:2]
    size <- bounds[, 3:4] - shift
    for (i in 1:4) {
        groups <- xml2::xml_find_all(slides[[i]], "//p:grpSp/p:grpSpPr/a:xfrm", ns)
        for (g in groups) {
            off <- xml2::xml_find_first(g, "a:off", ns)
            ch <- xml2::xml_find_first(g, "a:chOff", ns)
            ex <- xml2::xml_find_first(g, "a:ext", ns)
            ce <- xml2::xml_find_first(g, "a:chExt", ns)
            stopifnot(identical(xml2::xml_attrs(off), xml2::xml_attrs(ch)), identical(xml2::xml_attrs(ex), xml2::xml_attrs(ce)))
            for (n in list(off, ch)) xml2::xml_set_attrs(n, setNames(sprintf("%.0f", shift[i, ]), c("x", "y")))
            for (n in list(ex, ce)) xml2::xml_set_attrs(n, setNames(sprintf("%.0f", size[i, ]), c("cx", "cy")))
        }
    }
    pad <- c(left = 0.38, top = 0.45, right = 0.08, bottom = 0.08) * 360000
    px <- sum(pad[c("left", "right")])
    py <- sum(pad[c("top", "bottom")])
    h <- size[1, 2] + py
    w <- c(size[1, 1], size[1, 2] * size[2, 1]/size[2, 2]) + px
    left <- sum(w)
    bottom <- (left - px) * size[3, 2]/size[3, 1] + py
    height <- h + bottom
    right <- (height - py) * size[4, 1]/size[4, 2] + px
    layout <- data.frame(x = c(0, w[1], 0, left + 0.25 * 360000), y = c(0, 0, h, 0), w = c(w, left, right), h = c(h, 
        h, bottom, height))
    base <- file.path(tmp, "base")
    unzip(src[1], exdir = base)
    file <- file.path(base, "ppt/slides/slide1.xml")
    doc <- xml2::read_xml(file)
    tree <- xml2::xml_find_first(doc, "//p:spTree", ns)
    xml2::xml_remove(xml2::xml_children(tree)[-(1:2)])
    relfile <- file.path(base, "ppt/slides/_rels/slide1.xml.rels")
    rel <- xml2::read_xml(relfile)
    for (i in 1:4) {
        s <- slides[[i]]
        rr <- xml2::read_xml(file.path(dirs[i], "ppt/slides/_rels/slide1.xml.rels"))
        for (r in xml2::xml_children(rr)) {
            if (grepl("/slideLayout$", xml2::xml_attr(r, "Type"))) 
                next
            stopifnot(grepl("/image$", xml2::xml_attr(r, "Type")))
            old <- xml2::xml_attr(r, "Id")
            id <- paste0("rId", i, "_", old)
            target <- xml2::xml_attr(r, "Target")
            name <- paste0(letters[i], "_", basename(target))
            dir.create(file.path(base, "ppt/media"), showWarnings = FALSE)
            stopifnot(file.copy(file.path(dirs[i], "ppt/slides", target), file.path(base, "ppt/media", name), overwrite = TRUE))
            xml2::xml_set_attr(r, "Id", id)
            xml2::xml_set_attr(r, "Target", paste0("../media/", name))
            xml2::xml_add_child(rel, r)
            nodes <- xml2::xml_find_all(s, paste0("//*[@r:embed='", old, "']"), ns)
            xml2::xml_set_attr(nodes, "r:embed", id)
        }
        at <- layout[i, ]
        at$x <- at$x + pad["left"]
        at$y <- at$y + pad["top"]
        at$w <- at$w - px
        at$h <- at$h - py
        group <- xml2::read_xml(sprintf(paste0("<p:grpSp xmlns:p=\"%s\" xmlns:a=\"%s\"><p:nvGrpSpPr><p:cNvPr id=\"1\" name=\"panel_%s\"/>", 
            "<p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr><p:grpSpPr><a:xfrm>", "<a:off x=\"%.0f\" y=\"%.0f\"/><a:ext cx=\"%.0f\" cy=\"%.0f\"/>", 
            "<a:chOff x=\"%.0f\" y=\"%.0f\"/><a:chExt cx=\"%.0f\" cy=\"%.0f\"/>", "</a:xfrm></p:grpSpPr></p:grpSp>"), 
            ns["p"], ns["a"], letters[i], at$x, at$y, at$w, at$h, shift[i, 1], shift[i, 2], size[i, 1], size[i, 
                2]))
        nodes <- xml2::xml_children(xml2::xml_find_first(s, "//p:spTree", ns))[-(1:2)]
        for (node in nodes) xml2::xml_add_child(group, node)
        xml2::xml_add_child(tree, group)
    }
    for (i in 1:4) {
        tag <- xml2::read_xml(sprintf(paste0("<p:sp xmlns:p=\"%s\" xmlns:a=\"%s\"><p:nvSpPr><p:cNvPr id=\"1\" name=\"label_%s\"/>", 
            "<p:cNvSpPr txBox=\"1\"/><p:nvPr/></p:nvSpPr><p:spPr><a:xfrm>", "<a:off x=\"%.0f\" y=\"%.0f\"/><a:ext cx=\"180000\" cy=\"144000\"/></a:xfrm>", 
            "<a:prstGeom prst=\"rect\"><a:avLst/></a:prstGeom><a:noFill/><a:ln><a:noFill/></a:ln></p:spPr>", "<p:txBody><a:bodyPr lIns=\"0\" rIns=\"0\" tIns=\"0\" bIns=\"0\" anchor=\"t\" wrap=\"none\"/>", 
            "<a:lstStyle/><a:p><a:pPr algn=\"l\"/><a:r><a:rPr sz=\"1000\" b=\"1\">", "<a:solidFill><a:srgbClr val=\"000000\"/></a:solidFill><a:latin typeface=\"Arial\"/>", 
            "<a:cs typeface=\"Arial\"/></a:rPr><a:t>%s</a:t></a:r></a:p></p:txBody></p:sp>"), ns["p"], ns["a"], 
            letters[i], layout$x[i] + 28800, layout$y[i] + 28800, letters[i]))
        xml2::xml_add_child(tree, tag)
    }
    ids <- xml2::xml_find_all(doc, "//p:cNvPr", ns)
    xml2::xml_set_attr(ids, "id", as.character(seq_along(ids)))
    xml2::write_xml(doc, file)
    xml2::write_xml(rel, relfile)
    pfile <- file.path(base, "ppt/presentation.xml")
    pres <- xml2::read_xml(pfile)
    xml2::xml_set_attrs(xml2::xml_find_first(pres, "//p:sldSz", ns), c(cx = sprintf("%.0f", layout$x[4] + right), 
        cy = sprintf("%.0f", height)))
    xml2::write_xml(pres, pfile)
    target <- file.path(tmp, "Fig2.pptx")
    zip::zipr(target, list.files(base, recursive = TRUE, all.files = TRUE, no.. = TRUE), root = base, mode = 'mirror')
    stopifnot(file.copy(target, file.path(out, "Fig2.pptx"), overwrite = TRUE))
    invisible(layout/360000)
}

fig2a_draw <- function(data) {
    stopifnot(inherits(data, "move_fig2a_data"), length(data$labs) == 13L)
    labs <- data$labs
    d <- data$d
    h <- data$h
    font <- "Arial"
    cell <- 0.6
    side <- length(labs) * cell
    x0 <- 2.8
    y0 <- 0.8
    bx <- x0 + side + 0.1
    bw <- 2.15
    W <- bx + bw + 0.25
    H <- y0 + side + 2.55
    ramp <- scales::gradient_n_pal(c("#8F65B7", "#FFFFFF", "#00A087"), values = c(0, 0.5, 1))
    colour <- function(v) if (is.na(v)) 
        "#FFFFFF"
    else ramp((max(-1, min(1, v)) + 1)/2)
    text <- function(label, x, y, size = 6, just = c("left", "top"), rot = 0) {
        grid::grid.text(label, x = grid::unit(x, "cm"), y = grid::unit(H - y, "cm"), just = just, rot = rot, gp = grid::gpar(fontfamily = font, 
            fontsize = size, col = "#222222"))
    }
    rect <- function(x, y, w, h, fill, border = NA, lwd = 0.35) {
        grid::grid.rect(x = grid::unit(x, "cm"), y = grid::unit(H - y, "cm"), width = grid::unit(w, "cm"), height = grid::unit(h, 
            "cm"), just = c("left", "top"), gp = grid::gpar(fill = fill, col = border, lwd = lwd))
    }
    line <- function(x1, y1, x2, y2, col = "#333333", lwd = 0.45) {
        grid::grid.lines(x = grid::unit(c(x1, x2), "cm"), y = grid::unit(H - c(y1, y2), "cm"), gp = grid::gpar(col = col, 
            lwd = lwd))
    }
    draw <- function() {
        grid::grid.newpage()
        grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
        lx <- x0 + 3.15
        lw <- side - 3.15
        text(expression("Phenotypic " * italic(r) * " / LDSC " * italic(r)[g]), x0, 0.2, 6)
        for (k in 1:300) rect(lx + (k - 1) * lw/300, 0.16, lw/300, 0.14, ramp((k - 0.5)/300))
        for (v in c(-1, -0.5, 0, 0.5, 1)) {
            tg <- grid::textGrob(sprintf("%.1f", v), gp = grid::gpar(fontfamily = font, fontsize = 6))
            th <- grid::convertHeight(grid::grobHeight(tg), "cm", valueOnly = TRUE)
            text(sprintf("%.1f", v), lx + (v + 1) * lw/2, y0 - 0.1 - th, 6, c("centre", "top"))
        }
        for (k in seq_len(nrow(d))) {
            z <- d[k, ]
            xx <- x0 + (z$col - 1) * cell
            yy <- y0 + (z$row - 1) * cell
            rect(xx, yy, cell, cell, "white", "#D9D9D9")
            if (z$row == z$col) 
                next
            if (z$half == "Phenotypic") 
                rect(xx, yy, cell, cell, colour(z$val), "white")
            if (z$half == "LDSC" && is.finite(z$val)) 
                grid::grid.circle(x = grid::unit(xx + cell/2, "cm"), y = grid::unit(H - yy - cell/2, "cm"), r = grid::unit(cell * 
                  0.45 * sqrt(min(abs(z$val), 1)), "cm"), gp = grid::gpar(fill = colour(z$val), col = NA))
            if (is.finite(z$val)) 
                text(z$label, xx + cell/2, yy + cell/2, 6, c("centre", "centre"))
        }
        lim <- max(h$hi) * 1.36
        ticks <- seq(0, floor(lim/0.05) * 0.05, 0.05)
        for (v in ticks) line(bx + bw * v/lim, y0, bx + bw * v/lim, y0 + side, "#E3E3E3", 0.3)
        for (j in seq_along(labs)) {
            yy <- y0 + (j - 0.5) * cell
            text(labs[j], x0 - 0.12, yy, 7, c("right", "centre"))
            text(labs[j], x0 + (j - 0.5) * cell, y0 + side + 0.1, 7, c("right", "centre"), 45)
            rect(bx, yy - cell * 0.31, bw * h$h2[j]/lim, cell * 0.62, "#01A6AC")
            lo <- bx + bw * h$lo[j]/lim
            hi <- bx + bw * h$hi[j]/lim
            line(lo, yy, hi, yy, "#015661")
            line(lo, yy - 0.06, lo, yy + 0.06, "#015661")
            line(hi, yy - 0.06, hi, yy + 0.06, "#015661")
            text(h$text[j], hi + 0.07, yy, 6, c("left", "centre"))
        }
        ticktop <- y0 + side + 0.1
        for (v in ticks) text(sprintf("%.2f", v), bx + bw * v/lim, ticktop, 6, c("centre", "top"))
        th <- grid::convertHeight(grid::grobHeight(grid::textGrob("0.00", gp = grid::gpar(fontfamily = font, fontsize = 6))), 
            "cm", valueOnly = TRUE)
        text(expression(italic(h)[SNP]^2), bx + bw/2, ticktop + th + 0.1, 7, c("centre", "top"))
    }
    attr(draw, "size_cm") <- c(width = W, height = H)
    draw
}

run_fig2a <- function(file = file.path(out, "fig2a.pptx"), data = plot_data$a) {
    if (!grepl("\\.pptx$", file, ignore.case = TRUE)) 
        stop("file must end in .pptx")
    if (!dir.exists(dirname(file))) 
        stop("Output directory must already exist")
    draw <- fig2a_draw(data)
    size <- attr(draw, "size_cm")
    draw_ppt(file, draw, size[1]/2.54, size[2]/2.54)
    invisible(data)
}

run_fig2a(data = plot_data$a)

run_fig2b()

run_fig2c()

plot_interaction()

run_fig2()

finish("Fig2")

