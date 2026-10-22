# 1. Figure 3
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

region_map_labels <- function(a, cen, i, ctx) {
    z <- merge(cen, a[id == i & context == ctx, .(code, beta)], by = "code")
    keep <- which(abs(z$beta) >= log(1.1))
    keep <- unique(c(keep, which.min(z$beta), which.max(z$beta)))
    rank <- order(abs(z$beta), decreasing = TRUE)
    if (length(keep) > 8) 
        keep <- head(rank[rank %in% keep], 8)
    if (length(keep) < 6) 
        keep <- head(unique(c(keep, rank)), 6)
    stopifnot(length(keep) >= 6, length(keep) <= 8)
    z[keep]
}

region_label_positions <- function(z, ch = 12.5) {
    z <- data.table::copy(z)
    z[, `:=`(py, ch * (1 - y/1230000))]
    data.table::setorder(z, py)
    z[, `:=`(shown, vapply(name, function(s) paste(strwrap(s, width = 19), collapse = "\n"), character(1)))]
    h <- 0.42 * lengths(strsplit(z$shown, "\n")) + 0.22
    yy <- pmax(z$py - 0.12, 3.55)
    if (nrow(z) > 1) 
        for (k in 2:nrow(z)) yy[k] <- max(yy[k], yy[k - 1] + h[k - 1])
    if (tail(yy + h, 1) > ch - 0.15) {
        yy[nrow(z)] <- ch - 0.15 - tail(h, 1)
        if (nrow(z) > 1) 
            for (k in (nrow(z) - 1):1) yy[k] <- min(yy[k], yy[k + 1] - h[k])
    }
    stopifnot(min(yy) > 2.8)
    z[, `:=`(label_y_cm = yy, label_height_cm = h)]
    z
}

region_figure <- function(d = region_read(), geometry = region_geometry()) {
    a <- data.table::copy(d$area)
    t <- data.table::copy(d$tests)
    poly <- geometry$poly
    cen <- geometry$centroids
    map_w <- 4.55
    cell_h <- 8.75
    first_w <- 5.8
    gx <- 0.16
    gy <- 0.28
    outer <- 0.22
    main_gap <- 0.25
    W <- 2 * outer + first_w + main_gap + 4 * map_w + 3 * gx
    H <- 2 * outer + 3 * cell_h + 2 * gy
    cols <- c(origin = "#018571", destination = "#7545A0")
    pal <- list(origin = c("#003C30", "#018571", "#F7F5ED", "#A6611A", "#543005"), destination = c("#A61B29", "#D96A73", 
        "#FAF6FA", "#A27BC7", "#63328A"))
    span <- log(1.3)
    nodes <- scales::rescale(c(-span, -0.075, 0, 0.075, span), from = c(-span, span))
    label <- c(`All-cause mortality` = "All-cause mortality", `CVD mortality` = "CVD mortality", `CAD/MI` = "CAD", 
        Stroke = "Stroke", `Heart failure` = "Heart failure", `Atrial fibrillation` = "Atrial fibrillation", `Type 2 diabetes` = "Type 2 diabetes", 
        Depression = "Depression", `Anxiety disorder` = "Anxiety", Dementia = "Dementia")
    qtext <- function(x) ifelse(x < 0.001, formatC(x, format = "e", digits = 1), sprintf("%.3f", x))
    vtext <- function(x) ifelse(x < 1e-04, formatC(x, format = "e", digits = 1), sprintf("%.4f", x))
    selected <- t[q < 0.05]
    groups <- selected[, .(contexts = list(intersect(c("origin", "destination"), context)), maps = .N), by = .(id, 
        outcome)]
    groups <- data.table::rbindlist(list(groups[maps == 2][order(id)], groups[maps == 1][order(-id)]))
    stopifnot(nrow(selected) == 12L, nrow(groups) == 7L, sum(groups$maps) == 12L, sum(groups$maps == 2L) == 5L, 
        sum(groups$maps == 1L) == 2L)
    layout <- data.table::rbindlist(lapply(seq_len(nrow(groups)), function(k) {
        pair <- groups$maps[k] == 2L
        data.table::data.table(id = groups$id[k], outcome = groups$outcome[k], panel = letters[k + 1L], row = if (pair) 
            ceiling(k/2)
        else 3L, contexts = list(groups$contexts[[k]]), cols = list(if (pair) if (k%%2L) 1:2 else 3:4 else if (k == 
            6L) 3L else 4L))
    }))
    spread <- function(v, lo = 25000, hi = 975000, gap = 70000) {
        o <- order(v)
        y <- pmax(v[o], lo)
        if (length(y) > 1L) 
            for (k in 2:length(y)) y[k] <- max(y[k], y[k - 1L] + gap)
        if (tail(y, 1) > hi) {
            y[length(y)] <- hi
            if (length(y) > 1L) 
                for (k in (length(y) - 1L):1L) y[k] <- min(y[k], y[k + 1L] - gap)
        }
        ans <- numeric(length(v))
        ans[o] <- y
        ans
    }
    map_labels <- function(endpoint, ctx) {
        z <- region_map_labels(a, cen, endpoint, ctx)
        z[, `:=`(shown, vapply(name, function(x) paste(strwrap(x, width = 15), collapse = "\n"), character(1)))]
        data.table::setorder(z, x)
        z[, `:=`(side, rep(c("left", "right"), c(floor(.N/2), ceiling(.N/2))))]
        z[, `:=`(label_y, spread(y)), by = side]
        z[, `:=`(label_x = ifelse(side == "left", 9000, 651000), line_x = ifelse(side == "left", 120000, 540000), 
            label_h = ifelse(side == "left", 0, 1))]
        z
    }
    map_plot <- function(endpoint, ctx) {
        s <- t[id == endpoint & context == ctx]
        z <- merge(poly, a[id == endpoint & context == ctx, .(code, beta)], by = "code", all.x = TRUE, sort = FALSE)
        lab <- map_labels(endpoint, ctx)
        ann <- data.table::data.table(x = 647000, y = c(1210000, 1170000), text = c(if (ctx == "origin") "Birthplace" else "Recruitment", 
            paste("Effect variance", vtext(s$variance))))
        annq <- data.table::data.table(x = 647000, y = 1130000, text = paste0("italic(q) == ", qtext(s$q)))
        bar0 <- 485000
        bar1 <- 647000
        evx <- bar0 + (bar1 - bar0) * s$events/s$n
        nlab <- data.table::data.table(x = 647000, y = c(1045000, 1005000), text = c(paste("N", format(s$n, big.mark = ",")), 
            paste("Events", format(s$events, big.mark = ","))), colour = c("#333333", cols[[ctx]]))
        ggplot2::ggplot(z, ggplot2::aes(x, y, group = piece, subgroup = ring, fill = beta)) + ggplot2::geom_polygon(rule = "evenodd", 
            colour = NA) + ggplot2::geom_point(data = lab, ggplot2::aes(x, y), inherit.aes = FALSE, size = 0.35, 
            colour = "#333333") + ggplot2::geom_segment(data = lab, ggplot2::aes(x = x, y = y, xend = line_x, yend = label_y), 
            inherit.aes = FALSE, linewidth = 0.12, colour = "#767676", lineend = "butt") + ggplot2::geom_text(data = lab, 
            ggplot2::aes(x = label_x, y = label_y, label = shown, hjust = label_h), inherit.aes = FALSE, family = font, 
            size = pt(5), lineheight = 0.82, colour = "#242424") + ggplot2::geom_text(data = ann, ggplot2::aes(x, 
            y, label = text), inherit.aes = FALSE, hjust = 1, vjust = 1, family = font, size = pt(5.7), lineheight = 0.92, 
            colour = cols[[ctx]]) + ggplot2::geom_text(data = annq, ggplot2::aes(x, y, label = text), inherit.aes = FALSE, 
            hjust = 1, vjust = 1, family = font, size = pt(5.7), parse = TRUE, colour = cols[[ctx]]) + ggplot2::annotate("rect", 
            xmin = bar0, xmax = evx, ymin = 1070000, ymax = 1085000, fill = cols[[ctx]], colour = NA) + ggplot2::annotate("rect", 
            xmin = evx, xmax = bar1, ymin = 1070000, ymax = 1085000, fill = grDevices::adjustcolor(cols[[ctx]], 
                alpha.f = 0.22), colour = NA) + ggplot2::geom_text(data = nlab, ggplot2::aes(x, y, label = text, 
            colour = colour), inherit.aes = FALSE, hjust = 1, vjust = 1, family = font, size = pt(5), lineheight = 0.9, 
            show.legend = FALSE) + ggplot2::scale_colour_identity() + ggplot2::scale_fill_gradientn(colours = pal[[ctx]], 
            values = nodes, limits = c(-span, span), oob = scales::squish, na.value = "#DCDCDC", guide = "none") + 
            ggplot2::coord_fixed(xlim = c(-1000, 660000), ylim = c(0, 1230000), expand = FALSE, clip = "on") + 
            ggplot2::theme_void(base_family = font) + ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
    }
    summary_plot <- function() {
        v <- data.table::copy(t)
        v[, `:=`(outcome = factor(outcome, levels = rev(forest_out)), yoff = ifelse(context == "origin", 0.13, 
            -0.13), sig = q < 0.05, strength = pmin(-log10(pmax(q, 1e-12)), 10))]
        ggplot2::ggplot(v, ggplot2::aes(y = as.numeric(outcome) + yoff, colour = context)) + ggplot2::geom_vline(xintercept = 0, 
            colour = "#9A9A9A", linetype = "dashed", linewidth = 0.18) + ggplot2::geom_segment(ggplot2::aes(x = 0, 
            xend = variance, yend = as.numeric(outcome) + yoff), linewidth = 0.55, lineend = "round") + ggplot2::geom_point(ggplot2::aes(x = variance, 
            size = strength, shape = sig), stroke = 0.5) + ggplot2::scale_colour_manual(values = cols, guide = "none") + 
            ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") + ggplot2::scale_size_continuous(range = c(1.8, 
            4.3), guide = "none") + ggplot2::scale_y_continuous(breaks = seq_along(levels(v$outcome)), labels = levels(v$outcome), 
            limits = c(0.45, 10.55), expand = c(0, 0)) + ggplot2::scale_x_continuous(limits = c(-3e-04, 0.014), 
            breaks = c(0, 0.004, 0.008, 0.012), labels = c("0", ".004", ".008", ".012"), expand = c(0, 0)) + ggplot2::labs(x = "Effect variance", 
            y = NULL) + ggplot2::theme_classic(base_family = font, base_size = 6) + ggplot2::theme(axis.title.x = ggplot2::element_text(size = 7), 
            axis.text = ggplot2::element_text(size = 6, colour = "#222222"), axis.text.y = ggplot2::element_text(hjust = 1), 
            axis.line = ggplot2::element_line(linewidth = 0.14), axis.ticks = ggplot2::element_line(linewidth = 0.14), 
            panel.grid.major.y = ggplot2::element_line(colour = "#ECECEC", linewidth = 0.16), plot.margin = ggplot2::margin(4, 
                3, 3, 2))
    }
    place <- function(p, x, y, w, h) {
        grid::pushViewport(grid::viewport(x = grid::unit(x, "cm"), y = grid::unit(H - y, "cm"), width = grid::unit(w, 
            "cm"), height = grid::unit(h, "cm"), just = c("left", "top")))
        on.exit(grid::popViewport())
        grid::grid.draw(ggplot2::ggplotGrob(p))
    }
    txt <- function(s, x, y, size = 6, bold = FALSE, just = "left", col = "#222222", italic = FALSE) grid::grid.text(s, 
        x = grid::unit(x, "cm"), y = grid::unit(H - y, "cm"), just = c(just, "top"), gp = grid::gpar(fontfamily = font, 
            fontsize = size, fontface = if (bold) 
                "bold"
            else if (italic) 
                "italic"
            else "plain", col = col))
    rect <- function(x, y, w, h, col) grid::grid.rect(x = grid::unit(x, "cm"), y = grid::unit(H - y, "cm"), width = grid::unit(w, 
        "cm"), height = grid::unit(h, "cm"), just = c("left", "top"), gp = grid::gpar(fill = col, col = NA))
    point <- function(x, y, col, pch = 16, size = 2.2) grid::grid.points(grid::unit(x, "cm"), grid::unit(H - y, 
        "cm"), pch = pch, size = grid::unit(size, "mm"), gp = grid::gpar(col = col, fill = col, lwd = 0.5))
    legend_panel <- function(x, y, w, h) {
        point(x + 0.18, y + 0.28, cols[["origin"]], 16, 2.4)
        txt("Birthplace", x + 0.43, y + 0.14, 6)
        point(x + 2.85, y + 0.28, cols[["destination"]], 16, 2.4)
        txt("Recruitment", x + 3.1, y + 0.14, 6)
        point(x + 0.18, y + 0.91, "#333333", 16, 2.3)
        txt("q < 0.05", x + 0.43, y + 0.77, 6, italic = TRUE)
        point(x + 2.85, y + 0.91, "#333333", 1, 2.3)
        txt("q >= 0.05", x + 3.1, y + 0.77, 6, italic = TRUE)
        ramp <- function(p) scales::gradient_n_pal(p, values = nodes)
        top <- y + 1.9
        height <- h - 2.7
        for (j in seq_along(c("origin", "destination"))) {
            ctx <- c("origin", "destination")[j]
            xx <- x + 0.72 + (j - 1) * w/2
            txt(if (ctx == "origin") 
                "Birthplace area HR"
            else "Recruitment area HR", xx + 0.36, y + 1.48, 6, just = "center")
            for (k in 1:200) rect(xx, top + height - k * height/200, 0.34, height/200 + 0.005, ramp(pal[[ctx]])((k - 
                0.5)/200))
            at <- top + height * (1 - scales::rescale(log(c(0.8, 1, 1.3)), from = c(-span, span)))
            for (k in seq_along(at)) txt(c(".8", "1", "1.3")[k], xx + 0.55, at[k] - 0.06, 5.7)
        }
        rect(x + 0.08, y + h - 0.4, 0.26, 0.26, "#DCDCDC")
        txt("No estimate", x + 0.46, y + h - 0.45, 6)
    }
    tables <- list(panel_a = data.table::copy(t)[, `:=`(panel = "a", strength = pmin(-log10(pmax(q, 1e-12)), 10), 
        significant = q < 0.05)])
    labels <- list()
    for (j in seq_len(nrow(layout))) {
        g <- layout[j]
        tables[[paste0("panel_", g$panel)]] <- data.table::rbindlist(lapply(g$contexts[[1]], function(ctx) {
            s <- t[id == g$id & context == ctx]
            lab <- map_labels(g$id, ctx)
            labels[[paste(g$panel, ctx, sep = "_")]] <<- lab
            z <- merge(cen, a[id == g$id & context == ctx], by = "code", all.x = TRUE, sort = FALSE)
            z <- merge(z, lab[, .(code, label_x, label_y, line_x, shown)], by = "code", all.x = TRUE, sort = FALSE)
            z[, `:=`(panel = g$panel, id = g$id, outcome = g$outcome, context = ctx, panel_variance = s$variance, 
                panel_q = s$q, panel_n = s$n, panel_events = s$events, has_estimate = !is.na(beta), display_beta = pmin(pmax(beta, 
                  -span), span), is_label = !is.na(shown))]
            z
        }), fill = TRUE)
    }
    source_data <- list(panels = tables, selected = selected, labels = labels, all_area_estimates = a, tests = t, 
        legacy_joint = d$joint, method = plot_data$result$method, cache = d$cache_provenance, layout = list(grid = layout,
            width_cm = W, height_cm = H, map_width_cm = map_w, cell_height_cm = cell_h, gap_x_cm = gx, gap_y_cm = gy, 
            xlim = c(-1000, 660000), ylim = c(0, 1230000), palette = pal, colour_nodes = nodes, log_hr_limits = c(-span, 
                span), no_estimate = "#DCDCDC", label_rule = "HR >= 1.10 or HR <= 0.91; retain eight strongest if needed; fill to six labels"))
    xmap <- outer + first_w + main_gap
    xc <- xmap + (0:3) * (map_w + gx)
    yr <- outer + (0:2) * (cell_h + gy)
    draw <- function() {
        grid::grid.newpage()
        grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
        place(summary_plot(), outer, outer + 0.22, first_w, 2 * cell_h + gy - 0.22)
        txt("a", outer, outer, 10, TRUE)
        legend_panel(outer, yr[3], first_w, cell_h)
        for (j in seq_len(nrow(layout))) {
            g <- layout[j]
            cs <- g$cols[[1]]
            contexts <- g$contexts[[1]]
            for (k in seq_along(cs)) place(map_plot(g$id, contexts[k]), xc[cs[k]], yr[g$row] + 0.22, map_w, cell_h - 
                0.22)
            x0 <- xc[min(cs)]
            x1 <- xc[max(cs)] + map_w
            txt(g$panel, x0, yr[g$row], 10, TRUE)
            txt(unname(label[g$outcome]), (x0 + x1)/2, yr[g$row] + 0.01, 7, TRUE, "center")
        }
        invisible(NULL)
    }
    list(draw = draw, source_data = source_data, geometry = geometry, width_cm = W, height_cm = H, selected = selected, 
        plots = list(summary = summary_plot, map = map_plot))
}

r <- region_figure(plot_data$result, plot_data$geography)

draw_ppt(file.path(out, "Fig3.pptx"), r$draw, r$width_cm * cm, r$height_cm * cm)

ragg::agg_png(file.path(out, "Fig3.png"), width = r$width_cm, height = r$height_cm, units = "cm", res = 300)

r$draw()

dev.off()

finish("Fig3")

