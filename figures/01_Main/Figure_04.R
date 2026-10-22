# 1. Figure 4
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

run_fig4 <- function(r = NULL) {
    if (is.null(r)) 
        r <- plot_data$result
    cols <- c(UKHLS = "#B76585", HRS = "#39A6B1")
    for (n in setdiff(names(r), "Methods")) r[[n]][, `:=`(cohort, factor(cohort, names(cols)))]
    cs <- list(scale_colour_manual(values = cols), scale_fill_manual(values = cols))
    th <- theme0() + theme(legend.position = "none", axis.text = element_text(size = 6, colour = "black"), plot.margin = margin(4, 
        5, 4, 5), axis.title = element_text(size = 7), legend.background = element_blank(), legend.key = element_blank())
    cowplot::set_null_device(function(width, height) ragg::agg_capture(width = width, height = height, units = "in", 
        res = 96))
    device <- ragg::agg_capture(width = 24, height = 18, units = "cm", res = 96)
    on.exit(dev.off(), add = TRUE)
    w <- copy(r$Waves)
    w[, `:=`(x, wave + ifelse(cohort == "UKHLS", -0.18, 0.18))]
    a0 <- ggplot(w, aes(x, n, fill = cohort)) + geom_col(width = 0.32) + scale_fill_manual(values = cols) + th + 
        scale_x_continuous("Wave", limits = range(w$wave) + c(-0.6, 0.6), breaks = pretty(w$wave, 4), expand = c(0, 
            0)) + scale_y_continuous("N", labels = label_number(scale = 0.001, suffix = "k"), breaks = breaks_pretty(3), 
        expand = expansion(mult = c(0, 0.06))) + theme(axis.title.y.left = element_text(angle = 0, vjust = 0.5))
    s <- copy(r$Status)
    s[, `:=`(migrant, factor(mover, c(0, 1), c("Non-migrant", "Migrant")))]
    a1 <- ggplot(s, aes(migrant, n, fill = cohort)) + geom_col(position = position_dodge(0.75), width = 0.65) + 
        geom_text(aes(label = format(n, big.mark = ",", trim = TRUE)), position = position_dodge(0.75), vjust = -0.4, 
            size = 6/.pt, family = font) + scale_fill_manual(values = cols) + th + scale_y_continuous("Participants", 
        labels = label_number(big.mark = ","), expand = expansion(mult = c(0, 0.16))) + labs(x = NULL)
    a2 <- ggplot(r$Counts, aes(count, fraction, colour = cohort)) + geom_line(linewidth = 0.4) + geom_point(size = 1) + 
        scale_colour_manual(values = cols) + th + scale_x_continuous("MT", breaks = pretty(r$Counts$count, 5)) + 
        scale_y_continuous("Participants (%)", labels = label_percent(accuracy = 1))
    a3 <- ggplot(r$Density, aes(x, y, colour = cohort, fill = cohort)) + geom_area(alpha = 0.12, position = "identity", 
        colour = NA) + geom_line(linewidth = 0.45) + cs + th + coord_cartesian(expand = FALSE) + labs(x = "log(1 + MD)", 
        y = "Density")
    forest <- function(o) {
        d <- copy(r$Main[outcome == o])
        d[, `:=`(pos, match(exposure, c("MT", "MD", "Migrant")) + ifelse(cohort == "UKHLS", 0.12, -0.12))]
        d[, `:=`(mark, fifelse(status != "OK", "Warning", fifelse(q < 0.05, "q < 0.05", "q >= 0.05")))]
        p <- ggplot(d, aes(estimate, pos, colour = cohort)) + geom_vline(xintercept = if (o == "SES") 
            0
        else 1, colour = "grey60", linetype = 2, linewidth = 0.18) + geom_segment(aes(x = lo, xend = hi, yend = pos), 
            linewidth = 0.35) + geom_point(aes(shape = mark), size = 1.9, stroke = 0.45) + scale_shape_manual(values = c(`q < 0.05` = 16, 
            `q >= 0.05` = 1, Warning = 4)) + scale_colour_manual(values = cols) + th + scale_y_continuous(NULL, 
            breaks = 1:3, labels = c("MT (per SD)", "MD (per SD)", "Migrant"), limits = c(0.4, 3.85)) + labs(x = if (o == 
            "SES") 
            paste0(intToUtf8(946), " (95% CI)")
        else "HR (95% CI)")
        xr <- range(c(d$lo, d$hi))
        xpos <- function(f) if (o == "SES") 
            xr[1] + diff(xr) * f
        else exp(log(xr[1]) + diff(log(xr)) * f)
        p <- p + annotate("text", x = xpos(0), y = 3.6, label = if (o == "SES") 
            "SES"
        else "All-cause mortality", hjust = 0, family = font, size = 7/.pt)
        if (o == "SES") 
            p <- p + annotate("point", x = xpos(c(0.02, 0.52)), y = 0.6, shape = c(16, 1), size = 1.5) + annotate("text", 
                x = xpos(c(0.07, 0.57)), y = 0.6, label = c("q < 0.05", "q >= 0.05"), hjust = 0, family = font, 
                size = 6/.pt)
        if (o == "Mortality") 
            p <- p + scale_x_log10()
        p
    }
    curve <- function(o, e, individual = FALSE) {
        d <- copy(if (individual) 
            r$Within[exposure == e]
        else r$Curves[outcome == o & exposure == e])
        d[, `:=`(xx, if (e == "MD") 
            log1p(x)
        else x)]
        bk <- if (e == "MD") 
            c(0, 1, 10, 100, 1000, 10000)
        else pretty(d$x, 5)
        bk <- bk[bk >= min(d$x) & bk <= max(d$x)]
        p <- ggplot(d, aes(xx, estimate, colour = cohort, fill = cohort)) + geom_ribbon(aes(ymin = lo, ymax = hi), 
            alpha = 0.12, colour = NA) + geom_line(linewidth = 0.5) + cs + th + scale_x_continuous(if (e == "MD") 
            "MD (km; log scale)"
        else "MT", breaks = if (e == "MD") 
            log1p(bk)
        else bk, labels = as.character(bk), expand = expansion(mult = c(0.055, 0.025))) + labs(y = if (o == "Mortality") 
            "HR"
        else if (individual) 
            "Individual SES change (SD)"
        else "SES difference (SD)")
        p <- p + geom_hline(yintercept = if (o == "Mortality") 
            1
        else 0, linetype = 2, linewidth = 0.18, colour = "grey60")
        if (e == "MT") 
            p <- p + geom_point(size = 0.85)
        test <- unique(d[, .(cohort, p_overall, p_nonlinear)])
        val <- function(v) if (v < 0.001) 
            "< 0.001"
        else paste("=", sprintf("%.3f", v))
        test[, `:=`(label, paste0("paste(italic(p), ' overall ", vapply(p_overall, val, character(1)), "; ', italic(p), ' non-linear ", 
            vapply(p_nonlinear, val, character(1)), "')"))]
        yr <- range(c(d$lo, d$hi, if (o == "Mortality") 1 else 0))
        test[, `:=`(xx = min(d$xx) + 0.01 * diff(range(d$xx)), yy = yr[2] + diff(yr) * ifelse(cohort == "UKHLS", 
            0.24, 0.14))]
        p + geom_text(data = test, aes(xx, yy, label = label, colour = cohort), inherit.aes = FALSE, parse = TRUE, 
            hjust = 0, size = 6/.pt, family = font) + scale_y_continuous(expand = expansion(mult = c(0.04, 0.04)))
    }
    a <- cowplot::plot_grid(a0, a1, a2, a3, nrow = 1, align = "h", axis = "tb")
    b <- cowplot::plot_grid(forest("SES"), forest("Mortality"), curve("SES", "MD"), curve("SES", "MT"), nrow = 1, 
        align = "h", axis = "tb")
    d <- cowplot::plot_grid(curve("Mortality", "MD"), curve("Mortality", "MT"), curve("SES", "MD", TRUE), curve("SES", 
        "MT", TRUE), nrow = 1, align = "h", axis = "tb")
    fig <- cowplot::ggdraw() + cowplot::draw_plot(a, 0.018, 0.668, 0.975, 0.283) + cowplot::draw_plot(b, 0.018, 
        0.356, 0.975, 0.283) + cowplot::draw_plot(d, 0.018, 0.044, 0.975, 0.283)
    for (i in 1:5) fig <- fig + cowplot::draw_label(c("a", "b", "c", "d", "e")[i], x = c(0.005, 0.005, 0.493, 0.005, 
        0.493)[i], y = c(0.956, 0.644, 0.644, 0.332, 0.332)[i], hjust = 0, vjust = 1, size = 10, fontfamily = font, 
        fontface = "bold")
    legend <- data.table(cohort = names(cols), x = c(0, 1))
    key <- ggplot(legend) + geom_rect(aes(xmin = x, xmax = x + 0.12, ymin = 0.31, ymax = 0.69, fill = cohort)) + 
        geom_segment(aes(x = x + 0.18, xend = x + 0.4, y = 0.5, yend = 0.5, colour = cohort), linewidth = 0.5) + 
        geom_point(aes(x = x + 0.29, y = 0.5, colour = cohort), size = 1.2) + geom_text(aes(x = x + 0.46, y = 0.5, 
        label = cohort), hjust = 0, size = 7/.pt, family = font) + cs + xlim(-0.02, 1.95) + ylim(0, 1) + theme_void() + 
        theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))
    fig <- fig + cowplot::draw_plot(key, 0.33, 0.962, 0.34, 0.033)
    stage <- tempfile(fileext = ".pptx")
    on.exit(unlink(stage), add = TRUE)
    write_ppt(stage, list(fig), 24 * cm, 18 * cm)
    if (!file.copy(stage, file.path(out, "Fig4.pptx"), overwrite = TRUE)) 
        stop("Close Fig4.pptx before replacing it")
    invisible(fig)
}

run_fig4(plot_data$result)

finish("Fig4")

