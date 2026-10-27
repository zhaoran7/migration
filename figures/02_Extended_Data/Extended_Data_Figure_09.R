# 1. Extended Data Figure 09
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

run_ed9 <- function() {
    device <- cowplot::set_null_device("agg")
    on.exit(cowplot::set_null_device(device), add = TRUE)
    r <- plot_data$result
    cols <- c(UKHLS = "#B76585", HRS = "#39A6B1")
    th <- ed_theme() + theme(legend.position = "none", plot.margin = margin(3, 4, 3, 4))
    n <- copy(r$Effects[group %in% c("Core", "Adjusted")])
    n[, `:=`(exposure = factor(exposure, c("Migrant", "MD", "MT")), group = factor(group, c("Core", "Adjusted")), 
        cohort = factor(cohort, names(cols)), outcome = factor(outcome, c("SES", "Mortality"), c("SES", "All-cause mortality")))]
    n[, `:=`(colour, factor(paste(cohort, group), c("UKHLS Core", "UKHLS Adjusted", "HRS Core", "HRS Adjusted")))]
    a <- ggplot(n, aes(exposure, n, fill = colour)) + geom_col(position = position_dodge(0.7), width = 0.6) + facet_grid(outcome ~ 
        cohort, scales = "free_y") + scale_fill_manual(values = c(`UKHLS Core` = "#ECD8E2", `UKHLS Adjusted` = "#B76585", 
        `HRS Core` = "#D5E4E8", `HRS Adjusted` = "#39A6B1")) + scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), 
        expand = expansion(mult = c(0, 0.1))) + labs(x = NULL, y = "N", fill = NULL) + ed_theme() + theme(axis.title.y = element_text(angle = 0, 
        vjust = 0.5, margin = margin(r = 4)), strip.background = element_blank(), strip.text = element_text(size = 6), 
        legend.position = "bottom") + guides(fill = guide_legend(nrow = 2, byrow = TRUE))
    forest <- function(o) {
        d <- copy(r$Effects[group == "Adjusted" & outcome == o])
        d[, `:=`(cohort, factor(cohort, names(cols)))]
        d[, `:=`(pos, match(exposure, c("MT", "MD", "Migrant")) + ifelse(cohort == "UKHLS", 0.12, -0.12))]
        p <- ggplot(d, aes(estimate, pos, colour = cohort)) + geom_vline(xintercept = if (o == "SES") 
            0
        else 1, colour = "grey60", linetype = 2, linewidth = 0.18) + geom_linerange(aes(xmin = lower, xmax = upper), 
            orientation = "y", linewidth = 0.3) + geom_point(aes(shape = p < 0.05), size = 1.9, stroke = 0.4) + 
            scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = expression(italic(p) >= 0.05, italic(p) < 
                0.05)) + scale_colour_manual(values = cols) + scale_y_continuous(NULL, breaks = 1:3, labels = c("MT (per SD)", 
            "MD (per SD)", "Migrant"), limits = c(0.4, 3.5)) + labs(x = if (o == "SES") 
            expression(beta)
        else "HR", title = if (o == "SES") 
            "SES"
        else "All-cause mortality", colour = NULL, shape = NULL) + th + theme(legend.position = "bottom", legend.box = "vertical") + 
            guides(colour = guide_legend(nrow = 1, order = 1), shape = if (o == "SES") 
                guide_legend(nrow = 1, order = 2, override.aes = list(colour = "black"))
            else "none")
        if (o == "Mortality") 
            p <- p + scale_x_log10()
        p
    }
    curve <- function(o, e, individual = FALSE) {
        d <- copy(if (individual) 
            r$Within[exposure == e]
        else r$Curves[outcome == o & exposure == e])
        d[, `:=`(cohort, factor(cohort, names(cols)))]
        d[, `:=`(xx, if (e == "MD") 
            log1p(x)
        else x)]
        bk <- if (e == "MD") 
            c(0, 1, 10, 100, 1000, 10000)
        else pretty(d$x, 5)
        bk <- bk[bk >= min(d$x) & bk <= max(d$x)]
        p <- ggplot(d, aes(xx, estimate, colour = cohort, fill = cohort)) + geom_hline(yintercept = if (o == "Mortality") 
            1
        else 0, colour = "grey60", linetype = 2, linewidth = 0.18) + geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, 
            colour = NA) + geom_line(linewidth = 0.4) + scale_colour_manual(values = cols) + scale_fill_manual(values = cols) + 
            scale_x_continuous(if (e == "MD") 
                "MD (km; log scale)"
            else "MT", breaks = if (e == "MD") 
                log1p(bk)
            else bk, labels = as.character(bk)) + labs(y = if (o == "Mortality") 
            "HR"
        else if (individual) 
            "Individual SES change (SD)"
        else "SES difference (SD)") + th
        if (e == "MT") 
            p <- p + geom_point(size = 0.9)
        test <- unique(d[, .(cohort, p_overall, p_nonlinear)])
        test[, `:=`(label, paste0("paste(italic(p)[overall], \" \",", vapply(p_overall, function(v) if (v < 0.001) 
            "\"< .001\""
        else paste0("\"= ", sprintf("%.3f", v), "\""), character(1)), ", \"; \", italic(p)[nonlinear], \" \",", 
            vapply(p_nonlinear, function(v) if (v < 0.001) 
                "\"< .001\""
            else paste0("\"= ", sprintf("%.3f", v), "\""), character(1)), ")"))]
        yr <- range(c(d$lo, d$hi, if (o == "Mortality") 1 else 0), na.rm = TRUE)
        test[, `:=`(xx = min(d$xx) + 0.015 * diff(range(d$xx)), yy = yr[2] + diff(yr) * ifelse(cohort == "UKHLS", 
            0.27, 0.14))]
        p + geom_text(data = test, aes(xx, yy, label = label, colour = cohort), inherit.aes = FALSE, parse = TRUE, 
            hjust = 0, family = font, size = pt(6)) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.06)))
    }
    b <- cowplot::plot_grid(forest("SES"), forest("Mortality"), nrow = 1, align = "h")
    c <- cowplot::plot_grid(curve("SES", "MD"), curve("SES", "MT"), nrow = 1, align = "h")
    d <- cowplot::plot_grid(curve("Mortality", "MD"), curve("Mortality", "MT"), nrow = 1, align = "h")
    e <- cowplot::plot_grid(curve("SES", "MD", TRUE), curve("SES", "MT", TRUE), nrow = 1, align = "h")
    boxes <- rbind(data.table(x = 0.3, y = 0.3, w = 23.4, h = 4.1), data.table(x = rep(c(0.3, 12.2), 2), y = rep(c(4.8, 
        10.9), each = 2), w = 11.5, h = 5.7))
    ed_save("ED9", list(a, b, c, d, e), boxes, 24, 16.9)
}

run_ed9()

finish("ED9")

