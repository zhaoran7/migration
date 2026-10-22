# 1. Extended Data Figure 7
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

run_ed7 <- function() {
    r <- plot_data$result
    plots <- list()
    ref <- plot_data$main_result$forest
    baseline <- rbind(ref$main[term %in% c("MD", "Migration cross TTWA"), .(outcome, exposure = term, n, event, 
        reference = HR)], ref$tdi[term == paste("Per SD", dtdi_lab), .(outcome, exposure = dtdi_lab, n, event, 
        reference = HR)])
    check <- merge(r$Effects[group == "Baseline"], baseline, by = c("outcome", "exposure"), suffixes = c("", "_main"))
    check[, `:=`(difference, estimate - reference)]
    stopifnot(nrow(check) == 30, all(check$n == check$n_main), all(abs(check$difference) < 1e-07))
    groups <- c("Baseline", "RINT", "2 years", "5 years")
    cols <- c(Baseline = "#344C62", RINT = "#8F65B7", `2 years` = "#7F9FB8", `5 years` = "#D45C69")
    for (e in unname(med_exp)) {
        d <- copy(r$Effects[exposure == e])
        d[, `:=`(outcome = factor(ed_out(outcome), levels = rev(ed_out(forest_out))), group = factor(group, levels = groups))]
        plots[[length(plots) + 1]] <- ggplot(d, aes(estimate, outcome, colour = group, group = group)) + geom_vline(xintercept = 1, 
            colour = "#AAAAAA", linewidth = 0.14) + geom_linerange(aes(xmin = lower, xmax = upper), orientation = "y", 
            position = position_dodge(0.6), linewidth = 0.2) + geom_point(aes(shape = q < 0.05), position = position_dodge(0.6), 
            size = 1.35, stroke = 0.3) + scale_x_log10() + scale_colour_manual(values = cols, labels = c(Baseline = "Baseline", 
            RINT = "MD RINT", `2 years` = "2 years", `5 years` = "5 years")) + scale_shape_manual(values = c(`FALSE` = 1, 
            `TRUE` = 16), guide = "none") + labs(x = "HR", y = NULL, title = e, colour = NULL) + ed_theme() + guides(colour = guide_legend(nrow = 2, 
            byrow = TRUE))
    }
    cross <- rbind(r$Effects[group == "Baseline" & exposure == "Migration cross TTWA"], r$TTWA, fill = TRUE)
    cross[, `:=`(outcome = factor(outcome, levels = forest_out), group = factor(group, levels = c("Baseline", "MD adjusted")), 
        sig = q < 0.05)]
    cross <- cross[, .(outcome, term = as.character(group), HR = estimate, lo = lower, hi = upper, sig)]
    plots[[length(plots) + 1]] <- forest_panel(cross, c(Baseline = "#8F65B7", `MD adjusted` = "#4B2E68")) + guides(shape = "none")
    for (e in unname(med_exp[c("md_z", "dtdi_z")])) {
        d <- copy(r$Curves[exposure == e])
        d[, `:=`(outcome, factor(ed_out(outcome), levels = ed_out(forest_out)))]
        plots[[length(plots) + 1]] <- ggplot(d, aes(value, estimate)) + geom_hline(yintercept = 1, colour = "#AAAAAA", 
            linewidth = 0.14) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = ed_col[e], alpha = 0.18) + 
            geom_line(colour = ed_col[e], linewidth = 0.3) + facet_wrap(~outcome, ncol = 5, scales = "free_y") + 
            labs(x = paste(e, "(SD)"), y = "HR") + ed_theme() + theme(strip.background = element_blank(), strip.text = element_text(size = 6), 
            panel.spacing = grid::unit(2, "mm"))
    }
    boxes <- rbind(ed_boxes(3, 1, 23, 7.8), data.table(x = 0.3, y = c(8.2, 12.4, 17.9), w = 22.4, h = c(3.7, 5, 
        5)))
    ed_save("ED7", plots, boxes, 23, 23.2)
}

run_ed7()

finish("ED7")

