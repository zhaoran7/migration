# 1. Supplementary Figure 4
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

supp_birth <- function() {
    r <- plot_data$supp
    d <- copy(r$adjustment)
    colors <- c(Primary = "#39A6B1", `Birthplace TDI adjusted` = "#8F65B7")
    labels <- setNames(c("All-cause death", "CVD death", "CAD", "Stroke", "HF", "AF", "T2D", "Depression", "Anxiety", 
        "Dementia"), forest_out)
    panel <- function(a, ses = FALSE) {
        a[, `:=`(adjustment, factor(adjustment, levels = names(colors)))]
        a[, `:=`(row, factor(if (ses) 
            exposure
        else unname(labels[outcome]), levels = if (ses) 
            rev(c("MD", dtdi_lab))
        else rev(unname(labels))))]
        p <- ggplot(a, aes(estimate, row, colour = adjustment, group = adjustment)) + geom_vline(xintercept = if (ses) 
            0
        else 1, colour = "#999999", linewidth = 0.14) + geom_linerange(aes(xmin = lower, xmax = upper), orientation = "y", 
            position = position_dodge(0.5), linewidth = 0.28) + geom_point(aes(shape = q < 0.05), position = position_dodge(0.5), 
            size = 1.65, stroke = 0.4) + scale_colour_manual(values = colors, name = NULL) + scale_shape_manual(values = c(`FALSE` = 1, 
            `TRUE` = 16), guide = "none") + labs(x = if (ses) 
            expression(italic(beta) ~ "(95% CI)")
        else "HR (95% CI)", y = NULL) + theme0() + theme(legend.position = "bottom", legend.background = element_blank(), 
            legend.key = element_blank(), plot.margin = margin(4, 4, 4, 4))
        if (!ses) 
            p <- p + scale_x_log10()
        p
    }
    a <- panel(copy(d[outcome == "SES"]), TRUE)
    b <- panel(copy(d[outcome != "SES" & exposure == "MD"])) + labs(title = "MD")
    c <- panel(copy(d[outcome != "SES" & exposure == dtdi_lab])) + labs(title = dtdi_lab)
    p <- a/(b | c) + plot_layout(heights = c(1, 3), guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = "bottom", 
        plot.tag = element_text(family = font, size = 10, face = "bold"))
    supp_save("S4", list(p), 19, 15)
    fig2_text(file.path(out, "S4.pptx"))
    invisible(r$adjustment)
}

supp_birth()

