# 1. Extended Data Figure 08
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

run_ed8 <- function() {
    device <- cowplot::set_null_device("agg")
    on.exit(cowplot::set_null_device(device), add = TRUE)
    r <- plot_data$result
    r$Bounds <- plot_data$bounds
    make <- function(d, cols, show_x, title) {
        d <- copy(d)
        d[, `:=`(term = group, sig = is.finite(FDR) & FDR < 0.05)]
        forest_panel(d, cols, show_x) + labs(y = "HR per SD MD", title = title) + theme(plot.title = element_text(size = 7), 
            legend.background = element_blank(), legend.key = element_blank(), legend.position = "bottom") + guides(shape = "none")
    }
    a <- cowplot::plot_grid(make(r$Strata[strat == "tdi"], pal$tdi, FALSE, "MD by TDI group") + theme(legend.position = "none"), 
        ed_tdi_key(), ncol = 1, rel_heights = c(3.05, 0.65))
    sex <- copy(r$Extra[strat == "sex"])
    sex[, `:=`(group, fifelse(group == "Female", "Women", "Men"))]
    b <- make(sex, c(Women = "#D45C69", Men = "#2667AD"), FALSE, "MD by sex")
    c <- make(r$Extra[strat == "age"], setNames(c("#EFA88A", "#D8815E", "#B95837", "#8F3026"), c("40-49", "50-59", 
        "60-69", "70+")), TRUE, "MD by age")
    main <- plot_data$selected[, .(interaction, outcome)]
    selected <- r$Tests[FDR < 0.05][!main, on = .(interaction, outcome)]
    selected <- selected[order(match(interaction, c("MD x density change", "MD x age", "MD x sex", "dTDI x age", 
        "dTDI x sex")), match(outcome, forest_out))]
    curves <- lapply(seq_len(nrow(selected)), function(i) {
        s <- selected[i]
        type <- if (s$interaction == "MD x density change") 
            "density"
        else if (grepl("age$", s$interaction)) 
            "age"
        else "sex"
        int_curve(r$Effects, r$Tests, s$interaction, s$outcome, ed_out(s$outcome), type, r$Bounds)
    })
    d <- cowplot::plot_grid(plotlist = curves, ncol = 4, align = "hv")
    dh <- ceiling(length(curves)/4) * 4.2
    boxes <- rbind(data.table(x = 0.3, y = c(0.3, 4.3, 8.3), w = 22.4, h = 3.7), data.table(x = 0.3, y = 12.4, 
        w = 22.4, h = dh))
    ed_save("ED8", list(a, b, c, d), boxes, 23, 12.7 + dh)
}

run_ed8()

finish("ED8")

