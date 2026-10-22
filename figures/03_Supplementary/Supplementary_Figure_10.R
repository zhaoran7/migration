# 1. Supplementary Figure 10
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

supp_spatial <- function() {
    cowplot::set_null_device(function(width, height) ragg::agg_capture(width = width, height = height, units = "in", 
        res = 96))
    atlases <- c(cs8 = "Human CS8 embryo", dlpfc = "Human adult DLPFC")
    rows <- list()
    points <- list()
    plots <- list()
    for (atlas in names(atlases)) {
        xy <- copy(plot_data$coordinates[[atlas]])
        a <- rbindlist(lapply(c("MD", "DS"), function(k) {
            d <- copy(plot_data$tissues[[atlas]][[k]])
            d[, `:=`(Trait, if (k == "MD") 
                "MD"
            else dtdi_lab)]
            d
        }))
        a[, `:=`(Atlas = atlases[[atlas]], Bonferroni = p.adjust(p_cauchy, "bonferroni")), by = Trait]
        order <- a[, .(p = min(p_cauchy)), by = annotation][order(p)]$annotation
        a[, `:=`(xx = match(Trait, c("MD", dtdi_lab)), yy = length(order) - match(annotation, order) + 1)]
        p1 <- ggplot(a, aes(xx, yy)) + geom_tile(aes(fill = pmin(15, -log10(p_cauchy))), colour = "#BFC3C7", linewidth = 0.15, 
            width = 0.96, height = 0.96) + geom_text(aes(label = ifelse(Bonferroni < 0.05, "*", "")), size = pt(6), 
            family = font) + scale_x_continuous(breaks = 1:2, labels = c("MD", dtdi_lab), position = "top") + scale_y_continuous(breaks = seq_along(order), 
            labels = rev(order)) + scale_fill_gradientn(colours = enrich_pal, limits = c(0, 15), guide = "none") + 
            coord_fixed(expand = FALSE, clip = "off") + theme_void(base_family = font) + theme(axis.text = element_text(size = 6, 
            colour = "black"), axis.text.y = element_text(hjust = 1, margin = margin(r = 4)), plot.margin = margin(12, 
            8, 5, 5))
        maps <- lapply(c("MD", "DS"), function(k) {
            d <- copy(plot_data$spatial[[atlas]][[k]])
            d <- merge(xy, d, by = "spot", sort = FALSE)
            d[, `:=`(Atlas = atlases[[atlas]], Trait = if (k == "MD") 
                "MD"
            else dtdi_lab)]
            points[[length(points) + 1]] <<- d
            ggplot(d, aes(x, -y, colour = pmin(15, -log10(pmax(p, 1e-300))))) + geom_point(size = if (atlas == 
                "cs8") 
                1.35
            else 1, stroke = 0, alpha = 1) + scale_colour_gradientn(colours = enrich_pal, limits = c(0, 15), guide = "none") + 
                coord_equal(expand = FALSE) + labs(title = if (k == "MD") 
                "MD"
            else dtdi_lab) + theme_void(base_family = font) + theme(plot.title = element_text(size = 7, hjust = 0.5), 
                plot.margin = margin(4, 4, 4, 4))
        })
        rows[[atlas]] <- a
        plots[[atlas]] <- cowplot::plot_grid(p1, maps[[1]], maps[[2]], nrow = 1, rel_widths = c(0.36, 0.32, 0.32), 
            labels = c(if (atlas == "cs8") "a" else "b", "", ""), label_size = 10, label_fontface = "bold", label_fontfamily = font)
    }
    legend <- ggplot(data.table(x = seq(0, 15, length.out = 500)), aes(x, 1, fill = x)) + geom_tile() + scale_fill_gradientn(colours = enrich_pal, 
        limits = c(0, 15), name = expression(paste("-", log[10], "(", italic(p), ")"))) + guides(fill = guide_colorbar(direction = "horizontal", 
        title.position = "top", barwidth = unit(5, "cm"), barheight = unit(0.18, "cm"), ticks = FALSE, display = "rectangles", 
        nbin = 256)) + theme_void(base_family = font) + theme(legend.position = "bottom", legend.text = element_text(size = 6), 
        legend.title = element_text(size = 7), legend.ticks = element_blank())
    p <- cowplot::plot_grid(plots$cs8, plots$dlpfc, cowplot::get_legend(legend), ncol = 1, rel_heights = c(1, 1, 
        0.13))
    supp_save("S10", p, 18, 17)
    a <- rbindlist(rows)
    d <- rbindlist(points)
    print(a[, .(Tissues = .N, Significant = sum(Bonferroni < 0.05)), by = .(Atlas, Trait)])
    invisible(list(tissues = a, spots = d))
}

supp_spatial()

