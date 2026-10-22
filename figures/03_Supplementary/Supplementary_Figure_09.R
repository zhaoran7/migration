# 1. Supplementary Figure 9
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

cells <- plot_data$cells

plots <- list()

for (atlas in c("Heart", "Liver")) {
    xy <- copy(plot_data$umap[[atlas]])
    d <- copy(cells[Atlas == atlas])
    labels <- unique(c(d[order(-Cells), head(unique(Cell), 8)], d[Bonferroni < 0.05, Cell]))
    cols <- setNames(c("#3F7687", "#D45C69", "#8F65B7", "#C49A6C", "#62A58C", "#7F9FB8", "#C79BAD", "#A6A95A", 
        "#7865A7", "#D4A84C", "#67A9CF", "#A88462")[seq_along(labels)], labels)
    xy[, `:=`(type, ifelse(label %in% labels, label, "Other"))]
    plots[[length(plots) + 1]] <- ggplot(xy, aes(x, y, colour = type)) + geom_point(size = 0.16, alpha = 0.7) + 
        coord_equal() + theme_void(base_family = font) + scale_colour_manual(values = c(cols, Other = "#D4D6D8"), 
        name = NULL, labels = function(x) vapply(x, function(v) paste(strwrap(v, 25), collapse = "\n"), "")) + 
        guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1), ncol = 1)) + labs(title = atlas) + 
        theme(plot.title = element_text(size = 7), legend.position = "right", legend.text = element_text(size = 6), 
            legend.key.height = unit(2.7, "mm"), legend.key.width = unit(2, "mm"), plot.margin = margin(4, 4, 4, 
                4))
    d[, `:=`(Cell, factor(Cell, levels = rev(unique(d[order(-Z)]$Cell))))]
    plots[[length(plots) + 1]] <- ggplot(d, aes(Trait, Cell, fill = Z)) + geom_tile(colour = "#C8CDD0", linewidth = 0.15) + 
        geom_text(aes(label = ifelse(Bonferroni < 0.05, "*", "")), size = pt(6), family = font) + scale_fill_gradient(low = "#F2FAFA", 
        high = "#168A96", name = "scDRS Z") + coord_equal() + labs(x = NULL, y = NULL) + theme_supp() + theme(axis.line = element_blank(), 
        axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "bottom") + 
        guides(fill = guide_colorbar(barwidth = unit(2.8, "cm"), barheight = unit(0.12, "cm"), ticks = FALSE, title.position = "top"))
}

supp_save("S9", tag_supp(plots, ncol = 2, rel_widths = c(1.2, 1)), 18, max(14, 0.23 * max(cells[, uniqueN(Cell), 
    by = Atlas]$V1) * 2 + 3))

