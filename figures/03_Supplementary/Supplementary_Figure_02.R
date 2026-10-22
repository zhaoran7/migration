# 1. Supplementary Figure 2
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

s <- plot_data$supp$base

result <- plot_data$supp

co <- c(UKHLS = "#D45C69", HRS = "#39A6B1")

loc <- c(Birthplace = "#168D83", Recruitment = "#9A4077")

r <- s$ses

vn <- c(edu.years = "Education years", income = "Income", camsis = "CAMSIS")

a <- copy(r$available)[, `:=`(variable, factor(vn[variable], levels = rev(vn)))]

p1 <- ggplot(a, aes(n, variable)) + geom_col(fill = "#6F9EAA", width = 0.62) + geom_text(aes(label = scales::comma(n)), 
    hjust = -0.1, size = pt(6), family = font) + scale_x_continuous(labels = scales::label_number(scale = 0.001, 
    suffix = "k"), expand = expansion(mult = c(0, 0.4))) + labs(x = "Participants", y = NULL) + th()

a <- copy(r$cor)

setnames(a, 1:3, c("v1", "v2", "r"))

a[, `:=`(v1 = factor(vn[as.character(v1)], levels = vn), v2 = factor(vn[as.character(v2)], levels = rev(vn)))]

p2 <- ggplot(a, aes(v1, v2, fill = r)) + geom_tile(colour = "white", linewidth = 0.3) + geom_text(aes(label = sprintf("%.2f", 
    r)), size = pt(6), family = font) + scale_fill_gradient(low = "#F4FAFB", high = "#3F7687", limits = c(0, 1), 
    name = "Pearson r") + coord_equal() + labs(x = NULL, y = NULL) + th() + theme(axis.line = element_blank(), 
    axis.ticks = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1)) + guides(fill = guide_colorbar(barwidth = unit(2.8, 
    "cm"), barheight = unit(0.12, "cm"), ticks = FALSE, title.position = "top"))

a <- copy(r$loading)[, `:=`(variable, factor(vn[variable], levels = rev(vn)))]

p3 <- ggplot(a, aes(loading, variable)) + geom_segment(aes(x = 0, xend = loading, yend = variable), colour = "#6F9EAA", 
    linewidth = 0.35) + geom_point(size = 1.8, colour = "#3F7687") + geom_text(aes(label = sprintf("%.3f", loading)), 
    hjust = -0.3, size = pt(6), family = font) + scale_x_continuous(limits = c(0, 0.8), expand = expansion(mult = c(0, 
    0.04))) + labs(x = "PC1 loading", y = NULL) + th()

p4 <- ggplot(r$variance, aes(component, variance)) + geom_col(fill = "#9DC0C8", width = 0.6) + geom_text(aes(label = scales::percent(variance, 
    accuracy = 0.1)), vjust = -0.5, size = pt(6), family = font) + scale_y_continuous(labels = scales::percent, 
    limits = c(0, 0.65), expand = expansion(mult = c(0, 0))) + labs(x = NULL, y = "Variance explained") + th()

cp <- plot_data$result$ses$comp

left <- tag(list(p1, p2, p3, p4))

right <- tag(list(fig2b_components(cp)), labels = "e")

supp_save("S2", cowplot::plot_grid(left, right, nrow = 1, rel_widths = c(2, 1)), 19, 11.8)

