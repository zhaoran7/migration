# 1. Supplementary Figure 3
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

r <- plot_data$supp$base$geo

a <- r$nearest

p1 <- ggplot(a, aes(X, Y)) + geom_step(colour = "#7F9FB8", linewidth = 0.35) + scale_x_continuous(breaks = log1p(c(0, 
    1, 10, 100)), labels = c("0", "1", "10", "100")) + scale_y_continuous(labels = scales::percent, limits = c(0, 
    1), expand = expansion(mult = c(0, 0.01))) + labs(x = "Nearest TDI reference (km)", y = "Cumulative proportion") + 
    theme_supp()

p2 <- ggplot(r$tdi) + geom_rect(aes(xmin = Xmin, xmax = Xmax, ymin = Ymin, ymax = Ymax, fill = MT)) + geom_vline(xintercept = r$cut, 
    colour = "#777777", linewidth = 0.18, linetype = 2) + geom_hline(yintercept = r$cut, colour = "#777777", linewidth = 0.18, 
    linetype = 2) + scale_fill_gradient(low = "#E8F1F6", high = "#2166AC", trans = "log10", name = "Participants") + 
    coord_fixed() + labs(x = "Birthplace TDI", y = "Recruitment TDI") + theme_supp() + guides(fill = guide_colorbar(barwidth = unit(2.8, 
    "cm"), barheight = unit(0.12, "cm"), ticks = FALSE, title.position = "top"))

a <- copy(r$cor)

setnames(a, 1:3, c("v1", "v2", "rho"))

vn <- c(md_z = "MD", dtdi_z = dtdi_lab, density_z = "Density change", cross_ttwa = "Cross TTWA")

a[, `:=`(v1 = factor(vn[as.character(v1)], levels = vn), v2 = factor(vn[as.character(v2)], levels = rev(vn)))]

p3 <- ggplot(a, aes(v1, v2, fill = rho)) + geom_tile(colour = "white", linewidth = 0.3) + geom_text(aes(label = sprintf("%.2f", 
    rho)), size = pt(6), family = font) + scale_fill_gradient2(low = "#2667AD", mid = "white", high = "#D45C69", 
    limits = c(-1, 1), name = "Spearman ρ") + coord_equal() + labs(x = NULL, y = NULL) + theme_supp() + theme(axis.line = element_blank(), 
    axis.ticks = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1)) + guides(fill = guide_colorbar(barwidth = unit(2.8, 
    "cm"), barheight = unit(0.12, "cm"), ticks = FALSE, title.position = "top"))

a <- copy(r$routes)

a[, `:=`(percent, n/sum(n)), by = tdi4]

a[, `:=`(tdi4 = factor(route(tdi4), levels = rev(route(tdi_lev))), urban_route = factor(route(urban_route), levels = route(c("Urban to Urban", 
    "Urban to Rural", "Rural to Urban", "Rural to Rural"))))]

p4 <- ggplot(a, aes(urban_route, tdi4, fill = percent)) + geom_tile(colour = "white", linewidth = 0.3) + geom_text(aes(label = scales::percent(percent, 
    accuracy = 0.1)), size = pt(6), family = font) + scale_fill_gradient(low = "#FBF4EA", high = "#C4834C", name = "Within TDI group") + 
    coord_equal() + labs(x = NULL, y = NULL) + theme_supp() + theme(axis.line = element_blank(), axis.ticks = element_blank(), 
    axis.text.x = element_text(angle = 35, hjust = 1)) + guides(fill = guide_colorbar(barwidth = unit(2.8, "cm"), 
    barheight = unit(0.12, "cm"), ticks = FALSE, title.position = "top"))

a <- copy(plot_data$supp$matching$effects[outcome != "SES" & term == dtdi_lab])

a[, `:=`(outcome, factor(ed_out(outcome), levels = rev(ed_out(forest_out))))]

col <- c(`All participants` = "#D45C69", `Match <= 1 km` = "#3F7687")

p5 <- ggplot(a, aes(estimate, outcome, colour = sample, group = sample)) + geom_vline(xintercept = 1, colour = "#AAAAAA", 
    linewidth = 0.14) + geom_linerange(aes(xmin = lower, xmax = upper), orientation = "y", position = position_dodge(0.55), 
    linewidth = 0.3) + geom_point(aes(shape = q < 0.05), position = position_dodge(0.55), size = 1.7, stroke = 0.4) + 
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") + scale_colour_manual(values = col, 
    labels = c("All participants", "Match ≤1 km"), name = NULL) + labs(x = expression("HR per SD " * Delta * 
    "TDI"), y = NULL) + theme_supp()

b <- copy(plot_data$supp$matching$effects[outcome == "SES" & term != "MD"])

lev <- c(dtdi_lab, "Low to high vs High to low", "High to high vs Low to low", "High to low vs High to high", "Low to high vs Low to low")

b[, `:=`(term, factor(route(term), levels = rev(route(lev))))]

p6 <- ggplot(b, aes(estimate, term, colour = sample, group = sample)) + geom_vline(xintercept = 0, colour = "#AAAAAA", 
    linewidth = 0.14) + geom_linerange(aes(xmin = lower, xmax = upper), orientation = "y", position = position_dodge(0.55), 
    linewidth = 0.3) + geom_point(aes(shape = q < 0.05), position = position_dodge(0.55), size = 1.7, stroke = 0.4) + 
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") + scale_colour_manual(values = col, 
    labels = c("All participants", "Match ≤1 km"), name = NULL) + labs(x = expression("SES " * italic(beta)), 
    y = NULL) + theme_supp() + scale_y_discrete(labels = function(v) gsub(" vs ", "\nvs ", v, fixed = TRUE))

supp_save("S3", tag_supp(list(p1, p2, p3, p4, p5, p6), rel_heights = c(1, 1, 1.15)), 20, 24)

