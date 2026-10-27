# 1. Supplementary Figure 05
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

s <- plot_data$supp$base

result <- plot_data$supp

co <- c(UKHLS = "#D45C69", HRS = "#39A6B1")

loc <- c(Birthplace = "#168D83", Recruitment = "#9A4077")

a <- copy(s$regional$area)

a[, `:=`(location = fifelse(context == "origin", "Birthplace", "Recruitment"), outcome0 = factor(outcome(outcome), 
    levels = outcome(forest_out)), precision = 2 * 1.96 * se)]

p1 <- ggplot(a, aes(n, colour = location)) + stat_ecdf(geom = "step", linewidth = 0.3) + scale_x_log10(labels = scales::label_number(big.mark = ",")) + 
    scale_y_continuous(labels = scales::percent) + cols(loc) + labs(x = "Participants per TTWA", y = "Cumulative share of TTWAs") + 
    th()

q <- a[, .(Areas = .N, adequate = sum(n >= 100 & events >= 10)), by = .(outcome0, location)]

q[, `:=`(fraction, adequate/Areas)]

p2 <- ggplot(q, aes(outcome0, fraction, colour = location, group = location)) + geom_linerange(aes(ymin = 0, ymax = fraction), 
    position = position_dodge(0.5), linewidth = 0.3) + geom_point(position = position_dodge(0.5), size = 1.5) + 
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) + scale_x_discrete(labels = function(x) gsub(" ", 
    "\n", x, fixed = TRUE)) + cols(loc) + labs(x = NULL, y = "Areas meeting support threshold") + th() + theme(axis.text.x = element_text(angle = 45, 
    hjust = 1, size = 6))

p3 <- ggplot(a, aes(events, precision, colour = location)) + geom_point(size = 0.55, alpha = 0.65) + scale_x_continuous(trans = scales::transform_log1p(), 
    breaks = c(0, 100, 10000)) + cols(loc) + facet_wrap(~outcome0, ncol = 5, scales = "free_y") + labs(x = "Events per TTWA", 
    y = "Conditional 95% CI width (log HR)") + th() + theme(legend.position = "none", strip.text = element_text(size = 6))

upper <- tag(list(p1, p2))

supp_save("S5", cowplot::plot_grid(upper, tag(list(p3), labels = "c"), ncol = 1, rel_heights = c(1, 1.35)), 18, 
    16)

