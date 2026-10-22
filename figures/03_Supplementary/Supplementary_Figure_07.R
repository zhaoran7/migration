# 1. Supplementary Figure 7
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

cols <- c(UKHLS = "#D45C69", HRS = "#39A6B1")

p1 <- ggplot(plot_data$a, aes(Observed.Intervals, Fraction, colour = Cohort)) + geom_line(linewidth = 0.3) + geom_point(size = 1) + 
    scale_x_continuous(breaks = c(1, 4, 8, 12, 16)) + scale_y_continuous(labels = scales::percent) + scale_colour_manual(values = cols, 
    name = NULL) + labs(x = "Observed distance intervals", y = "Participants") + theme_supp()

p2 <- ggplot(plot_data$b, aes(Recorded.Fraction, Cumulative.Share, colour = Cohort)) + geom_step(linewidth = 0.35) + 
    scale_x_continuous(labels = scales::percent, limits = c(0, 1)) + scale_y_continuous(labels = scales::percent, 
    limits = c(0, 1)) + scale_colour_manual(values = cols, name = NULL) + labs(x = "Intervals with distance recorded", 
    y = "Cumulative proportion") + theme_supp()

density <- function(d, label) ggplot(d, aes(Years, Density, colour = Cohort)) + geom_line(linewidth = 0.35) + scale_colour_manual(values = cols, 
    name = NULL) + labs(x = label, y = "Density") + theme_supp()

supp_save("S7", tag_supp(list(p1, p2, density(plot_data$c, "Observed mortality follow-up (years)"), density(plot_data$d, 
    "Follow-up interval (years)"))), 18, 12)

