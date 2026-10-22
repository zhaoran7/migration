# 1. Extended Data Figure 6
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

ed_prs("ED6", c("density_shift_q", "urban4"))

finish("ED6")

