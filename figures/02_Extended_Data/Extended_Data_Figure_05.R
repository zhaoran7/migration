# 1. Extended Data Figure 05
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

ed_prs("ED5", c("dtdi_q", "tdi4"))

finish("ED5")

