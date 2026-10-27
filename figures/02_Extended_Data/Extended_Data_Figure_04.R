# 1. Extended Data Figure 04
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

ed_prs("ED4", c("md_q", "ttwa2"))

finish("ED4")

