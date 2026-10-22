# 1. Fixed variables and model definitions
font <- "Arial"

cm <- 1/2.54

end <- as.IDate("2023-04-01")

dtdi_lab <- paste0(intToUtf8(916), "TDI")

tdi_lev <- c("Low to low", "High to low", "Low to high", "High to high")

tdi_pair <- data.table(term = c("High to low vs High to high", "High to low vs Low to low", "Low to high vs Low to low"), 
    a = c("High to low", "High to low", "Low to high"), b = c("High to high", "Low to low", "Low to low"))

trait <- c(md = "MD", dtdi = dtdi_lab)

expv <- c(md = "md_z", dtdi = "dtdi_z")

forest_out <- c("All-cause mortality", "CVD mortality", "CAD/MI", "Stroke", "Heart failure", "Atrial fibrillation", 
    "Type 2 diabetes", "Depression", "Anxiety disorder", "Dementia")

cov0 <- c("age", "sex", "center", "gen_batch", paste0("PC", 1:5), "btdi_z")

cov1 <- c(cov0, "bmi", "smoke", "alcohol")

ukb_var <- c(md = "md_z", dtdi = "dtdi_z", ttwa = "cross_ttwa")

ukb_lab <- c(md = "MD", dtdi = dtdi_lab, ttwa = "Migration cross TTWA")

pal <- list(md = c("#D8F0EA", "#ACDFD2", "#01A6AC", "#015661"), dtdi = c("#DCEAF6", "#AFCFEA", "#6FA1D1", 
    "#2667AD"), main = c(md = "#C18543", dtdi = "#B2182B"), tdi = setNames(c("#00A087", "#5F99D2", "#3C5488", 
    "#FF6663"), tdi_lev), rb = c("#2166AC", "#F7F7F7", "#B2182B"), density = "#7F9FB8", hex = c("#D5E4E8", 
    "#9DC0C8", "#6F9EAA", "#3F7687"))

fig2_cov0 <- setdiff(cov0, "btdi_z")

fig2_cov1 <- setdiff(cov1, "btdi_z")

tdi_cov0 <- fig2_cov0

tdi_cov1 <- fig2_cov1

fig2_model <- "joint"

region_cov <- c("md_z", "dtdi_z", setdiff(fig2_cov1, "btdi_z"))

urban_levels <- c("Urban to urban", "Urban to rural", "Rural to urban", "Rural to rural")

med_exp <- c(md_z = "MD", dtdi_z = dtdi_lab, cross_ttwa = "Migration cross TTWA")

med_view <- unname(med_exp[c("md_z", "dtdi_z")])

med_col <- c(MD = "#39A6B1", setNames("#D45C69", dtdi_lab), `Migration cross TTWA` = "#2667AD")

med_node <- med_col[med_view]

med_flow <- setNames(c("#98D3D9", "#F0C9CD"), med_view)

med_alpha <- setNames(c(0.7, 1), med_view)

med_out <- setNames(rep("#7F9FB8", length(forest_out)), forest_out)

med_model <- paste0(fig2_model, "_omics_z_scan_v4_qc")

enrich_pal <- c("#3B44A4", "#4F8DC0", "#BFE5ED", "#FCE7A4", "#FDA65A", "#E43D30")

int_lab <- setNames(c("All-cause\nmortality", "CVD\nmortality", "CAD", "Stroke", "HF", "AF", "T2D", "Depression", 
    "Anxiety", "Dementia", "SES\n(β)"), c(forest_out, "SES"))

region_version <- "20260922.current-model.1"

region_module_version <- "20260922.compact-grid.1"

ed_col <- c(MD = "#39A6B1", setNames("#D45C69", dtdi_lab), `Migration cross TTWA` = "#8F65B7", UKHLS = "#39A6B1", 
    HRS = "#D45C69")

