# 1. Spatial and single-cell result collection
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

enrich_data <- function() {
    base <- file.path(root, "analysis/enrich/pc40")
    coord <- file.path(root, "analysis/enrich")
    sample <- "E16.5_E1S1.MOSTA"
    atlas <- data.table(id = c("brain", "immune", "vascular", "fat"), Atlas = c("Brain", "Immune", "Vascular", 
        "Adipose"), annotation = c("cell_type", "Manually_curated_celltype", "cell_type", "cell_type"))
    tissue <- rbindlist(lapply(c("MD", "DS"), function(k) {
        file <- file.path(base, "gsmap", sample, "cauchy_combination", paste0(sample, "_", k, ".Cauchy.csv.gz"))
        d <- fread(file)
        d[, .(Trait = if (k == "MD") 
            "MD"
        else dtdi_lab, Tissue = annotation, p = p_cauchy, Median_p = p_median, Source = file)]
    }))
    tissue[, `:=`(q = p.adjust(p, "BH"), Bonferroni = p.adjust(p, "bonferroni")), by = Trait]
    tissue[, `:=`(score = -log10(pmax(p, 1e-300)), sig = Bonferroni < 0.05)]
    spatial <- rbindlist(lapply(c("MD", "DS"), function(k) {
        xy <- fread(file.path(coord, "gsmap/official", paste0(k, ".official.csv")))
        setnames(xy, 1, "spot")
        xy[, `:=`(spot, as.character(spot))]
        d <- fread(file.path(base, "gsmap", sample, "spatial_ldsc", paste0(sample, "_", k, ".csv.gz")))
        d[, `:=`(spot, as.character(spot))]
        n <- nrow(d)
        d <- merge(d, xy[, .(spot, x = sx, y = sy, Tissue = annotation)], by = "spot")
        stopifnot(nrow(d) == n, !anyDuplicated(d$spot))
        d[, .(Trait = if (k == "MD") 
            "MD"
        else dtdi_lab, spot, x, y, Tissue, Beta = beta, SE = se, Z = z, p, score = -log10(pmax(p, 1e-300)))]
    }))
    cells <- rbindlist(lapply(seq_len(nrow(atlas)), function(i) {
        rbindlist(lapply(c("MD", "DS"), function(k) {
            file <- file.path(base, "group", atlas$id[i], paste0(k, ".scdrs_group.", atlas$annotation[i]))
            d <- fread(file)
            d[, .(Atlas = atlas$Atlas[i], Trait = if (k == "MD") 
                "MD"
            else dtdi_lab, Cell = group, Cells = n_cell, Controls = n_ctrl, p = assoc_mcp, Z = assoc_mcz, Heterogeneity_p = hetero_mcp, 
                Heterogeneity_Z = hetero_mcz, Cells_FDR_005 = n_fdr_0.05, Source = file)]
        }))
    }))
    cells[, `:=`(q = p.adjust(p, "BH"), Bonferroni = p.adjust(p, "bonferroni")), by = .(Atlas, Trait)]
    cells[, `:=`(score = Z, sig = Bonferroni < 0.05)]
    select <- cells[, .(sig = any(sig), p = min(p), Cells = max(Cells)), by = .(Atlas, Cell)]
    select[Cells >= 20, `:=`(rank, frank(p, ties.method = "first")), by = Atlas]
    select[, `:=`(Display, sig | (!is.na(rank) & rank <= 3))]
    cells <- merge(cells, select[, .(Atlas, Cell, Display)], by = c("Atlas", "Cell"), all.x = TRUE)
    cells[is.na(Display), `:=`(Display, FALSE)]
    umap <- rbindlist(lapply(seq_len(nrow(atlas)), function(i) {
        d <- fread(file.path(coord, "scdrs/coord", paste0(atlas$id[i], ".umap.tsv.gz")))
        ids <- fread(file.path(base, "group", atlas$id[i], "MD.score.gz"), select = 1)[[1]]
        d <- d[cell %in% ids]
        stopifnot(nrow(d) == length(ids))
        d[, .(Atlas = atlas$Atlas[i], Cell = label, x = umap1, y = umap2)]
    }))
    counts <- umap[, .(Cells = .N, x = median(x), y = median(y)), by = .(Atlas, Cell)]
    list(tissue = tissue, spatial = spatial, cells = cells, umap = umap, counts = counts, atlas = atlas, methods = data.table(Method = c("gsMap", 
        "scDRS"), GWAS = rep("migration_distance_pc40 / delta_tdi_ukb_cond, PC1-40 GWAS", 2), Adjustment = rep("MD: age, age2, sex, center, gen_batch, PC1-PC40; Delta TDI: the same plus Birthplace TDI, using the original UKB TDI definition", 
        2), Reference = c("Mouse E16.5 E1S1 MOSTA", "Human brain, immune, vascular and adipose atlases"), Test = c("Tissue-level Cauchy combination", 
        "Cell-type Monte Carlo association, 1000 controls"), Correction = c("Bonferroni across 25 tissues per trait", 
        "Bonferroni across all cell types within each atlas and trait"), Display = c("All tissues; spatial colours capped at -log10(p)=15", 
        "Heatmap: scDRS assoc_mcz Z score, sequential pale-to-deep teal; stars: Bonferroni p<0.05. All significant cell types plus top 3 types per atlas with >=20 cells; UMAP includes scored cells and colours show cell identity"), 
        Genes = c(NA_integer_, 1000L)))
}

if (sys.nframe() == 0L) {
    saveRDS(enrich_data(), file.path(out, "enrich.rds"), compress = "gzip")
}

