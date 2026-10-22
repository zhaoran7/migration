# 1. Supplementary Figure 8
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

supp_protein_map <- function(d, lev) {
    n <- uniqueN(d$feature)
    rows <- rev(levels(d$feature))
    cols <- c("#8F65B7", "#C18543", "#39A6B1", "#D45C69")
    side <- min(0.325, 9.25/n)
    top <- 10.45
    xpos <- 1.82 + (seq_along(lev) - 0.5) * side
    ypos <- top - (seq_len(n) - 0.5) * side
    limit <- max(abs(d$beta), na.rm = TRUE)
    shade <- scales::col_numeric(c("#3F7687", "white", "#D45C69"), c(-limit, limit))
    children <- list()
    add <- function(g) children[[length(children) + 1L]] <<- g
    text <- function(label, x, y, size = 6.5, just = "left", ...) grid::textGrob(label, x = grid::unit(x, "native"), 
        y = grid::unit(y, "native"), just = just, gp = grid::gpar(fontfamily = font, fontsize = size, col = "#222222", 
            ...))
    for (i in seq_along(rows)) {
        add(text(rows[i], 1.72, ypos[i], just = "right"))
        for (j in seq_along(lev)) {
            a <- d[as.character(feature) == rows[i] & as.character(comparison) == lev[j]]
            stopifnot(nrow(a) == 1L)
            add(grid::rectGrob(x = grid::unit(xpos[j], "native"), y = grid::unit(ypos[i], "native"), width = grid::unit(side, 
                "native"), height = grid::unit(side, "native"), gp = grid::gpar(fill = shade(a$beta), col = "#D1D4D6", 
                lwd = 0.4)))
            if (a$bonferroni < 0.05) 
                add(text("*", xpos[j], ypos[i] - 0.015, just = "centre"))
        }
    }
    base <- top - n * side - 0.23
    for (j in seq_along(lev)) {
        add(grid::pointsGrob(x = grid::unit(xpos[j], "native"), y = grid::unit(base, "native"), pch = 16, size = grid::unit(1.8, 
            "mm"), gp = grid::gpar(col = cols[j])))
        yy <- 8.95 - (j - 1) * 1.47
        add(grid::pointsGrob(x = grid::unit(4.05, "native"), y = grid::unit(yy, "native"), pch = 16, size = grid::unit(1.8, 
            "mm"), gp = grid::gpar(col = cols[j])))
        label <- gsub(" vs ", "\nvs ", lev[j], fixed = TRUE)
        add(text(label, 4.3, yy - 0.02, size = 6))
    }
    xx <- seq(4.05, 7.65, length.out = 256)
    add(grid::rectGrob(x = grid::unit(xx, "native"), y = grid::unit(2.3, "native"), width = grid::unit(diff(xx)[1] * 
        1.02, "native"), height = grid::unit(0.16, "native"), gp = grid::gpar(fill = shade(seq(-limit, limit, length.out = 256)), 
        col = NA)))
    ticks <- pretty(c(-limit, limit), n = 4)
    ticks <- ticks[ticks >= -limit & ticks <= limit]
    add(text(expression(italic(beta)), 4.05, 2.72, size = 7))
    for (v in ticks) add(text(format(v, trim = TRUE, nsmall = 2), 4.05 + (v + limit)/(2 * limit) * 3.6, 1.99, just = "centre"))
    add(text("*  Bonferroni significant", 4.05, 1.46, size = 6))
    g <- grid::grobTree(children = do.call(grid::gList, children), vp = grid::viewport(xscale = c(0, 9), yscale = c(0, 
        11), clip = "off"))
    cowplot::ggdraw() + cowplot::draw_grob(g)
}

s <- plot_data$supp$base

result <- plot_data$supp

r <- result

co <- c(UKHLS = "#D45C69", HRS = "#39A6B1")

loc <- c(Birthplace = "#168D83", Recruitment = "#9A4077")

r$pattern <- plot_data$supp$pattern

a <- copy(r$omics$Protein$contrasts)

lev <- supp_pairs()$comparison

a[, `:=`(comparison, factor(comparison, levels = lev))]

sig <- a[bonferroni < 0.05]

counts <- a[, .(Negative = sum(bonferroni < 0.05 & beta < 0), Positive = sum(bonferroni < 0.05 & beta > 0)), by = comparison]

counts <- melt(counts, id.vars = "comparison", variable.name = "Direction", value.name = "n")

counts[, `:=`(comparison, factor(route(comparison), levels = rev(route(lev))))]

p1 <- ggplot(counts, aes(n, comparison, fill = Direction)) + geom_col(width = 0.6) + geom_text(aes(label = ifelse(n > 
    0, n, "")), position = position_stack(vjust = 0.5), size = pt(6), family = font, colour = "white") + fill(c(Negative = "#3F7687", 
    Positive = "#D45C69")) + labs(x = "Significant proteins", y = NULL) + th()

d <- a[Chromosome %in% c(as.character(1:22), "X") & is.finite(Start)]

d[, `:=`(chr, match(Chromosome, c(as.character(1:22), "X")))]

en <- d[, .(length = max(Start) + 1e+06), by = chr][order(chr)]

en[, `:=`(offset, shift(cumsum(length), fill = 0))]

d <- merge(d, en[, .(chr, offset)], by = "chr", sort = FALSE)

d[, `:=`(pos = Start + offset, score = sign(beta) * (-log10(pmax(p, 1e-300))), significant = bonferroni < 0.05)]

label <- d[significant == TRUE][order(p)][, head(.SD, 5), by = comparison]

d[, `:=`(comparison, factor(gsub(" vs ", "\nvs ", route(comparison), fixed = TRUE), levels = gsub(" vs ", "\nvs ", 
    route(lev), fixed = TRUE)))]

label[, `:=`(comparison, factor(gsub(" vs ", "\nvs ", route(comparison), fixed = TRUE), levels = levels(d$comparison)))]

threshold <- -log10(0.05/(4 * 2910))

p2 <- ggplot(d, aes(pos, score)) + geom_hline(yintercept = 0, linewidth = 0.12) + geom_hline(yintercept = c(-threshold, 
    threshold), colour = "#999999", linetype = 2, linewidth = 0.12) + geom_point(data = d[significant == FALSE], 
    colour = "#CFD4D7", size = 0.4, alpha = 0.65) + geom_point(data = d[significant == TRUE], aes(colour = beta > 
    0), size = 0.7) + scale_colour_manual(values = c(`TRUE` = "#D45C69", `FALSE` = "#3F7687"), guide = "none") + 
    ggrepel::geom_text_repel(data = label, aes(label = Gene), family = font, size = pt(6), max.overlaps = Inf, 
        segment.size = 0.1, box.padding = 0.15, seed = 5, min.segment.length = 0) + scale_x_continuous(breaks = en$offset + 
    en$length/2, labels = ifelse(en$chr == 23, "X", en$chr), guide = guide_axis(n.dodge = 2)) + facet_wrap(~comparison, 
    ncol = 2, scales = "free_y") + labs(x = "Chromosome", y = expression(paste("Signed -", log[10], "(", italic(p), 
    ")"))) + th() + theme(legend.position = "none", panel.spacing = unit(4, "mm"))

features <- unique(a[bonferroni < 0.05][order(p)][, head(.SD, 6), by = comparison]$feature)

small <- r$pattern$map[, .(n = .N), by = pattern][n <= 12, pattern]

features <- unique(c(features, r$pattern$map[pattern %in% small, feature]))

mat <- dcast(a[feature %in% features], feature ~ comparison, value.var = "beta")

M <- as.matrix(mat[, -1])

rownames(M) <- mat$feature

features <- rownames(M)[hclust(dist(M))$order]

b <- a[feature %in% features]

b[, `:=`(feature, factor(feature, levels = rev(features)))]

b[, `:=`(comparison, factor(route(comparison), levels = route(lev)))]

p3 <- supp_protein_map(b, route(lev))

b <- copy(r$pattern$profiles)

b[, `:=`(group, factor(route(group), levels = route(tdi_lev)))]

pp <- lapply(sort(unique(b$pattern)), function(k) {
    d <- b[pattern == k]
    m <- d[, .(shape = mean(shape)), by = group]
    ggplot(d, aes(group, shape, group = feature)) + geom_line(colour = if (k == 1) 
        "#D45C69"
    else "#3F7687", alpha = 0.12, linewidth = 0.18) + geom_line(data = m, aes(group = 1), linewidth = 0.5, colour = "#222222") + 
        geom_point(data = m, aes(group = 1), size = 1.2, colour = "#222222") + labs(x = NULL, y = "Scaled group mean", 
        title = paste0("Pattern ", k, " (", uniqueN(d$feature), " proteins)")) + th() + theme(axis.text.x = element_text(angle = 35, 
        hjust = 1), plot.margin = margin(4, 4, 6, 4))
})

p4 <- cowplot::plot_grid(plotlist = pp, ncol = 1)

top <- tag(list(p1), labels = "a")

middle <- tag(list(p2), labels = "b")

bottom <- cowplot::plot_grid(p4, p3, nrow = 1, rel_widths = c(1, 1), labels = c("c", "d"), label_size = 10, label_fontfamily = font, 
    label_fontface = "bold", hjust = 0, vjust = 1, label_x = 0.008, label_y = 0.995)

supp_save("S8", cowplot::plot_grid(top, middle, bottom, ncol = 1, rel_heights = c(3.1, 8.3, 11)), 18, 22.4)

