# 1. Figure 5
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

med_pick <- function(x) if ("selected" %in% names(x)) x[selected == TRUE] else x[0]

med_show <- function(x, max_nodes = 40L) {
    x <- med_pick(copy(x))[exposure %in% med_view]
    if (!nrow(x)) 
        return(x)
    rank <- x[, .(prop = max(proportion), p = min(p_indirect_bonf)), by = feature]
    setorder(rank, -prop, p, feature)
    md <- rank[feature %in% x[exposure == "MD", feature], feature]
    ids <- head(c(md, setdiff(rank$feature, md)), max_nodes)
    x[feature %in% ids]
}

med_protein_plot <- function(x) {
    d <- copy(x)[Chromosome %in% c(as.character(1:22), "X") & is.finite(Start)]
    d[, `:=`(chr, match(Chromosome, c(as.character(1:22), "X")))]
    ends <- d[, .(length = max(Start) + 1e+06), by = chr][order(chr)]
    ends[, `:=`(offset, shift(cumsum(length), fill = 0))]
    d <- merge(d, ends[, .(chr, offset)], by = "chr", sort = FALSE)
    d[, `:=`(pos = Start + offset, score = sign(beta) * (-log10(pmax(p, 1e-300))), sig = p_bonf < 0.05)]
    lapply(med_view, function(e) {
        s <- d[exposure == e]
        count <- if (e == dtdi_lab) 
            9L
        else 5L
        labels <- s[sig == TRUE & !is.na(Gene) & nzchar(Gene)][order(p)]
        labels <- unique(labels, by = "Gene")[, head(.SD, count), by = .(direction = sign(beta))]
        threshold <- -log10(0.05/(uniqueN(x$feature) * length(med_exp)))
        ggplot(s, aes(pos, score)) + geom_hline(yintercept = 0, linewidth = 0.1, colour = "black") + geom_hline(yintercept = c(-threshold, 
            threshold), linetype = "dashed", linewidth = 0.1, colour = "#808080") + geom_point(data = s[sig == 
            FALSE], colour = "#CDD2D5", size = 0.45, alpha = 0.65) + geom_point(data = s[sig == TRUE], colour = med_col[[e]], 
            size = 0.8, alpha = 0.85) + ggrepel::geom_text_repel(data = labels, aes(label = Gene), size = pt(6), 
            family = font, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.08, box.padding = 0.12, 
            point.padding = 0.08, seed = 11, max.iter = 15000, show.legend = FALSE) + scale_x_continuous(breaks = ends$offset + 
            ends$length/2, labels = ifelse(ends$chr == 23, "X", ends$chr), expand = expansion(mult = c(0.002, 0.005)), 
            guide = guide_axis(n.dodge = 2)) + scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) + annotate("text", 
            x = min(s$pos) + 0.015 * diff(range(s$pos)), y = max(s$score) + 0.12 * diff(range(s$score)), label = e, 
            hjust = 0, family = font, size = pt(7), colour = med_col[[e]]) + labs(x = "Chromosome", y = expression(paste("Signed -", 
            log[10], "(", italic(p), ")"))) + theme0() + theme(axis.text.x = element_text(size = 6), plot.margin = margin(1, 
            2, 1, 2))
    })
}

med_mri_plot <- function(x) {
    d <- copy(x)[exposure %in% med_view]
    d <- d[feature %in% d[p_bonf < 0.05, unique(feature)]]
    categories <- c("Cortical area", "Cortical thickness", "Cortical volume", "Subcortical volume", "White matter FA", 
        "White matter MD")
    colors <- setNames(c("#CDE8DF", "#DCD5ED", "#F4E5BA", "#D2E5F1", "#F0D8D1", "#DCE9C9"), categories)
    d[, `:=`(category, ifelse(`MRI class` == "cortex", paste("Cortical", Measure), ifelse(`MRI class` == "subcortex", 
        "Subcortical volume", paste("White matter", Measure))))]
    rows <- unique(d[, .(feature, category, Region, Hemisphere)])
    rows[, `:=`(order0, match(category, categories))]
    setorder(rows, order0, Region, Hemisphere)
    nr <- nrow(rows)
    names0 <- c(caudalmiddlefrontal = "Caudal middle frontal", middletemporal = "Middle temporal", parsopercularis = "Pars opercularis", 
        superiortemporal = "Superior temporal", transversetemporal = "Transverse temporal", GlobalMeanThickness = "Global mean thickness", 
        bankssts = "Banks of superior temporal sulcus", isthmuscingulate = "Isthmus cingulate", inferiorparietal = "Inferior parietal", 
        TotalSurface = "Total surface area", superiorfrontal = "Superior frontal")
    rows[, `:=`(label, fifelse(Region %in% names(names0), unname(names0[Region]), Region))]
    rows[, `:=`(label, paste0(toupper(substr(label, 1, 1)), substr(label, 2, nchar(label))))]
    rows[, `:=`(y = nr - seq_len(nr) + 1, label = paste0(label, ifelse(Hemisphere %in% c("left", "right"), paste0(" (", 
        toupper(substr(Hemisphere, 1, 1)), ")"), "")))]
    d[, `:=`(xx = match(exposure, med_view) + 3.7, yy = rows$y[match(feature, rows$feature)], score = fifelse(p_bonf < 
        0.05, beta/se, NA_real_))]
    blocks <- rows[, .(ymin = min(y) - 0.5, ymax = max(y) + 0.5), by = category]
    blocks[, `:=`(label, gsub(" ", "\n", category))]
    blocks[grepl("White matter", category), `:=`(label, paste("White matter", sub("White matter ", "", category), 
        sep = "\n"))]
    limit <- ceiling(max(abs(d$score), na.rm = TRUE))
    head <- data.table(x = seq_along(med_view) + 3.7, label = med_view)
    p <- ggplot(d, aes(xx, yy)) + geom_tile(aes(fill = score), width = 1, height = 1, colour = "#BFC3C7", linewidth = 0.08) + 
        geom_rect(data = blocks, aes(xmin = 0.05, xmax = 4.05, ymin = ymin, ymax = ymax), inherit.aes = FALSE, 
            fill = unname(colors[blocks$category]), colour = "#BFC3C7", linewidth = 0.08) + geom_text(data = blocks, 
        aes(x = 2.05, y = (ymin + ymax)/2, label = label), inherit.aes = FALSE, family = font, size = pt(6), lineheight = 0.95) + 
        geom_text(data = rows, aes(x = 6.5, y = y, label = label), inherit.aes = FALSE, hjust = 0, family = font, 
            size = pt(6)) + geom_text(data = head, aes(x = x, y = nr + 1.65, label = label), inherit.aes = FALSE, 
        family = font, size = pt(6), angle = 90, colour = unname(med_col[head$label])) + scale_fill_gradient2(low = med_col[["MD"]], 
        mid = "white", high = med_col[[dtdi_lab]], midpoint = 0, limits = c(-limit, limit), na.value = "white", 
        name = "Z score", guide = guide_colorbar(barwidth = unit(2.8, "cm"), barheight = unit(0.12, "cm"), ticks = FALSE, 
            title.position = "top", display = "rectangles", nbin = 128)) + coord_fixed(ratio = 1, xlim = c(0, 18.4), 
        ylim = c(0.5, nr + 3), expand = FALSE, clip = "off") + theme_void(base_family = font) + theme(legend.position = "bottom", 
        legend.title = element_text(size = 7), legend.text = element_text(size = 6), legend.margin = margin(0, 
            0, 0, 0), plot.margin = margin(1, 2, 0, 2))
    p
}

med_plot <- function(x, middle) {
    x <- med_show(x)
    if (!nrow(x)) 
        return(ggplot() + annotate("text", 0.5, 0.5, label = "No Bonferroni-significant paths with >5% mediated proportion", 
            family = font, size = pt(7)) + theme_void())
    x[, `:=`(mid, fifelse(!is.na(Gene) & nzchar(Gene), Gene, `Display name`))]
    x[, `:=`(exposure, factor(exposure, levels = med_view))]
    own <- x[, .(mid = first(mid), value = sum(proportion), share = sum(proportion[exposure == med_view[1]])/sum(proportion), 
        target = weighted.mean(match(outcome, forest_out), proportion)), by = feature]
    setorder(own, -share, target, mid)
    ex <- x[, .(value = sum(proportion)), by = exposure][order(exposure)]
    od <- x[, .(value = sum(proportion)), by = outcome][order(match(outcome, forest_out))]
    od[, `:=`(label, unname(int_lab[outcome]))]
    gap <- 0.006
    minnode <- 0.015
    scale <- uniroot(function(k) sum(pmax(k * own$value, minnode)) + gap * (nrow(own) - 1) - 0.95, c(0, 1))$root
    place <- function(d, minimum = minnode) {
        d[, `:=`(height, pmax(value * scale, minimum))]
        space <- if (nrow(d) > 1) 
            (0.95 - sum(d$height))/(nrow(d) - 1)
        else 0
        d[, `:=`(y, 0.985 - cumsum(height) + height/2 - (seq_len(.N) - 1) * space)]
        d
    }
    own <- place(own)
    ex <- place(ex)
    space <- (0.95 - sum(ex$height))/(nrow(ex) + 1)
    ex[, `:=`(y, 0.985 - space - cumsum(height) + height/2 - (seq_len(.N) - 1) * space)]
    od <- place(od, 0.016 * lengths(strsplit(od$label, "\n", fixed = TRUE)))
    x[, `:=`(ei = match(exposure, ex$exposure), mi = match(feature, own$feature), oi = match(outcome, od$outcome), 
        color = unname(med_flow[as.character(exposure)]))]
    left <- x[, .(value = sum(proportion), color = first(color)), by = .(ei, mi)]
    port <- function(d, by, node, column) {
        d[, `:=`(c(column), node$y[.BY[[1]]] + value * scale/2 - cumsum(value) * scale + sum(value) * scale/2), 
            by = by]
        d
    }
    setorder(left, ei, mi)
    left <- port(left, "ei", ex, "from")
    setorder(left, mi, ei)
    left <- port(left, "mi", own, "to")
    x[, `:=`(value, proportion)]
    setorder(x, mi, ei, oi)
    x <- port(x, "mi", own, "from")
    setorder(x, oi, ei, mi)
    x <- port(x, "oi", od, "to")
    ribbon <- function(d, xx, prefix) rbindlist(lapply(seq_len(nrow(d)), function(i) {
        half <- d$value[i] * scale/2
        data.table(id = paste0(prefix, i), x = rep(xx, each = 2), y = c(d$from[i] - half, d$from[i] + half, d$to[i] - 
            half, d$to[i] + half), value = d$value[i], color = d$color[i], part = prefix, ei = d$ei[i], mi = d$mi[i])
    }))
    edge <- rbind(ribbon(left, c(0.07, 0.4675), "a"), ribbon(x, c(0.5025, 0.9), "b"))
    edge[, `:=`(id, factor(id, levels = unique(id[order(-value)])))]
    edge[, `:=`(opacity, unname(med_alpha[as.character(ex$exposure[ei])]))]
    key <- data.table(value = c(0.1, 0.2, 0.4), xx = c(0.55, 0.71, 0.87))
    node <- rbindlist(list(ex[, .(xmin = 0.035, xmax = 0.07, y, height, label = as.character(exposure), side = "left", 
        fill = unname(med_node[as.character(exposure)]))], own[, .(xmin = 0.4675, xmax = 0.5025, y, height, label = mid, 
        side = "middle", fill = "#B5B5B5")], od[, .(xmin = 0.9, xmax = 0.935, y, height, label, side = "right", 
        fill = unname(med_out[outcome]))]))
    node[, `:=`(border, grDevices::rgb(t(grDevices::col2rgb(fill)) * 0.49, maxColorValue = 255))]
    ggplot() + ggforce::geom_diagonal_wide(data = edge, aes(x, y, group = id, fill = color, alpha = opacity), n = 80, 
        strength = 0.5, colour = NA, show.legend = FALSE) + geom_rect(data = node, aes(xmin = xmin, xmax = xmax, 
        ymin = y - height/2, ymax = y + height/2, fill = fill, colour = border), alpha = 0.9, linewidth = 0.05, 
        show.legend = FALSE) + geom_text(data = node[side == "middle"], aes(x = xmax + 0.007, y = y, label = label), 
        hjust = 0, family = font, size = pt(6), lineheight = 1) + geom_text(data = node[side == "left"], aes(x = xmin - 
        0.007, y = y, label = label), hjust = 1, family = font, size = pt(6)) + geom_text(data = node[side == "right"], 
        aes(x = xmax + 0.007, y = y, label = label), hjust = 0, family = font, size = pt(6), lineheight = 1) + 
        annotate("text", x = 0.71, y = -0.024, label = "Mediation proportion", family = font, size = pt(6)) + geom_rect(data = key, 
        aes(xmin = xx - 0.035, xmax = xx + 0.035, ymin = -0.043 - value * scale/2, ymax = -0.043 + value * scale/2), 
        fill = "#B5B5B5", colour = NA) + geom_text(data = key, aes(x = xx, y = -0.065, label = scales::percent(value, 
        accuracy = 1)), family = font, size = pt(6)) + scale_fill_identity() + scale_alpha_identity() + scale_colour_identity() + 
        coord_cartesian(xlim = c(-0.04, 1.1), ylim = c(-0.088, 1.01), expand = FALSE, clip = "off") + theme_void(base_family = font) + 
        theme(plot.margin = margin(2, 3, 2, 3))
}

med_figure <- function(r) {
    device <- cowplot::set_null_device("agg")
    on.exit(cowplot::set_null_device(device), add = TRUE)
    a <- med_protein_plot(r$Protein$association)
    b <- med_mri_plot(r$MRI$association)
    c <- med_plot(r$Protein$mediation, "Protein")
    width <- height <- 20.3
    pad <- c(x = 0.33, y = 0.24)
    gap <- 0.3
    inner <- width - 0.6
    top <- 4
    bottom <- height - 0.6 - top - gap
    split <- 6.3
    y <- (pad[["y"]] + bottom + gap)/height
    p <- cowplot::ggdraw() + cowplot::draw_plot(a[[1]], x = pad[["x"]]/width, y = y, width = (inner - gap)/2/width, 
        height = top/height) + cowplot::draw_plot(a[[2]], x = (pad[["x"]] + (inner + gap)/2)/width, y = y, width = (inner - 
        gap)/2/width, height = top/height) + cowplot::draw_plot(b, x = pad[["x"]]/width, y = pad[["y"]]/height, 
        width = split/width, height = bottom/height) + cowplot::draw_plot(c, x = (pad[["x"]] + split + 0.15)/width, 
        y = pad[["y"]]/height, width = (inner - split - 0.15)/width, height = bottom/height) + cowplot::draw_text(c("a", 
        "b", "c"), x = (pad[["x"]] + c(0, 0, split + 0.15))/width, y = (pad[["y"]] + bottom + 0.08 + c(gap + top, 
        0, 0))/height, hjust = 0, vjust = 1, family = font, fontface = "bold", size = 10) + theme(plot.background = element_rect(fill = "white", 
        colour = NA))
    list(plot = p, width = width, height = height)
}

r <- med_figure(plot_data$result)

write_ppt(file.path(out, "Fig5.pptx"), list(r$plot), r$width * cm, r$height * cm)

ragg::agg_png(file.path(out, "Fig5.png"), width = r$width, height = r$height, units = "cm", res = 600)

print(r$plot)

dev.off()

finish("Fig5")

