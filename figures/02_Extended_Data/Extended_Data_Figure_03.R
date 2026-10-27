# 1. Extended Data Figure 03
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

ed3_qq <- function(d, label, color) {
    ggplot(d, aes(expected, observed)) + geom_ribbon(aes(ymin = null_lower, ymax = null_upper), fill = "#E4E7E9", 
        colour = NA) + geom_abline(slope = 1, intercept = 0, linewidth = 0.14, colour = "#777777") + geom_point(size = 0.45, 
        colour = color, stroke = 0) + annotate("text", x = 0, y = Inf, label = label, hjust = -0.15, vjust = 1.2, 
        family = font, size = pt(7), colour = color) + scale_x_continuous(limits = c(0, ceiling(max(d$expected))), 
        expand = expansion(mult = 0)) + scale_y_continuous(limits = c(0, ceiling(max(d$observed, d$null_upper))), 
        expand = expansion(mult = 0)) + labs(x = expression(Expected ~ -log[10](italic(p))), y = expression(Observed ~ 
        -log[10](italic(p)))) + theme0() + theme(aspect.ratio = 1, plot.margin = margin(2, 4, 2, 2))
}

ed3_zz <- function(d, id, color) {
    lim <- d$metrics$limit
    label <- paste0("italic(r)[g]==", sprintf("%.3f", d$metrics$rg), "~\"(SE ", sprintf("%.3f", d$metrics$rg_se), 
        ")\"")
    ggplot(d$bins, aes(x, y)) + geom_tile(aes(fill = log10(variants + 1)), width = d$bins$width[1], height = d$bins$width[1]) + 
        scale_fill_gradient(low = "#E4EAED", high = "#687B85", guide = "none") + geom_abline(slope = 1, intercept = 0, 
        linewidth = 0.14, linetype = "dashed", colour = "#777777") + geom_hline(yintercept = 0, linewidth = 0.12, 
        colour = "#D7DADD") + geom_vline(xintercept = 0, linewidth = 0.12, colour = "#D7DADD") + geom_point(data = d$signals, 
        aes(x = z_reference, y = z_current), inherit.aes = FALSE, size = 0.5, stroke = 0, colour = color, alpha = 0.8) + 
        annotate("text", x = -lim, y = lim, label = label, parse = TRUE, hjust = -0.1, vjust = 1.1, family = font, 
            size = pt(6)) + annotate("text", x = -lim, y = lim * 0.81, label = paste0(format(d$metrics$matched, 
        big.mark = ","), " SNPs"), hjust = -0.1, vjust = 1.1, family = font, size = pt(6)) + coord_fixed(xlim = c(-lim, 
        lim), ylim = c(-lim, lim), expand = FALSE) + labs(x = if (id == "md") 
        expression(italic(Z)[MD] ~ "(Furuya 2023)")
    else expression(italic(Z)[TDI] ~ "(Recruitment)"), y = if (id == "md") 
        expression(italic(Z)[MD] ~ "(this study)")
    else expression(italic(Z)[Delta * TDI] ~ "(this study)")) + theme0() + theme(plot.margin = margin(2, 4, 2, 
        2))
}

run_ed3 <- function() {
    r <- ed3_manhattan()
    q <- plot_data$qc
    z <- plot_data$comparison
    plots <- list(ed3_qq(q$qq$md, "MD", "#01A6AC"), ed3_qq(q$qq$dtdi, dtdi_lab, "#B2182B"), ed3_zz(z$md, "md", 
        "#01A6AC"), ed3_zz(z$dtdi, "dtdi", "#B2182B"))
    draw <- function() {
        w <- 20.35
        h <- 17.45
        plot <- function(p, x, y, width, height) {
            grid::pushViewport(grid::viewport(x = x/w, y = 1 - y/h, width = width/w, height = height/h, just = c("left", 
                "top")))
            grid::grid.draw(if (inherits(p, "patchwork")) 
                patchwork::patchworkGrob(p)
            else ggplotGrob(p))
            grid::popViewport()
        }
        text <- function(label, x, y, size = 6, just = "left", color = "#222222", face = "plain") grid::grid.text(label, 
            x = x/w, y = 1 - y/h, just = just, gp = grid::gpar(fontfamily = font, fontsize = size, col = color, 
                fontface = face))
        plot(r$plot, 0.35, 0.45, 19.65, 11)
        for (i in 1:4) plot(plots[[i]], 0.35 + (i - 1) * 5, 12, 4.65, 5.1)
        text("a", 0.1, 0.15, 10, face = "bold")
        for (i in 1:4) text(letters[i + 1], 0.1 + (i - 1) * 5, 11.65, 10, face = "bold")
    }
    file <- file.path(out, "ED3.pptx")
    draw_ppt(file, draw, 20.35 * cm, 17.45 * cm)
    tmp <- tempfile("ed3_")
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    unzip(file, exdir = tmp)
    slide <- file.path(tmp, "ppt/slides/slide1.xml")
    doc <- xml2::read_xml(slide)
    nodes <- xml2::xml_find_all(doc, "//*[local-name()='latin' or local-name()='ea' or local-name()='cs'][@typeface='Symbol']")
    for (node in nodes) xml2::xml_set_attr(node, "typeface", font)
    xml2::write_xml(doc, slide)
    stage <- tempfile(fileext = ".pptx")
    on.exit(unlink(stage), add = TRUE)
    zip::zipr(stage, list.files(tmp, recursive = TRUE, all.files = TRUE, no.. = TRUE), root = tmp, mode = 'mirror')
    stopifnot(file.copy(stage, file, overwrite = TRUE))
    invisible(NULL)
    print(q$metrics[, .(exposure, n, variants, invalid, cojo_signals, lambda_all, intercept)])
    invisible(list(file = file, qc = q$metrics, comparison = z))
}

ed3_manhattan <- function() {
    sig <- 5e-08
    text_pt <- 72.27/25.4
    chr_len <- plot_data$manhattan$chromosomes
    chr_len[, `:=`(CHR, gsub("chr", "", V1))]
    chr_len <- chr_len[CHR %in% as.character(1:22)]
    chr_len[, `:=`(CHR = as.integer(CHR), V2 = as.numeric(V2))]
    data.table::setorder(chr_len, CHR)
    stopifnot(identical(chr_len$CHR, 1:22))
    chr_len[, `:=`(offset = cumsum(data.table::shift(V2, fill = 0)), center = cumsum(data.table::shift(V2, fill = 0)) + 
        V2/2)]
    offset <- setNames(chr_len$offset, chr_len$CHR)
    xmax <- max(chr_len$offset + chr_len$V2)
    place_labels <- function(x, side) {
        data.table::setorder(x, x)
        min_gap <- xmax * 0.0055
        cluster_gap <- xmax * 0.003
        left <- xmax * 0.006
        right <- xmax * 0.985
        label_x <- pmin(pmax(x$x, left), right)
        if (length(label_x) > 1) {
            for (i in 2:length(label_x)) label_x[i] <- max(label_x[i], label_x[i - 1] + min_gap)
            if (tail(label_x, 1) > right) {
                label_x <- label_x - (tail(label_x, 1) - right)
                label_x[1] <- max(label_x[1], left)
                for (i in 2:length(label_x)) label_x[i] <- max(label_x[i], label_x[i - 1] + min_gap)
            }
        }
        x[, `:=`(cluster, cumsum(c(TRUE, diff(CHROM) != 0 | diff(x) > cluster_gap)))]
        base <- x[, .(base_abs = max(y) + 0.46), by = cluster]
        x <- base[x, on = "cluster"]
        x[, `:=`(label_x = label_x, label_y = side * base_abs, line_end = side * (base_abs - 0.12), elbow_y = side * 
            pmax(y + 0.24, base_abs - 0.2), hjust = ifelse(side > 0, 0, 1), vjust = 0.5)]
        x[Gene_Symbol %in% c("CD164", "LCORL"), `:=`(label_x, pmin(label_x + xmax * 0.004, right))]
        x
    }
    plot_trait <- function(d, lab, trait, cols, side) {
        d <- data.table::copy(d)
        d[, `:=`(yplot = side * y, grp = paste0(trait, "_", col))]
        if (nrow(lab)) {
            lab <- data.table::copy(lab)
            lab[, `:=`(ypoint, side * y)]
            lab <- place_labels(lab, side)
            yext <- max(nchar(lab$Gene_Symbol)) * 0.27
            ylim <- if (side > 0) 
                c(3, max(lab$label_y, d$yplot) + yext + 0.6)
            else c(min(lab$label_y, d$yplot) - yext - 0.6, -3)
        }
        else {
            ylim <- if (side > 0) 
                c(3, max(d$yplot) + 0.8)
            else c(min(d$yplot) - 0.8, -3)
        }
        br <- seq(3, ceiling(max(abs(ylim))), by = 3)
        p <- ggplot2::ggplot(d, ggplot2::aes(x, yplot, color = grp)) + ggplot2::geom_point(size = 0.95, alpha = 0.76, 
            stroke = 0) + ggplot2::geom_hline(yintercept = side * -log10(sig), linetype = "dashed", linewidth = 0.08, 
            color = "#777777")
        if (nrow(lab)) 
            p <- p + ggplot2::geom_segment(data = lab, ggplot2::aes(x = x, xend = x, y = ypoint, yend = elbow_y), 
                inherit.aes = FALSE, linewidth = 0.045, color = "#000000") + ggplot2::geom_segment(data = lab, 
                ggplot2::aes(x = x, xend = label_x, y = elbow_y, yend = line_end), inherit.aes = FALSE, linewidth = 0.045, 
                color = "#000000") + ggplot2::geom_text(data = lab, ggplot2::aes(x = label_x, y = label_y, label = Gene_Symbol, 
                hjust = hjust, vjust = vjust), inherit.aes = FALSE, family = font, size = 6/text_pt, angle = 90, 
                color = "#222222")
        p + ggplot2::annotate("text", x = xmax * 0.985, y = side * -log10(sig) + 0.45, label = "5e-08", hjust = 1, 
            size = 6/text_pt, family = font, color = "#666666") + ggplot2::annotate("text", x = xmax * 0.985, y = if (side > 
            0) 
            ylim[2] - 0.8
        else ylim[1] + 0.8, label = trait, hjust = 1, family = font, size = 7/text_pt, color = "#333333") + ggplot2::scale_color_manual(values = cols) + 
            ggplot2::scale_x_continuous(limits = c(0, xmax), expand = ggplot2::expansion(mult = c(0.003, 0.003))) + 
            ggplot2::scale_y_continuous(limits = ylim, breaks = if (side > 0) 
                br[br <= ylim[2]]
            else -br[br <= abs(ylim[1])], labels = function(z) abs(z), expand = ggplot2::expansion(mult = c(0.01, 
                0.01))) + ggplot2::labs(x = NULL, y = expression(-log[10](italic(p)))) + ggplot2::coord_cartesian(clip = "off") + 
            ggplot2::theme_classic(base_family = font, base_size = 6) + ggplot2::theme(legend.position = "none", 
            text = ggplot2::element_text(family = font, face = "plain", color = "#333333"), axis.title = ggplot2::element_text(size = 6), 
            axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size = 6), axis.line.x = ggplot2::element_blank(), 
            axis.ticks.x = ggplot2::element_blank(), axis.line.y = ggplot2::element_line(linewidth = 0.07, color = "#000000"), 
            axis.ticks.y = ggplot2::element_line(linewidth = 0.06, color = "#000000"), plot.margin = ggplot2::margin(if (side > 
                0) 
                2
            else 0.2, 4, if (side > 0) 
                0.2
            else 2, 2))
    }
    plot_chr <- function() {
        ggplot2::ggplot(chr_len, ggplot2::aes(center, 0, label = CHR)) + ggplot2::geom_text(family = font, size = 6/text_pt, 
            color = "#333333") + ggplot2::scale_x_continuous(limits = c(0, xmax), expand = ggplot2::expansion(mult = c(0.003, 
            0.003))) + ggplot2::coord_cartesian(clip = "off") + ggplot2::theme_void(base_family = font) + ggplot2::theme(plot.margin = ggplot2::margin(0, 
            4, 0, 2))
    }
    md <- plot_data$manhattan$points$md
    ds <- plot_data$manhattan$points$dtdi
    md_lab <- plot_data$manhattan$labels$md
    ds_lab <- plot_data$manhattan$labels$dtdi
    p <- patchwork::wrap_plots(list(plot_trait(md, md_lab, "MD", c(MD_1 = "#01A6AC", MD_0 = "#8F6B49"), 1), plot_chr(), 
        plot_trait(ds, ds_lab, "ΔTDI", c(ΔTDI_1 = "#2166AC", ΔTDI_0 = "#B2182B"), -1)), ncol = 1, heights = c(1, 
        0.075, 1))
    list(plot = p, points = list(md = md, dtdi = ds), labels = list(md = md_lab, dtdi = ds_lab))
}

run_ed3()

finish("ED3")

