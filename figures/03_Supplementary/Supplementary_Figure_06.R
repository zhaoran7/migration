# 1. Supplementary Figure 6
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

supp_region <- function() {
    folder <- out
    dir.create(folder, showWarnings = FALSE)
    r <- plot_data$result
    geo <- plot_data$geography
    context <- c(origin = "Birthplace", destination = "Recruitment")
    palette <- list(origin = c("#018571", "#FBF8EF", "#91642E"), destination = c("#B52D55", "#FAF6FA", "#754598"))
    maps <- lapply(names(context), function(v) {
        a <- copy(r$areas)[context == v]
        a[n < 200, `:=`(beta, NA_real_)]
        d <- merge(geo$poly, a, by = "code", all.x = TRUE, sort = FALSE)
        z <- merge(geo$centroids, a, by = "code")[is.finite(beta)]
        z <- z[unique(c(head(order(beta), 3), head(order(-beta), 3)))]
        t <- r$tests[context == v]
        ggplot(d) + geom_polygon(aes(x, y, fill = beta, group = piece, subgroup = ring), rule = "evenodd", colour = NA) + 
            ggrepel::geom_text_repel(data = z, aes(x = x, y = y, label = name), inherit.aes = FALSE, family = font, 
                size = pt(6), box.padding = 0.3, point.padding = 0.1, min.segment.length = 0, segment.colour = "#52616D", 
                segment.size = 0.15, seed = 11, max.overlaps = Inf) + scale_fill_gradientn(colours = palette[[v]], 
            limits = c(-0.5, 0.5), oob = scales::squish, na.value = "#E3E5E7", name = "SES deviation (SD)", breaks = c(-0.5, 
                0, 0.5)) + coord_equal(xlim = c(-10000, 690000), ylim = c(0, 1225000), expand = FALSE) + labs(title = paste0(context[[v]], 
            " | regional SD = ", sprintf("%.3f", t$SD), "; VPC = ", sprintf("%.2f", 100 * t$VPC), "%")) + guides(fill = guide_colourbar(barwidth = unit(3.7, 
            "cm"), barheight = unit(0.18, "cm"), ticks = FALSE, title.position = "top", label.position = "bottom")) + 
            theme_void(base_size = 7, base_family = font) + theme(plot.title = element_text(size = 7, hjust = 0.5, 
            colour = "#25313C"), legend.position = "bottom", legend.title = element_text(size = 7), legend.text = element_text(size = 6), 
            plot.margin = margin(3, 3, 0, 3))
    })
    slopes <- rbindlist(lapply(r$slopes, `[[`, "metrics"))
    slopes[, `:=`(context, factor(context, levels = c("destination", "origin"), labels = c("Recruitment", "Birthplace")))]
    slopes[, `:=`(colour, ifelse(context == "Birthplace", "#267B62", "#A04D85"))]
    curves <- lapply(c("md_z", "dtdi_z"), function(v) {
        ggplot(slopes[exposure == v], aes(y = context, x = mean_beta, colour = colour)) + geom_vline(xintercept = 0, 
            linetype = "dashed", linewidth = 0.25, colour = "#969A9D") + geom_segment(aes(x = beta_low, xend = beta_high, 
            yend = context), linewidth = 1.3, alpha = 0.3) + geom_point(size = 2) + scale_colour_identity() + labs(x = expression("SES " * 
            beta * " per SD"), y = NULL, title = if (v == "md_z") 
            "MD"
        else expression(Delta * "TDI")) + scale_y_discrete(expand = expansion(add = 0.65)) + theme_classic(base_family = font, 
            base_size = 7) + theme(axis.text = element_text(size = 6), axis.line = element_line(linewidth = 0.2), 
            axis.ticks = element_line(linewidth = 0.2), plot.title = element_text(size = 7), panel.grid.major.x = element_line(linewidth = 0.15, 
                colour = "#EEF0F2"), plot.margin = margin(3, 3, 3, 3))
    })
    plots <- c(maps, curves)
    w <- 19
    h <- 16.5
    boxes <- data.table(x = c(0.2, 9.65, 0.2, 9.65), y = c(0.4, 0.4, 12.35, 12.35), w = c(9.1, 9.1, 9.1, 9.1), 
        h = c(11.6, 11.6, 3.9, 3.9))
    tags <- data.table(label = letters[1:4], x = c(0.16, 9.61, 0.16, 9.61), y = c(0.16, 0.16, 12.14, 12.14))
    draw <- function() {
        grid::grid.newpage()
        grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
        for (i in 1:4) {
            b <- boxes[i]
            print(plots[[i]], newpage = FALSE, vp = grid::viewport(x = (b$x + b$w/2)/w, y = 1 - (b$y + b$h/2)/h, 
                width = b$w/w, height = b$h/h))
            grid::grid.text(tags$label[i], x = tags$x[i]/w, y = 1 - tags$y[i]/h, just = c("left", "top"), gp = grid::gpar(fontfamily = font, 
                fontsize = 10, fontface = "bold"))
        }
    }
    base <- file.path(folder, "S6")
    draw_ppt(paste0(base, ".pptx"), draw, w * cm, h * cm)
    tmp <- tempfile("ses_math_")
    dir.create(tmp)
    unzip(paste0(base, ".pptx"), exdir = tmp)
    file <- file.path(tmp, "ppt/slides/slide1.xml")
    d <- xml2::read_xml(file)
    math <- xml2::xml_find_all(d, "//*[@typeface='Symbol']")
    xml2::xml_set_attr(math, "typeface", "Arial")
    beta <- xml2::xml_find_all(d, "//*[local-name()='r'][*[local-name()='t' and text()='β']]/*[local-name()='rPr']")
    xml2::xml_set_attr(beta, "i", "1")
    xml2::write_xml(d, file)
    zip::zipr(paste0(base, ".pptx"), list.files(tmp, recursive = TRUE, all.files = TRUE, no.. = TRUE), root = tmp, 
        mode = "mirror")
    unlink(tmp, recursive = TRUE)
    for (res in c(600, 150)) {
        ragg::agg_png(paste0(base, if (res == 150) 
            "_preview"
        else "", ".png"), width = w, height = h, units = "cm", res = res, background = "white")
        draw()
        dev.off()
    }
    svglite::svglite(paste0(base, ".svg"), width = w * cm, height = h * cm)
    draw()
    dev.off()
    invisible(list(result = r, tags = tags, boxes = boxes))
}

supp_region()

