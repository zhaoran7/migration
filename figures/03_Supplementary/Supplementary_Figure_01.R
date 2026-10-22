# 1. Supplementary Figure 1
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

supp_flow_plot <- function(r) {
    d <- rbindlist(lapply(seq_len(nrow(r$links)), function(i) {
        a <- r$links[i]
        data.table(id = a$id, quartile = a$quartile, x = c(a$axis + 0.16, a$axis + 0.16, a$axis + 1 - 0.16, a$axis + 
            1 - 0.16), y = c(a$left - a$height, a$left, a$right, a$right - a$height))
    }))
    node <- copy(r$nodes)
    label <- function(v) {
        v <- gsub(" to ", "\n→ ", v, fixed = TRUE)
        v <- sub("Within TTWA", "Within\nTTWA", v, fixed = TRUE)
        sub("Cross TTWA", "Cross\nTTWA", v, fixed = TRUE)
    }
    node[, `:=`(label, label(group))]
    node[, `:=`(small, (top - bottom) < 0.04)]
    node[, `:=`(label_y, ifelse(small, -0.035, (bottom + top)/2))]
    node[small == TRUE, `:=`(label, gsub("\n", " ", label, fixed = TRUE))]
    title <- unique(node[, .(axis, stage)])
    title[, `:=`(stage, gsub(" ", "\n", stage, fixed = TRUE))]
    ggplot() + ggforce::geom_diagonal_wide(data = d, aes(x, y, group = id, fill = quartile), alpha = 0.2, colour = NA, 
        strength = 0.5) + geom_rect(data = node, aes(xmin = axis - 0.16, xmax = axis + 0.16, ymin = bottom, ymax = top), 
        fill = "white", colour = "#B9BEC2", linewidth = 0.14) + geom_segment(data = node[small == TRUE], aes(x = axis, 
        xend = axis, y = bottom, yend = -0.015), linewidth = 0.12, colour = "black") + geom_text(data = node, aes(axis, 
        label_y, label = label), family = font, size = pt(6), lineheight = 0.95) + geom_text(data = title, aes(axis, 
        1.1, label = stage), family = font, size = pt(7), lineheight = 0.95) + scale_fill_manual(values = setNames(unname(pal$tdi), 
        paste0("Q", 1:4)), guide = "none") + coord_cartesian(xlim = c(0.5, 8.5), ylim = c(-0.07, 1.19), expand = FALSE, 
        clip = "off") + theme_void(base_family = font) + theme(plot.margin = margin(8, 6, 6, 6))
}

supp_save("S1", supp_flow_plot(plot_data$result), 18, 11)

