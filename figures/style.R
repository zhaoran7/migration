suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(patchwork)
    library(cowplot)
    library(scales)
    library(officer)
    library(rvg)
    library(grid)
    library(ggrepel)
})

args <- commandArgs(TRUE)

script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))

folder <- dirname(script)

id <- tools::file_path_sans_ext(basename(script))
id <- sub("^Extended_Data_Figure_", "ED", id)
id <- sub("^Supplementary_Figure_", "S", id)
id <- sub("^Figure_", "Fig", id)
id <- sub("^([A-Za-z]+)0+([1-9][0-9]*)$", "\\1\\2", id)

input <- if (length(args)) args[1] else file.path(folder, "source", paste0(id, "_source_data.xlsx"))

out <- if (length(args) > 1) args[2] else file.path(getwd(), "output")

dir.create(out, recursive = TRUE, showWarnings = FALSE)

out <- normalizePath(out, mustWork = TRUE)

cowplot::set_null_device("agg")

proj <- file.path(dirname(dirname(R.home())),'share/proj')
if(!nzchar(Sys.getenv('PROJ_DATA'))&&file.exists(file.path(proj,'proj.db')))
    Sys.setenv(PROJ_DATA=proj,PROJ_LIB=proj)

font <- "Arial"

cm <- 1/2.54

dtdi_lab <- paste0(intToUtf8(916), "TDI")

tdi_lev <- c("Low to low", "High to low", "Low to high", "High to high")

forest_out <- c("All-cause mortality", "CVD mortality", "CAD/MI", "Stroke", "Heart failure", "Atrial fibrillation", 
    "Type 2 diabetes", "Depression", "Anxiety disorder", "Dementia")

ukb_lab <- c(md = "MD", dtdi = dtdi_lab, ttwa = "Migration cross TTWA")

pal <- list(md = c("#D8F0EA", "#ACDFD2", "#01A6AC", "#015661"), dtdi = c("#DCEAF6", "#AFCFEA", "#6FA1D1", "#2667AD"), 
    main = c(md = "#C18543", dtdi = "#B2182B"), tdi = setNames(c("#00A087", "#5F99D2", "#3C5488", "#FF6663"), tdi_lev), 
    rb = c("#2166AC", "#F7F7F7", "#B2182B"), density = "#7F9FB8", hex = c("#D5E4E8", "#9DC0C8", "#6F9EAA", "#3F7687"))

urban_levels <- c("Urban to urban", "Urban to rural", "Rural to urban", "Rural to rural")

med_exp <- c(md_z = "MD", dtdi_z = dtdi_lab, cross_ttwa = "Migration cross TTWA")

med_view <- unname(med_exp[c("md_z", "dtdi_z")])

med_col <- c(MD = "#39A6B1", setNames("#D45C69", dtdi_lab), `Migration cross TTWA` = "#2667AD")

med_node <- med_col[med_view]

med_flow <- setNames(c("#98D3D9", "#F0C9CD"), med_view)

med_alpha <- setNames(c(0.7, 1), med_view)

med_out <- setNames(rep("#7F9FB8", length(forest_out)), forest_out)

enrich_pal <- c("#3B44A4", "#4F8DC0", "#BFE5ED", "#FCE7A4", "#FDA65A", "#E43D30")

int_lab <- setNames(c("All-cause\nmortality", "CVD\nmortality", "CAD", "Stroke", "HF", "AF", "T2D", "Depression", 
    "Anxiety", "Dementia", "SES\n(β)"), c(forest_out, "SES"))

ed_col <- c(MD = "#39A6B1", setNames("#D45C69", dtdi_lab), `Migration cross TTWA` = "#8F65B7", UKHLS = "#39A6B1", 
    HRS = "#D45C69")

source_read <- function(file) {
    idx <- data.table::as.data.table(openxlsx::read.xlsx(file, "Contents", check.names = FALSE))
    read <- function(id) {
        a <- idx[ID == id]
        if (a$Kind == "null") 
            return(NULL)
        if (a$Kind == "sf") 
            return(sf::st_read(file.path(dirname(file), a$Sheet), quiet = TRUE))
        meta <- if (is.na(a$Metadata) || !nzchar(a$Metadata)) 
            list()
        else jsonlite::fromJSON(a$Metadata, simplifyVector = FALSE)
        if (a$Kind == "list") {
            kids <- idx[Parent == id, ID]
            x <- lapply(kids, read)
            if (!is.null(meta$names)) 
                names(x) <- unlist(meta$names)
            if (!is.na(a$Class) && nzchar(a$Class)) 
                class(x) <- strsplit(a$Class, "|", fixed = TRUE)[[1]]
            return(x)
        }
        x <- openxlsx::read.xlsx(file, sheet = a$Sheet, check.names = FALSE, skipEmptyRows = FALSE, skipEmptyCols = FALSE, 
            na.strings = "<NA>")
        if (nrow(x) < meta$nrow) 
            x <- x[seq_len(meta$nrow), , drop = FALSE]
        names(x) <- unlist(meta$columns)
        for (j in seq_along(x)) {
            info <- meta$types[[j]]
            cl <- unlist(info$class)
            v <- x[[j]]
            if (is.character(v)) 
                v[v == "<EMPTY>"] <- ""
            if ("factor" %in% cl) 
                v <- factor(v, levels = unlist(info$levels), ordered = "ordered" %in% cl)
            else if (any(c("Date", "IDate") %in% cl)) {
                v <- as.Date(as.numeric(v), origin = "1970-01-01")
                if ("IDate" %in% cl) 
                  v <- data.table::as.IDate(v)
            }
            else if (info$type == "integer") 
                v <- as.integer(v)
            else if (info$type == "double") 
                v <- as.numeric(v)
            else if (info$type == "logical") 
                v <- as.logical(v)
            else if (info$type == "character") 
                v <- as.character(v)
            if (names(x)[j] %in% unlist(meta$list_columns)) 
                v <- lapply(v, jsonlite::fromJSON)
            x[[j]] <- v
        }
        if (a$Kind == "vector") {
            v <- x[[1]]
            if (!is.null(meta$names)) 
                names(v) <- unlist(meta$names)
            return(v)
        }
        if (a$Kind == "matrix") {
            x <- as.matrix(x)
            if (!is.null(meta$rownames)) 
                rownames(x) <- unlist(meta$rownames)
            return(x)
        }
        data.table::as.data.table(x)
    }
    read(1L)
}

theme0 <- function(size = 6) {
    theme_classic(base_family = font, base_size = size) + theme(text = element_text(family = font, face = "plain", 
        colour = "#222222"), axis.title = element_text(size = 7), axis.text = element_text(size = 6), axis.line = element_line(linewidth = 0.14), 
        axis.ticks = element_line(linewidth = 0.14), legend.title = element_text(size = 6), legend.text = element_text(size = 6), 
        plot.title = element_text(size = 7, face = "plain"), plot.margin = margin(2, 2, 2, 2))
}

pt <- function(x) x/ggplot2::.pt

ppt_size <- function(file, w, h) {
    tmp <- tempfile("pptx_")
    dir.create(tmp)
    unzip(file, exdir = tmp)
    xml <- file.path(tmp, "ppt/presentation.xml")
    x <- readLines(xml, warn = FALSE)
    w <- format(round(w * 914400), scientific = FALSE, trim = TRUE)
    h <- format(round(h * 914400), scientific = FALSE, trim = TRUE)
    x <- gsub("(<p:sldSz[^>]*cx=\")[0-9]+(\"[^>]*cy=\")[0-9]+(\")", paste0("\\1", w, "\\2", h, "\\3"), x)
    writeLines(x, xml, useBytes = TRUE)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    unlink(file)
    zip::zipr(file, list.files(tmp, recursive = TRUE, all.files = TRUE, no.. = TRUE), root = tmp, mode = "mirror")
}

write_ppt <- function(file, plots, w, h) {
    stage <- tempfile(fileext = ".pptx")
    on.exit(unlink(stage), add = TRUE)
    ppt <- officer::read_pptx()
    for (p in plots) {
        ppt <- officer::add_slide(ppt, layout = "Blank", master = "Office Theme")
        ppt <- officer::ph_with(ppt, rvg::dml(ggobj = p), officer::ph_location(left = 0, top = 0, width = w, height = h))
    }
    print(ppt, target = stage)
    ppt_size(stage, w, h)
    if (grepl("^fig2[a-d]\\.pptx$", basename(file))) 
        fig2_text(stage)
    if (!file.copy(stage, file, overwrite = TRUE)) 
        stop("Close the presentation before replacing: ", file)
    invisible(file)
}

draw_ppt <- function(file, draw, w, h) {
    stage <- tempfile(fileext = ".pptx")
    on.exit(unlink(stage), add = TRUE)
    ppt <- officer::add_slide(officer::read_pptx(), layout = "Blank", master = "Office Theme")
    ppt <- officer::ph_with(ppt, rvg::dml(code = draw()), officer::ph_location(left = 0, top = 0, width = w, height = h))
    print(ppt, target = stage)
    ppt_size(stage, w, h)
    if (grepl("^fig2[a-d]\\.pptx$", basename(file))) 
        fig2_text(stage)
    if (!file.copy(stage, file, overwrite = TRUE)) 
        stop("Close the presentation before replacing: ", file)
    invisible(file)
}

fig2_text <- function(x) {
    if (is.character(x)) {
        tmp <- tempfile("labels_")
        dir.create(tmp)
        on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
        unzip(x, exdir = tmp)
        for (f in list.files(file.path(tmp, "ppt/slides"), "^slide[0-9]+\\.xml$", full.names = TRUE)) {
            doc <- xml2::read_xml(f)
            fig2_text(doc)
            xml2::write_xml(doc, f)
        }
        target <- tempfile(fileext = ".pptx")
        on.exit(unlink(target), add = TRUE)
        zip::zipr(target, list.files(tmp, recursive = TRUE, all.files = TRUE, no.. = TRUE), root = tmp, mode = 'mirror')
        stopifnot(file.copy(target, x, overwrite = TRUE))
        return(invisible(x))
    }
    ns <- c(p = "http://schemas.openxmlformats.org/presentationml/2006/main", a = "http://schemas.openxmlformats.org/drawingml/2006/main")
    for (s in xml2::xml_find_all(x, "//p:sp", ns)) {
        texts <- xml2::xml_find_all(s, ".//a:t", ns)
        old <- paste(xml2::xml_text(texts), collapse = "")
        new <- if (old == "TDI") 
            "Recruitment TDI"
        else old
        new <- gsub("Origin TDI", "Birthplace TDI", new, fixed = TRUE)
        new <- gsub("Destination TDI", "Recruitment TDI", new, fixed = TRUE)
        new <- gsub("Residence TDI", "Recruitment TDI", new, fixed = TRUE)
        if (new != old && length(texts) == 1) {
            tr <- xml2::xml_find_first(s, "p:spPr/a:xfrm", ns)
            off <- xml2::xml_find_first(tr, "a:off", ns)
            ext <- xml2::xml_find_first(tr, "a:ext", ns)
            xy <- as.numeric(xml2::xml_attrs(off)[c("x", "y")])
            w <- as.numeric(xml2::xml_attr(ext, "cx"))
            width <- systemfonts::string_width(c(old, new), family = "Arial", size = 7, res = 720)
            delta <- w * (width[2]/width[1] - 1)
            angle <- as.numeric(xml2::xml_attr(tr, "rot"))/60000 * pi/180
            if (!is.finite(angle)) 
                angle <- 0
            xy <- xy - c(delta/2, 0)
            if (old == "TDI") 
                xy <- xy - delta/2 * c(cos(angle), sin(angle))
            xml2::xml_set_attrs(off, setNames(sprintf("%.0f", xy), c("x", "y")))
            xml2::xml_set_attr(ext, "cx", sprintf("%.0f", w + delta))
            xml2::xml_set_text(texts[[1]], new)
        }
        for (r in xml2::xml_find_all(s, ".//a:r", ns)) {
            text <- xml2::xml_text(xml2::xml_find_first(r, "a:t", ns))
            if (text == intToUtf8(8722)) {
                fonts <- xml2::xml_find_all(r, ".//a:latin | .//a:cs | .//a:ea", ns)
                xml2::xml_set_attr(fonts, "typeface", "Arial")
            }
            if (!grepl("β|(?<![[:alnum:]])p(?![[:alnum:]])", text, perl = TRUE)) 
                next
            pieces <- strsplit(text, "(?<=β)|(?=β)|(?<![[:alnum:]])(?=p\\b)|(?<=\\bp)", perl = TRUE)[[1]]
            for (piece in pieces[nzchar(pieces)]) {
                copy <- xml2::xml_children(xml2::read_xml(paste0("<root xmlns:a=\"", ns["a"], "\">", as.character(r), 
                  "</root>")))[[1]]
                xml2::xml_set_text(xml2::xml_find_first(copy, "a:t", ns), piece)
                if (piece %in% c("p", "β")) 
                  xml2::xml_set_attr(xml2::xml_find_first(copy, "a:rPr", ns), "i", "1")
                if (piece == "β") {
                  fonts <- xml2::xml_find_all(copy, ".//a:latin | .//a:cs | .//a:ea", ns)
                  xml2::xml_set_attr(fonts, "typeface", "Arial")
                }
                xml2::xml_add_sibling(r, copy, .where = "before")
            }
            xml2::xml_remove(r)
        }
    }
    invisible(x)
}

ed_out <- function(x) {
    a <- c(`CAD/MI` = "CAD", `Type 2 diabetes` = "T2D", `Anxiety disorder` = "Anxiety")
    ifelse(x %in% names(a), a[x], x)
}

ed_theme <- function() theme0() + theme(plot.margin = margin(8, 4, 4, 6), legend.position = "bottom", legend.key.height = grid::unit(2.5, 
    "mm"), legend.key.width = grid::unit(3, "mm"), legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 
    0, 0, 0), legend.background = element_blank(), legend.key = element_blank())

ed_boxes <- function(cols, rows, w, h, gap = 0.35) {
    bw <- (w - 0.6 - (cols - 1) * gap)/cols
    bh <- (h - 0.6 - (rows - 1) * gap)/rows
    data.table(x = rep(0.3 + (0:(cols - 1)) * (bw + gap), rows), y = rep(0.3 + (0:(rows - 1)) * (bh + gap), each = cols), 
        w = bw, h = bh)
}

ed_tdi_key <- function() {
    routes <- list(`Low to low` = c(0.17, 0.17), `High to low` = c(0.46, 0.17), `Low to high` = c(0.17, 0.46), 
        `High to high` = c(0.46, 0.46))
    bands <- list()
    nodes <- list()
    t <- seq(0, 1, length.out = 65)
    for (i in seq_along(routes)) {
        x <- 1.05 + (i - 1) * 2.2
        for (g in names(routes)) {
            yy <- routes[[g]][1] + diff(routes[[g]]) * (3 * t^2 - 2 * t^3)
            xx <- x + 1.4 * (3 * (1 - t)^2 * t * 0.45 + 3 * (1 - t) * t^2 * 0.55 + t^3)
            width <- if (g == names(routes)[i]) 
                0.065
            else 0.015
            bands[[length(bands) + 1]] <- data.table(x = c(xx, rev(xx)), y = c(yy - width/2, rev(yy + width/2)), 
                id = paste(i, g), fill = if (g == names(routes)[i]) 
                  pal$tdi[g]
                else "#E5E7EB", order = as.integer(g == names(routes)[i]))
        }
        nodes[[i]] <- CJ(x = c(x, x + 1.4), y = c(0.17, 0.46))
    }
    b <- rbindlist(bands)[order(order)]
    n <- rbindlist(nodes)
    ggplot() + geom_polygon(data = b, aes(x, y, group = id, fill = I(fill)), colour = NA) + geom_rect(data = n, 
        aes(xmin = x - 0.012, xmax = x + 0.012, ymin = y - 0.07, ymax = y + 0.07), fill = "#727B83", colour = NA) + 
        annotate("text", x = 5.05, y = 0.77, label = "Birthplace TDI  →  Recruitment TDI", family = font, size = pt(6)) + 
        annotate("text", x = c(0.86, 0.86, 9.22, 9.22), y = rep(c(0.46, 0.17), 2), label = rep(c("High", "Low"), 
            2), hjust = c(1, 1, 0, 0), family = font, size = pt(6)) + coord_cartesian(xlim = c(0, 10.1), ylim = c(0, 
        0.92), expand = FALSE, clip = "off") + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
}

ed_save <- function(id, plots, boxes, w, h) {
    stopifnot(length(plots) == nrow(boxes), all(boxes$x + boxes$w <= w + 0.001), all(boxes$y + boxes$h <= h + 0.001))
    invisible(lapply(plots, ggplot_build))
    file <- file.path(out, paste0(id, ".pptx"))
    draw <- function() {
        grid::grid.newpage()
        grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
        for (i in seq_along(plots)) {
            b <- boxes[i, ]
            print(plots[[i]], newpage = FALSE, vp = grid::viewport(x = (b$x + b$w/2)/w, y = 1 - (b$y + b$h/2)/h, 
                width = b$w/w, height = b$h/h))
            grid::grid.text(letters[i], x = (b$x + 0.06)/w, y = 1 - (b$y + 0.04)/h, just = c("left", "top"), gp = grid::gpar(fontfamily = font, 
                fontsize = 10, fontface = "bold"))
        }
    }
    draw_ppt(file, draw, w * cm, h * cm)
    tmp <- tempfile("ed_font_")
    dir.create(tmp)
    unzip(file, exdir = tmp)
    for (f in list.files(file.path(tmp, "ppt/slides"), "^slide[0-9]+\\.xml$", full.names = TRUE)) {
        doc <- xml2::read_xml(f)
        nodes <- xml2::xml_find_all(doc, "//*[local-name()='latin' or local-name()='ea' or local-name()='cs'][@typeface='Symbol']")
        xml2::xml_set_attr(nodes, "typeface", "Arial")
        math <- xml2::xml_find_all(doc, "//*[local-name()='r'][*[local-name()='t' and (text()='β' or text()='p' or text()='q')]]/*[local-name()='rPr']")
        xml2::xml_set_attr(math, "i", "1")
        xml2::write_xml(doc, f)
    }
    stage <- tempfile(fileext = ".pptx")
    zip::zipr(stage, list.files(tmp, recursive = TRUE, all.files = TRUE, no.. = TRUE), root = tmp, mode = "mirror")
    stopifnot(file.copy(stage, file, overwrite = TRUE))
    unlink(c(stage, tmp), recursive = TRUE)
    ppt_finish(file)
    invisible(file)
}

supp_save <- function(id, plots, width, height, folder = out) {
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    if (!is.list(plots) || inherits(plots, "ggplot")) 
        plots <- list(plots)
    base <- file.path(folder, id)
    write_ppt(paste0(base, ".pptx"), plots, width * cm, height * cm)
    for (i in seq_along(plots)) {
        stem <- if (length(plots) == 1) 
            base
        else paste0(base, "_", letters[i])
        svglite::svglite(paste0(stem, ".svg"), width = width * cm, height = height * cm, system_fonts = list(sans = font))
        print(plots[[i]])
        dev.off()
        ragg::agg_png(paste0(stem, ".png"), width = width, height = height, units = "cm", res = 600, background = "white")
        print(plots[[i]])
        dev.off()
        ragg::agg_png(paste0(stem, "_preview.png"), width = width, height = height, units = "cm", res = 150, background = "white")
        print(plots[[i]])
        dev.off()
    }
    invisible(NULL)
}

plot_data <- source_read(input)

finish <- function(id) {
    file <- file.path(out, paste0(id, ".pptx"))
    layouts <- jsonlite::fromJSON(file.path(dirname(input), "layout.json"), simplifyVector = FALSE)
    a <- layouts[[id]]
    if (is.null(a) || !file.exists(file)) 
        return(invisible(file))
    tmp <- tempfile()
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE))
    unzip(file, exdir = tmp)
    ns <- c(p = "http://schemas.openxmlformats.org/presentationml/2006/main", a = "http://schemas.openxmlformats.org/drawingml/2006/main")
    path <- file.path(tmp, "ppt/presentation.xml")
    d <- xml2::read_xml(path)
    size <- xml2::xml_find_first(d, "//p:sldSz", ns)
    xml2::xml_set_attrs(size, c(cx = sprintf("%.0f", a$width), cy = sprintf("%.0f", a$height)))
    xml2::write_xml(d, path)
    for (path in list.files(file.path(tmp, 'ppt/slides'), '^slide[0-9]+\\.xml$', full.names = TRUE)) {
        d <- xml2::read_xml(path)
        symbol <- xml2::xml_find_all(d, "//*[@typeface='Symbol']")
        xml2::xml_set_attr(symbol, 'typeface', font)
        math <- xml2::xml_find_all(d, "//a:r[a:t='β' or a:t='p' or a:t='q']/a:rPr", ns)
        xml2::xml_set_attr(math, 'i', '1')
        xml2::write_xml(d, path)
    }
    if (is.character(a$transform) && length(a$transform) == 1L && nzchar(a$transform)) 
        for (path in list.files(file.path(tmp, "ppt/slides"), "^slide[0-9]+\\.xml$", full.names = TRUE)) {
            d <- xml2::read_xml(path)
            tree <- xml2::xml_find_first(d, "//p:spTree", ns)
            body <- xml2::xml_find_first(tree, "./p:grpSp[p:nvGrpSpPr/p:cNvPr[@name='figure']]", ns)
            if (inherits(body, "xml_missing")) {
                children <- as.list(xml2::xml_children(tree)[-(1:2)])
                g <- xml2::read_xml(sprintf("<p:grpSp xmlns:p=\"%s\" xmlns:a=\"%s\"><p:nvGrpSpPr><p:cNvPr id=\"1\" name=\"figure\"/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr>%s</p:grpSp>", 
                  ns[1], ns[2], a$transform))
                for (child in children) {
                  xml2::xml_add_child(g, child)
                  xml2::xml_remove(child)
                }
                xml2::xml_add_child(tree, g)
                nodes <- xml2::xml_find_all(d, "//p:cNvPr", ns)
                for (i in seq_along(nodes)) xml2::xml_set_attr(nodes[[i]], "id", as.character(i))
                xml2::write_xml(d, path)
            }
        }
    zip::zipr(file, list.files(tmp, recursive = TRUE, all.files = TRUE, no.. = TRUE), root = tmp, mode = "mirror")
    invisible(file)
}

theme_supp <- function() theme0() + theme(axis.title = element_text(size = 7), axis.text = element_text(size = 6), 
    axis.title.y = element_text(margin = margin(r = 3)), axis.title.x = element_text(margin = margin(t = 3)), strip.text = element_text(size = 7, 
        face = "plain"), strip.background = element_blank(), panel.spacing = unit(3, "mm"), legend.position = "bottom", 
    legend.key.height = unit(2.3, "mm"), legend.key.width = unit(3.5, "mm"), legend.background = element_blank(), 
    legend.key = element_blank(), plot.title = element_text(size = 7), plot.margin = margin(6, 4, 4, 6))

tag_supp <- function(plots, ncol = min(2, length(plots)), rel_heights = 1, rel_widths = 1, labels = letters[seq_along(plots)]) cowplot::plot_grid(plotlist = plots, 
    ncol = ncol, align = "hv", axis = "tblr", rel_heights = rel_heights, rel_widths = rel_widths, labels = labels, 
    label_size = 10, label_fontfamily = font, label_fontface = "bold", hjust = 0, vjust = 1, label_x = 0.008, label_y = 0.995)

route <- function(x) gsub(" to ", " → ", x, fixed = TRUE)

supp_pairs <- function() data.table(comparison = c("Low to high vs High to low", "High to high vs Low to low", 
    "High to low vs High to high", "Low to high vs Low to low"), a = c("Low to high", "High to high", "High to low", 
    "Low to high"), b = c("High to low", "Low to low", "High to high", "Low to low"))

direction_table <- function(r) {
    ref <- unique(r$results[, .(outcome, n, event)])
    ref[, `:=`(term = "urban_reference", HR = 1, lo = NA_real_, hi = NA_real_, p = NA_real_, FDR = NA_real_, sig = FALSE)]
    rbindlist(list(r$results, ref), fill = TRUE)
}

strat_spec <- function() list(tdi = list(title = "MD stratified by TDI group", var = "tdi_g", levels = tdi_lev, 
    colors = c("#79295E", "#A33D78", "#C65D94", "#DF8DB4"), drop = character()), work = list(title = "MD stratified by occupational physical demand", 
    var = "work_g", levels = c("Never/rarely", "Sometimes", "Usually", "Always"), colors = c("#176873", "#3D8F99", 
        "#72AFB7", "#A2CCD0"), drop = character()))

ppt_finish <- function(files, ...) for (f in files) finish(tools::file_path_sans_ext(basename(f)))

ed_prs <- function(id, keys) {
    r <- plot_data$prs
    names0 <- levels(r$bars$prs_label)
    labels <- ed_out(names0)
    titles <- c(md_q = "MD", dtdi_q = dtdi_lab, tdi4 = "TDI groups", urban4 = "Rural–urban routes", density_shift_q = "Density change", 
        ttwa2 = "Migration cross TTWA")
    palette <- function(key, g) {
        p <- switch(key, md_q = c("#D8F0EA", "#ACDFD2", "#39A6B1", "#015661"), dtdi_q = c("#2667AD", "#83B2D8", 
            "#E4A8AF", "#D45C69"), tdi4 = unname(pal$tdi), urban4 = c("#F0C09A", "#E99559", "#D4713B", "#AE4E25"), 
            density_shift_q = c("#DFE7ED", "#BFCFDB", "#9FB7C9", "#7F9FB8"), ttwa2 = c("#D8C5E6", "#8F65B7"))
        setNames(p, g)
    }
    plots <- list()
    for (key in keys) {
        b <- copy(r$bars[grouping == key])
        a <- copy(r$associations[grouping == key])
        stopifnot(nrow(b) > 0, nrow(a) > 0)
        lev <- unique(as.character(b$group))
        col <- palette(key, lev)
        b[, `:=`(group = factor(group, levels = lev), prs_label = factor(prs_label, levels = names0), lower = estimate - 
            1.96 * se, upper = estimate + 1.96 * se)]
        a[, `:=`(group, factor(group, levels = unique(group)))]
        lookup <- unique(b[, .(prs, prs_label)])
        a[, `:=`(prs_label, factor(as.character(lookup$prs_label[match(prs, lookup$prs)]), levels = rev(names0)))]
        p <- ggplot(b, aes(group, estimate, fill = group)) + geom_hline(yintercept = 0, colour = "#AAAAAA", linewidth = 0.14) + 
            geom_col(width = 0.72, colour = NA) + geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.14, 
            linewidth = 0.14, colour = "black") + facet_wrap(~prs_label, ncol = 3, scales = "free_y", labeller = as_labeller(setNames(labels, 
            names0))) + scale_fill_manual(values = col, labels = function(x) gsub(" to ", " → ", x, fixed = TRUE)) + 
            scale_y_continuous(breaks = scales::breaks_pretty(n = 3), expand = expansion(mult = c(0.08, 0.16))) + 
            labs(x = NULL, y = "PGS (SD)", title = titles[key], fill = NULL) + ed_theme() + theme(axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), axis.line.x = element_blank(), strip.background = element_blank(), 
            strip.text = element_text(size = 6, hjust = 0, face = "plain"), panel.spacing = grid::unit(3, "mm"), 
            plot.title = element_text(size = 7, hjust = 0.5), legend.text = element_text(size = 6), legend.key.width = grid::unit(3, 
                "mm")) + guides(fill = guide_legend(nrow = 2, byrow = TRUE))
        acol <- if (key %in% c("md_q", "dtdi_q", "density_shift_q")) 
            setNames(tail(col, 1), unique(a$group))
        else col
        pos <- position_dodge(width = 0.65)
        q <- ggplot(a, aes(beta, prs_label, colour = group, group = group)) + geom_vline(xintercept = 0, colour = "#AAAAAA", 
            linewidth = 0.14) + geom_linerange(aes(xmin = lower, xmax = upper), orientation = "y", position = pos, 
            linewidth = 0.32) + geom_point(aes(shape = q < 0.05), position = pos, size = 2.15, stroke = 0.42) + 
            scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), limits = c("FALSE", "TRUE"), labels = expression(italic(q) >= 
                0.05, italic(q) < 0.05), drop = FALSE) + scale_colour_manual(values = acol, labels = function(x) gsub(" to ", 
            " → ", x, fixed = TRUE)) + scale_y_discrete(labels = setNames(labels, names0)) + labs(x = expression(beta), 
            y = NULL, colour = NULL, shape = NULL) + ed_theme() + theme(plot.margin = margin(18, 4, 4, 2), legend.text = element_text(size = 6)) + 
            guides(colour = if (nlevels(a$group) > 1) 
                guide_legend(ncol = 1)
            else "none", shape = if (key == keys[1]) 
                guide_legend(nrow = 1, override.aes = list(colour = "black"))
            else "none")
        if (key == "tdi4") 
            p <- cowplot::plot_grid(p + theme(legend.position = "none"), ed_tdi_key(), ncol = 1, rel_heights = c(13.2, 
                0.9))
        plots[[key]] <- p
        plots[[paste0(key, "_assoc")]] <- q
    }
    order <- as.vector(rbind(keys, paste0(keys, "_assoc")))
    boxes <- data.table(x = rep(c(0.3, 19.7), 2), y = rep(c(0.3, 15), each = 2), w = rep(c(18.8, 6), 2), h = 14.1)
    ed_save(id, plots[order], boxes, 26, 29.4)
}

forest_panel <- function(d, col, show_x = TRUE) {
    lev <- names(col)
    lab <- c("All-cause death", "CVD death", "CAD", "Stroke", "HF", "AF", "T2D", "Depression", "Anxiety", "Dementia")
    d <- copy(d)
    d[, `:=`(outcome = factor(outcome, levels = forest_out), term = factor(term, levels = lev))]
    off <- setNames(seq(-0.28, 0.28, length.out = length(lev)), lev)
    d[, `:=`(xpos, as.numeric(outcome) + off[as.character(term)])]
    v <- c(d$lo, d$hi, d$HR, 1)
    v <- v[is.finite(v) & v > 0 & v < 3]
    lim <- range(v)
    pad <- max(0.02, diff(lim) * 0.08)
    lim <- c(lim[1] - pad, lim[2] + 2.7 * pad)
    ggplot(d, aes(xpos, HR, color = term)) + geom_hline(yintercept = 1, linewidth = 0.14, color = "grey45") + geom_linerange(data = d[is.finite(lo) & 
        is.finite(hi)], aes(ymin = lo, ymax = hi), linewidth = 0.24) + geom_point(aes(shape = sig), size = 1.65, 
        stroke = 0.45) + scale_color_manual(values = col, breaks = lev, drop = FALSE, name = NULL) + scale_shape_manual(values = c(`FALSE` = 1, 
        `TRUE` = 16), name = NULL, labels = c("Not significant", "FDR < 0.05")) + scale_x_continuous(NULL, breaks = seq_along(forest_out), 
        labels = lab) + scale_y_continuous("HR (95% CI)", breaks = breaks_pretty(4)) + coord_cartesian(ylim = lim) + 
        theme0() + guides(color = guide_legend(order = 1, nrow = 1, byrow = TRUE, override.aes = list(shape = 16, 
        linewidth = 0)), shape = guide_legend(order = 2, nrow = 1, override.aes = list(color = "#222222", size = 1.65))) + 
        theme(axis.text.x = if (show_x) 
            element_text(size = 6)
        else element_blank(), axis.ticks.x = if (show_x) 
            element_line(linewidth = 0.14)
        else element_blank(), legend.position = "inside", legend.position.inside = c(0.5, 0.99), legend.justification = c(0.5, 
            1), legend.box = "horizontal", legend.background = element_rect(fill = "white", color = NA), legend.key.width = unit(0.28, 
            "cm"), legend.key.height = unit(0.2, "cm"), legend.spacing.x = unit(0.1, "cm"), legend.box.spacing = unit(0.1, 
            "cm"), legend.margin = margin(0, 0, 0, 0), plot.margin = margin(1, 2, 1, 2))
}

int_curve <- function(e, t, name, outcome0, title, type, bounds = plot_data$bounds) {
    d <- copy(e[interaction == name & outcome == outcome0])
    if (type == "age") 
        d <- d[order %in% range(order)]
    lev <- switch(type, ses = c("-1 SD", "0", "+1 SD"), age = c("40", "70"), sex = c("Women", "Men"), work = strat_spec()$work$levels, 
        urban = c("Urban to Urban", "Urban to Rural", "Rural to Urban", "Rural to Rural"), density = paste0("Q", 
            1:4))
    d[, `:=`(level, factor(level, levels = lev))]
    cols <- switch(type, ses = setNames(c("#D7C4E3", "#A482BD", "#6E3F91"), lev), age = setNames(c("#EFA88A", "#B64A2F"), 
        lev), sex = setNames(c("#9FC5E8", "#315F8F"), lev), work = setNames(c("#176873", "#3D8F99", "#72AFB7", 
        "#A2CCD0"), lev), urban = setNames(c("#AD6F3B", "#D6AA78", "#E2B179", "#806044"), lev), density = setNames(c("#C3D1DE", 
        "#9BB4CA", "#6C90AE", "#3F6688"), lev))
    lg <- switch(type, ses = dtdi_lab, age = "Age", sex = "Sex", work = "Occupational\nphysical demand", urban = "Rural–urban route", 
        density = "Density change")
    ll <- if (type == "work") 
        c("Never", "Sometimes", "Usually", "Always")
    else if (type == "urban") 
        gsub(" to ", " → ", lev, fixed = TRUE)
    else lev
    v <- if (grepl("^MD", name)) 
        "md_z"
    else "dtdi_z"
    lim <- bounds[exposure == v]
    stopifnot(nrow(lim) == 1L, is.finite(lim$lower), is.finite(lim$upper), lim$lower < lim$upper)
    g <- d[, .(x = seq(lim$lower, lim$upper, length.out = 101)), by = .(level, beta, se, order)]
    linear <- unique(d$model) == "Linear"
    if (linear) 
        g[, `:=`(fit = beta * x, lo = beta * x - 1.96 * se * abs(x), hi = beta * x + 1.96 * se * abs(x))]
    else g[, `:=`(fit = exp(beta * x), lo = exp(beta * x - 1.96 * se * abs(x)), hi = exp(beta * x + 1.96 * se * 
        abs(x)))]
    span <- lim$upper - lim$lower
    right <- mean(g[x < lim$lower + 0.25 * span, fit], na.rm = TRUE) > mean(g[x > lim$upper - 0.25 * span, fit], 
        na.rm = TRUE)
    lp <- if (right) 
        c(0.05, 0.04)
    else c(0.92, 0.04)
    lj <- if (right) 
        c(0, 0)
    else c(1, 0)
    nr <- switch(type, ses = 3, work = 2, sex = 2, age = 1, urban = 4, density = 2)
    a <- t[interaction == name & outcome == outcome0][1]
    p <- ggplot(g, aes(x, fit, color = level, fill = level)) + geom_hline(yintercept = if (linear) 
        0
    else 1, linewidth = 0.12, color = "grey60") + geom_vline(xintercept = 0, linewidth = 0.12, color = "grey75") + 
        geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.08, color = NA) + geom_line(linewidth = 0.38) + scale_color_manual(values = cols, 
        labels = ll, name = lg) + scale_fill_manual(values = cols, guide = "none") + guides(color = guide_legend(nrow = nr, 
        byrow = TRUE, title.position = "top", title.hjust = 0.5)) + scale_x_continuous(breaks = scales::breaks_pretty(3), 
        expand = expansion(mult = c(0.04, 0.04))) + labs(x = paste0(if (grepl("^MD", name)) 
        "MD"
    else dtdi_lab, " (SD)"), y = if (linear) 
        "SES difference (SD)"
    else "HR") + annotate("text", x = if (right) 
        Inf
    else -Inf, y = Inf, label = title, hjust = if (right) 
        1.05
    else -0.05, vjust = 1.15, family = font, size = pt(7), color = "#222222") + annotate("text", x = if (right) 
        Inf
    else -Inf, y = Inf, label = int_p(a$p), parse = TRUE, hjust = if (right) 
        1.05
    else -0.05, vjust = 3.05, family = font, size = pt(6), color = "#222222") + scale_y_continuous(expand = expansion(mult = c(if (type %in% 
        c("work", "ses")) 0.48 else 0.32, 0.35))) + theme0(6) + theme(axis.title = element_text(size = 7), axis.text = element_text(size = 6), 
        axis.title.y = element_text(margin = margin(r = 0)), axis.text.y = element_text(margin = margin(r = 0)), 
        axis.ticks.length = unit(0.04, "cm"), legend.position = "inside", legend.position.inside = lp, legend.justification = lj, 
        legend.justification.inside = lj, legend.direction = "horizontal", legend.background = element_rect(fill = NA, 
            color = NA), legend.key = element_rect(fill = NA, color = NA), legend.box.background = element_rect(fill = NA, 
            color = NA), legend.title = element_text(size = 6), legend.text = element_text(size = 6), legend.key.width = unit(0.18, 
            "cm"), legend.key.height = unit(0.1, "cm"), legend.spacing.x = unit(0.03, "cm"), legend.margin = margin(0, 
            0, 0, 0), plot.margin = margin(1, 1, 1, 1))
    p
}

int_p <- function(x) {
    if (!is.finite(x)) 
        return("italic(p)[int] == NA")
    if (x < 0.001) 
        return("italic(p)[int] < .001")
    paste0("italic(p)[int] == ", sub("^0", "", sprintf("%.3f", x)))
}

fig2b_components <- function(cp) {
    cp <- copy(cp)
    cp[, `:=`(exposure = factor(exposure, levels = c("md", "dtdi", "ttwa"), labels = c("MD", dtdi_lab, "Migration\ncross TTWA")), 
        component = factor(component, levels = c("camsis", "income", "edu.years"), labels = c("CAMSIS", "Income", 
            "Education years")))]
    m <- max(abs(cp$estimate), na.rm = TRUE)
    ggplot(cp, aes(exposure, component, fill = estimate)) + geom_tile(color = "white", linewidth = 0.35) + geom_text(aes(label = sprintf("%.2f", 
        estimate), color = abs(estimate) > 0.62 * m), family = font, size = pt(6)) + scale_fill_gradient2(low = "#B07A4A", 
        mid = "white", high = "#5AA245", midpoint = 0, limits = c(-m, m), name = expression(italic(beta))) + scale_color_manual(values = c(`FALSE` = "#222222", 
        `TRUE` = "white"), guide = "none") + guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, 
        barwidth = unit(1.75, "cm"), barheight = unit(0.12, "cm"))) + labs(x = NULL, y = NULL) + coord_equal() + 
        theme0() + theme(axis.ticks = element_blank(), legend.position = "bottom", legend.margin = margin(0, 0, 
        0, 0), aspect.ratio = 1)
}

th <- function() theme0() + theme(axis.title = element_text(size = 7), axis.text = element_text(size = 6), axis.title.y = element_text(margin = margin(r = 3)), 
    axis.title.x = element_text(margin = margin(t = 3)), strip.text = element_text(size = 7, face = "plain"), strip.background = element_blank(), 
    panel.spacing = unit(3, "mm"), legend.position = "bottom", legend.key.height = unit(2.3, "mm"), legend.key.width = unit(3.5, 
        "mm"), legend.background = element_blank(), legend.key = element_blank(), legend.title = element_text(size = 6), 
    legend.text = element_text(size = 6), plot.title = element_text(size = 7), plot.margin = margin(6, 4, 4, 6))

cols <- function(x) scale_colour_manual(values = x, name = NULL)

fill <- function(x) scale_fill_manual(values = x, name = NULL)

tag <- function(plots, ncol = min(2, length(plots)), rel_heights = 1, rel_widths = 1, labels = letters[seq_along(plots)]) cowplot::plot_grid(plotlist = plots, 
    ncol = ncol, align = "hv", axis = "tblr", rel_heights = rel_heights, rel_widths = rel_widths, labels = labels, 
    label_size = 10, label_fontfamily = font, label_fontface = "bold", hjust = 0, vjust = 1, label_x = 0.008, label_y = 0.995)

outcome <- function(x) ed_out(x)

