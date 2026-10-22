# 1. Extended Data Figure 8
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

profile_layout <- function(r) {
    routes <- c("Urban to Urban", "Urban to Rural", "Rural to Urban", "Rural to Rural")
    ord <- unlist(lapply(routes, function(v) unlist(lapply(0:1, function(k) r$description[Route == v & Cross == 
        k, profile]))))
    d <- r$description[match(ord, profile)]
    boxes <- copy(r$boxes)
    boxes[, `:=`(column, match(profile, ord))]
    W <- 26.7
    H <- 25.15
    cell <- 0.77
    xx <- 2.85 + (0:23) * 0.84 + floor((0:23)/3) * 0.1 + floor((0:23)/6) * 0.18
    left <- xx[1] - cell/2
    right <- tail(xx, 1) + cell/2
    note_x <- right + 0.18
    delta <- dtdi_lab
    arrow <- intToUtf8(8594)
    minus <- intToUtf8(8722)
    times <- intToUtf8(215)
    ink <- "#25323A"
    teal <- "#168D83"
    blue <- "#2667AD"
    red <- "#D45C69"
    purple <- "#8F65B7"
    steel <- "#325375"
    light <- "#EFF1F3"
    route_fill <- route_stroke <- c("#F0C09A", "#E99559", "#D4713B", "#AE4E25")
    hr_pal <- colorRampPalette(c(blue, "#F8FAFC", red))(401)
    ses_pal <- colorRampPalette(c("#F3F8F5", "#A7D3C4", teal))(401)
    hr_range <- c(0.5, 2)
    ses_range <- range(pretty(r$results[outcome == "SES" & supported == TRUE, estimate], n = 4))
    vars <- c("md", "dtdi", "density_shift")
    box_panels <- list(md = list(top = 5.35, bottom = 7.4, limits = log10(c(0.5, 1000)), ticks = c(1, 10, 100, 
        1000), col = teal, title = "MD (km)", sub = "log scale"), dtdi = list(top = 7.95, bottom = 10, limits = c(-16, 
        16), ticks = c(-15, -10, -5, 0, 5, 10, 15), col = blue, title = delta, sub = ""), density_shift = list(top = 10.55, 
        bottom = 12.6, limits = c(-10, 10), ticks = c(-10, -5, 0, 5, 10), col = steel, title = "Density change", 
        sub = ""))
    labels <- c("SES", "All-cause death", "CVD death", "CAD", "Stroke", "Heart failure", "AF", "T2D", "Depression", 
        "Anxiety", "Dementia")
    cy <- c(13.65, 14.8 + (0:9) * 0.84)
    geometry <- list(text = list(), squares = list(), regional = list())
    risk <- TRUE
    regional <- FALSE
    joint <- NULL
    names <- c("routes", "ord", "d", "vars", "boxes", "W", "H", "cell", "xx", "left", "right", "note_x", "delta", 
        "arrow", "minus", "times", "ink", "teal", "blue", "red", "purple", "steel", "light", "route_fill", "route_stroke", 
        "hr_pal", "ses_pal", "hr_range", "ses_range", "box_panels", "labels", "cy", "geometry", "risk", "regional", 
        "r", "joint")
    r$plot <- list(state = mget(names, envir = environment()))
    r
}

profile_draw <- function(cache) {
    list2env(cache$plot$state, envir = environment())
    for (v in c("convertWidth", "gpar", "grid.bezier", "grid.circle", "grid.lines", "grid.newpage", "grid.polygon", 
        "grid.raster", "grid.rect", "grid.text", "popViewport", "pushViewport", "stringWidth", "unit", "viewport")) assign(v, 
        getExportedValue("grid", v), envir = environment())
    record <- FALSE
    fmt <- function(v) format(v, big.mark = ",", trim = TRUE, scientific = FALSE)
    sgn <- function(v, n = 2) {
        v[is.finite(v) & abs(v) < 0.5 * 10^-n] <- 0
        gsub("-", minus, sprintf(paste0("%.", n, "f"), v), fixed = TRUE)
    }
    colour <- function(value, type, limits = NULL) {
        if (!is.finite(value)) 
            return("#E8EBEE")
        t <- if (type == "Risk") 
            (value - limits[1])/diff(limits)
        else if (type == "HR") 
            (log(value) - log(hr_range[1]))/diff(log(hr_range))
        else (value - ses_range[1])/diff(ses_range)
        pal <- if (type %in% c("HR", "Risk")) 
            hr_pal
        else ses_pal
        pal[1 + round(max(0, min(1, t)) * (length(pal) - 1))]
    }
    {
        grid.newpage()
        txt <- function(s, x, y, size = 6.5, col = ink, just = "centre", face = "plain", maxwidth = NULL) {
            if (risk) 
                face <- "plain"
            grid.text(s, x = unit(x, "cm"), y = unit(H - y, "cm"), just = just, gp = gpar(fontfamily = "Arial", 
                fontsize = size, col = col, fontface = face, lineheight = 1.03))
            if (!is.null(maxwidth) && is.character(s)) {
                pushViewport(viewport(gp = gpar(fontfamily = "Arial", fontsize = size, fontface = face)))
                widths <- vapply(strsplit(s, "\n", fixed = TRUE)[[1]], function(a) convertWidth(stringWidth(a), 
                  "cm", valueOnly = TRUE), 0)
                popViewport()
                stopifnot(max(widths) <= maxwidth)
            }
            if (record) 
                geometry$text[[length(geometry$text) + 1L]] <<- list(label = as.character(s), x = x, y = y, size = size, 
                  just = just, face = face)
        }
        rect <- function(x, y, w, h, fill = "white", col = NA, lwd = 0.35) {
            grid.rect(x = unit(x, "cm"), y = unit(H - y, "cm"), width = unit(w, "cm"), height = unit(h, "cm"), 
                gp = gpar(fill = fill, col = col, lwd = lwd))
        }
        line <- function(x, y, col = "#CBD1D6", lwd = 0.4, lty = 1) grid.lines(x = unit(x, "cm"), y = unit(H - 
            y, "cm"), gp = gpar(col = col, lwd = lwd, lty = lty))
        curve <- function(x1, y1, x2, y2, col = "#C5CDD2", lwd = 0.5) {
            grid.bezier(x = unit(c(x1, x1, x2, x2), "cm"), y = unit(H - c(y1, (y1 + y2)/2, (y1 + y2)/2, y2), "cm"), 
                gp = gpar(col = col, fill = NA, lwd = lwd))
        }
        dot <- function(x, y, r = 0.035, fill = ink, col = NA, lwd = 0.3) grid.circle(x = unit(x, "cm"), y = unit(H - 
            y, "cm"), r = unit(r, "cm"), gp = gpar(fill = fill, col = col, lwd = lwd))
        diamond <- function(x, y, r = 0.065, col = ink) {
            grid.polygon(x = unit(x + c(0, r, 0, -r), "cm"), y = unit(H - y + c(r, 0, -r, 0), "cm"), gp = gpar(fill = col, 
                col = "white", lwd = 0.35))
        }
        txt(paste0("All participants (N = ", fmt(r$N), ")"), mean(xx), 0.64, 7)
        txt("Rural-urban\npathway", note_x, 1.48, 7, just = "left")
        txt("TTWA", note_x, 2.43, 7, just = "left")
        txt(paste0("k-means (k = 3)\nMD, ", delta, ", density"), note_x, 3.2, 6.5, just = "left", maxwidth = W - 
            note_x - 0.3)
        for (j in 1:4) {
            ii <- (j - 1) * 6 + 1:6
            mids <- mean(range(xx[ii]))
            lo <- min(xx[ii]) - cell/2
            hi <- max(xx[ii]) + cell/2
            curve(mean(xx), 0.8, mids, 1.19)
            rect(mids, 1.48, hi - lo, 0.52, route_fill[j], route_stroke[j], 0.45)
            txt(gsub(" to ", paste0(" ", arrow, " "), routes[j], fixed = TRUE), mids, 1.48, 7, if (j >= 3) 
                "white"
            else ink)
            for (k in 1:2) {
                ix <- ii[(k - 1) * 3 + 1:3]
                mid <- mean(range(xx[ix]))
                curve(mids, 1.74, mid, 2.22, route_stroke[j], 0.5)
                rect(mid, 2.43, diff(range(xx[ix])) + cell, 0.38, if (k == 1) 
                  "#F8F5FC"
                else "#F0EAF8", NA)
                txt(if (k == 1) 
                  "Same TTWA"
                else "Cross TTWA", mid, 2.43, 6.5, purple)
                for (i in ix) {
                  curve(mid, 2.64, xx[i], 3.3, "#AAB6C0", 0.42)
                  txt(ord[i], xx[i], 3.62, 7, maxwidth = 0.76)
                }
            }
        }
        txt("N", note_x, 4.52, 7, just = "left")
        for (i in 1:24) {
            share <- d$N[i]/r$N
            cx <- xx[i]
            yp <- 4.27
            radius <- 0.22
            dot(cx, yp, radius, light, "#D9DFE3", 0.25)
            a <- seq(pi/2, pi/2 - 2 * pi * share, length.out = max(8, ceiling(360 * share)))
            grid.polygon(x = unit(c(cx, cx + radius * cos(a), cx), "cm"), y = unit(c(H - yp, H - yp + radius * 
                sin(a), H - yp), "cm"), gp = gpar(fill = "#647F90", col = NA))
            txt(fmt(d$N[i]), cx, 4.8, 6, maxwidth = 0.78)
        }
        for (v in vars) {
            p <- box_panels[[v]]
            trans <- function(value) if (v == "md") 
                log10(value)
            else value
            ypos <- function(value) p$bottom - (trans(value) - p$limits[1])/diff(p$limits) * (p$bottom - p$top)
            for (cx in xx) line(rep(cx, 2), c(p$top, p$bottom), "#D4DADF", 0.36, 3)
            for (t in p$ticks) {
                y <- ypos(t)
                if (t == 0 && v != "md") 
                  line(c(left, right), c(y, y), "#BFC9D0", 0.4, 2)
                axis_x <- if (risk) 
                  left - 0.05
                else right
                line(axis_x + if (risk) 
                  c(-0.06, 0)
                else c(0, 0.06), c(y, y), "#6D7983", 0.4)
                txt(gsub("-", minus, as.character(t), fixed = TRUE), if (risk) 
                  axis_x - 0.14
                else right + 0.12, y, 6, just = if (risk) 
                  "right"
                else "left")
            }
            axis_x <- if (risk) 
                left - 0.05
            else right
            line(c(axis_x, axis_x), c(p$top, p$bottom), "#6D7983", 0.4)
            txt(p$title, note_x, mean(c(p$top, p$bottom)) - 0.12, 7, p$col, just = "left", maxwidth = W - note_x - 
                0.3)
            if (nzchar(p$sub)) 
                txt(p$sub, note_x, mean(c(p$top, p$bottom)) + 0.25, 6, just = "left")
            b <- boxes[variable == v][order(column)]
            for (i in 1:24) {
                z <- b[i]
                cx <- xx[i]
                col <- if (v == "dtdi") 
                  if (z$median < 0) 
                    blue
                  else red
                else p$col
                rgb0 <- as.numeric(col2rgb(col))/255
                fill <- if (v == "density_shift") 
                  col
                else rgb(t(0.55 * rgb0 + 0.45))
                whisker <- adjustcolor(col, alpha.f = 0.72)
                line(c(cx, cx), c(ypos(z$lower), ypos(z$upper)), whisker, 0.48)
                for (k in c(z$lower, z$upper)) line(cx + c(-0.11, 0.11), rep(ypos(k), 2), whisker, 0.48)
                hh <- abs(ypos(z$q3) - ypos(z$q1))
                rect(cx, mean(c(ypos(z$q1), ypos(z$q3))), 0.43, max(0.012, hh), fill, NA)
                lum <- sum(as.numeric(col2rgb(fill))/255 * c(0.2126, 0.7152, 0.0722))
                line(cx + c(-0.215, 0.215), rep(ypos(z$median), 2), if (lum < 0.48) 
                  "white"
                else ink, 0.72)
            }
        }
        txt("Most favourable\nprofile", note_x, 13.03, 7, just = "left", maxwidth = W - note_x - 0.3)
        for (j in seq_along(r$targets)) {
            txt(labels[j], left - 0.17, cy[j], 6.5, just = "right")
            z <- r$results[outcome == r$targets[j]]
            limits <- if (risk && j > 1) 
                range(z[supported == TRUE, estimate])
            else NULL
            for (i in 1:24) {
                a <- z[profile == ord[i]]
                value <- if (nrow(a) && isTRUE(a$supported)) 
                  a$estimate
                else NA_real_
                fill <- colour(value, if (j == 1) 
                  "SES"
                else if (risk) 
                  "Risk"
                else "HR", limits)
                rect(xx[i], cy[j], cell, cell, fill, "#D5DBDF", 0.24)
                optimum <- nrow(a) && isTRUE(a$lowest)
                if (optimum) {
                  rect(xx[i], cy[j], cell - 0.018, cell - 0.018, NA, "#17232C", 1.55)
                  rect(xx[i], cy[j], cell - 0.085, cell - 0.085, NA, "white", 0.7)
                  diamond(xx[i] + cell/2 - 0.075, cy[j] - cell/2 + 0.075, 0.065)
                }
                if (is.finite(value)) {
                  lum <- sum(as.numeric(col2rgb(fill))/255 * c(0.2126, 0.7152, 0.0722))
                  txt(sgn(value), xx[i], cy[j], if (optimum) 
                    7
                  else 6, if (lum < 0.5) 
                    "white"
                  else ink, face = if (optimum) 
                    "bold"
                  else "plain", maxwidth = 0.62)
                }
                if (record) 
                  geometry$squares[[length(geometry$squares) + 1L]] <<- list(profile = ord[i], outcome = r$targets[j], 
                    x = xx[i], y = cy[j], w = cell, h = cell, estimate = value, optimum = isTRUE(a$lowest))
            }
            a <- r$best[outcome == r$targets[j]]
            txt(a$profile, note_x, cy[j], 7, face = "bold", just = "left")
        }
        if (regional) {
            metrics <- joint$metrics
            rx <- c(left + 0.9, seq(left + 3.15, right - 0.65, length.out = 10))
            sd_ticks <- pretty(c(0, max(metrics[outcome == "SES", estimate]) * 1.15), n = 3)
            sd_ticks <- sd_ticks[sd_ticks >= 0]
            hr_ticks <- pretty(c(1, 1 + (max(metrics[outcome != "SES", estimate]) - 1) * 1.25), n = 3)
            hr_ticks <- hr_ticks[hr_ticks >= 1]
            sd_lim <- range(sd_ticks)
            hr_lim <- range(hr_ticks)
            region_col <- c(origin = "#018571", destination = "#9D3978")
            txt("Regional SD", rx[1], 23.48, 6.5)
            txt("Regional MHR", mean(range(rx[-1])), 23.48, 6.5)
            for (k in 1:2) {
                ctx <- c("origin", "destination")[k]
                top <- c(24.08, 27)[k]
                bottom <- top + 1.85
                yval <- function(value, ses = FALSE) {
                  lim <- if (ses) 
                    sd_lim
                  else hr_lim
                  bottom - (value - lim[1])/diff(lim) * (bottom - top)
                }
                for (ses in c(TRUE, FALSE)) {
                  ticks <- if (ses) 
                    sd_ticks
                  else hr_ticks
                  ax <- if (ses) 
                    left + 0.02
                  else left + 2.43
                  line(rep(ax, 2), c(top, bottom), "#6D7983", 0.4)
                  for (t in ticks) {
                    yy <- yval(t, ses)
                    line(ax + c(-0.05, 0), rep(yy, 2), "#6D7983", 0.4)
                    txt(sprintf(if (ses) 
                      "%.2f"
                    else "%.2f", t), ax - 0.1, yy, 6, just = "right")
                    if (!ses) 
                      line(c(ax, right), rep(yy, 2), "#E7EBEE", 0.3)
                  }
                }
                for (j in seq_along(r$targets)) {
                  z <- metrics[outcome == r$targets[j] & context == ctx]
                  stopifnot(nrow(z) == 2L, setequal(z$stage, c("without", "with")))
                  before <- z[stage == "without", estimate]
                  after <- z[stage == "with", estimate]
                  yy <- yval(c(before, after), j == 1)
                  px <- rx[j] + c(-0.16, 0.16)
                  line(px, yy, "#9AA8B2", 0.65)
                  dot(px[1], yy[1], 0.055, "#B6BEC5")
                  dot(px[2], yy[2], 0.065, region_col[[ctx]])
                  txt(sprintf("%.2f", after), px[2] + 0.02, yy[2] - 0.2, 6, region_col[[ctx]])
                  if (record) 
                    geometry$regional[[length(geometry$regional) + 1L]] <<- list(outcome = r$targets[j], context = ctx, 
                      x = rx[j], top = top, bottom = bottom, metric = if (j == 1) "Regional SD" else "MHR", before = before, 
                      after = after)
                }
                txt(if (k == 1) 
                  "Birthplace"
                else "Recruitment", note_x, mean(c(top, bottom)), 7, region_col[[ctx]], just = "left")
                line(rep(left + 1.95, 2), c(top - 0.05, bottom + 0.05), "#DBE1E5", 0.35)
            }
            short <- c("SES", "All-cause\ndeath", "CVD\ndeath", "CAD", "Stroke", "Heart\nfailure", "AF", "T2D", 
                "Depression", "Anxiety", "Dementia")
            for (j in seq_along(short)) txt(short[j], rx[j], 29.35, 6, maxwidth = 1.65)
            dot(11.7, 31.62, 0.055, "#B6BEC5")
            txt("Without profiles", 11.9, 31.62, 6.2, just = "left")
            dot(13.9, 31.62, 0.065, region_col[["origin"]])
            dot(14.09, 31.62, 0.065, region_col[["destination"]])
            txt("With 24 profiles", 14.3, 31.62, 6.2, just = "left")
        }
        ly <- if (regional) 
            30.55
        else 24
        bar <- function(pal, x, y, width, ticks, title) {
            txt(title, x, y - 0.35, 7, just = "left")
            grid.raster(as.raster(matrix(pal, nrow = 1)), x = unit(x + width/2, "cm"), y = unit(H - y, "cm"), width = unit(width, 
                "cm"), height = unit(0.23, "cm"), interpolate = TRUE)
            for (i in seq_along(ticks)) txt(ticks[i], x + (i - 1) * width/(length(ticks) - 1), y + 0.31, 6.2)
        }
        bar(ses_pal, 2.6, ly, 3.8, sgn(c(ses_range[1], mean(ses_range), ses_range[2])), if (risk) 
            "Adjusted SES (SD)"
        else expression(italic(beta)[SES]))
        bar(hr_pal, 7.7, ly, 3.8, if (risk) 
            c("Lower", "", "Higher")
        else sprintf("%.2f", c(0.5, 1, 2)), if (risk) 
            "10-year risk (%)"
        else paste0("HR (ref. ", r$reference, ")"))
        rect(13, ly, 0.33, 0.33, "white", "#17232C", 1.25)
        diamond(13.105, ly - 0.105, 0.055)
        txt("Most favourable\npoint estimate", 13.38, ly, 6.2, just = "left")
        rect(17, ly, 0.33, 0.33, "#E8EBEE", "#D5DBDF", 0.24)
        txt("Insufficient\ndata", 17.38, ly, 6.2, just = "left")
        keyx <- 20.95
        keyy <- ly
        line(c(keyx, keyx), c(keyy - 0.38, keyy + 0.38), "#64727D", 0.48)
        line(keyx + c(-0.1, 0.1), rep(keyy - 0.38, 2), "#64727D", 0.48)
        line(keyx + c(-0.1, 0.1), rep(keyy + 0.38, 2), "#64727D", 0.48)
        rect(keyx, keyy, 0.36, 0.4, "#BECBD3", NA)
        line(keyx + c(-0.18, 0.18), rep(keyy, 2), ink, 0.72)
        txt(paste0("Median; box: Q1-Q3\nWhiskers: 1.5 ", times, " IQR"), 21.42, keyy, 6.2, just = "left", maxwidth = 5.1)
    }
}

r <- plot_data$result

if (is.null(r$plot)) r <- profile_layout(r)

draw_ppt(file.path(out, "ED8.pptx"), function() profile_draw(r), r$plot$state$W * cm, r$plot$state$H * cm)

ragg::agg_png(file.path(out, "ED8.png"), width = r$plot$state$W, height = r$plot$state$H, units = "cm", res = 300)

profile_draw(r)

dev.off()

finish("ED8")

