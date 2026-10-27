# 1. Figure 03
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

# 1. Preserve the joint-profile classification and numerical outcomes
profile_layout <- function(cache) {
    r <- list(description = cache$description, boxes = cache$boxes,
        results = cache$outcomes, targets = c("SES", cache$targets), N = cache$N)
    r$best <- r$results[lowest == TRUE]
    routes <- c("Urban to Urban", "Urban to Rural", "Rural to Urban", "Rural to Rural")
    ord <- unlist(lapply(routes, function(v) unlist(lapply(0:1, function(k) r$description[Route == v & Cross ==
        k, profile]))))
    d <- r$description[match(ord, profile)]
    boxes <- copy(r$boxes)
    boxes[, `:=`(column, match(profile, ord))]
    W <- 23.35
    H <- 27.4
    cell <- 0.70
    xx <- 2.55 + (0:23) * 0.77 + floor((0:23)/3) * 0.08 + floor((0:23)/6) * 0.16
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
    box_panels <- list(md = list(top = 6.05, bottom = 9.15, limits = log10(c(0.5, 1000)), ticks = c(1, 10, 100,
        1000), col = teal, title = "MD (km)", sub = "log scale"), dtdi = list(top = 9.65, bottom = 12.75, limits = c(-16,
        16), ticks = c(-15, -10, -5, 0, 5, 10, 15), col = blue, title = delta, sub = ""), density_shift = list(top = 13.25,
        bottom = 16.35, limits = c(-10, 10), ticks = c(-10, -5, 0, 5, 10), col = steel, title = "Density change",
        sub = ""))
    labels <- c("SES", "All-cause\nmortality", "CVD\nmortality", "CAD", "Stroke", "HF", "AF", "T2D", "Depression",
        "Anxiety", "Dementia")
    cy <- c(17.1, 18.25 + (0:9) * 0.76)
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
    mismatch <- cache$points
    list2env(profile_layout(cache)$plot$state, envir = environment())
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
        tree <- TRUE
        yy <- function(y) if (tree) y * 1.16 else y
        txt <- function(s, x, y, size = 6.5, col = ink, just = "centre", face = "plain", maxwidth = NULL) {
            if (risk)
                face <- "plain"
            grid.text(s, x = unit(x, "cm"), y = unit(H - yy(y), "cm"), just = just, gp = gpar(fontfamily = "Arial",
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
            grid.rect(x = unit(x, "cm"), y = unit(H - yy(y), "cm"), width = unit(w, "cm"), height = unit(h, "cm"),
                gp = gpar(fill = fill, col = col, lwd = lwd))
        }
        line <- function(x, y, col = "#CBD1D6", lwd = 0.4, lty = 1) grid.lines(x = unit(x, "cm"), y = unit(H - yy(y), "cm"), gp = gpar(col = col, lwd = lwd, lty = lty))
        curve <- function(x1, y1, x2, y2, col = "#C5CDD2", lwd = 0.5) {
            grid.bezier(x = unit(c(x1, x1, x2, x2), "cm"), y = unit(H - yy(c(y1, (y1 + y2)/2, (y1 + y2)/2, y2)), "cm"),
                gp = gpar(col = col, fill = NA, lwd = lwd))
        }
        dot <- function(x, y, r = 0.035, fill = ink, col = NA, lwd = 0.3) grid.circle(x = unit(x, "cm"), y = unit(H - yy(y), "cm"), r = unit(r, "cm"), gp = gpar(fill = fill, col = col, lwd = lwd))
        diamond <- function(x, y, r = 0.065, col = ink) {
            grid.polygon(x = unit(x + c(0, r, 0, -r), "cm"), y = unit(H - y + c(r, 0, -r, 0), "cm"), gp = gpar(fill = col,
                col = "white", lwd = 0.35))
        }
        triangle <- function(x, y, fill) {
            grid.polygon(x = unit(x + c(-0.095, 0.095, -0.095), "cm"),
                y = unit(H - y - c(-0.095, -0.095, 0.095), "cm"),
                gp = gpar(fill = fill, col = "white", lwd = 0.6))
        }
        txt(paste0("All participants (N = ", fmt(r$N), ")"), mean(xx), 0.64, 7)
        txt("Rural-urban\npathway", 0.2, 1.48, 7, just = "left")
        txt("TTWA", 0.2, 2.43, 7, just = "left")
        txt(paste0("k-means (k = 3)\nMD, ", delta, ", density"), 0.2, 3.2, 6.5,
            just = "left", maxwidth = left - 0.25)
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
        txt("N", 0.2, 4.52, 7, just = "left")
        for (i in 1:24) {
            share <- d$N[i]/r$N
            cx <- xx[i]
            yp <- 4.27
            radius <- 0.22
            dot(cx, yp, radius, light, "#D9DFE3", 0.25)
            a <- seq(pi/2, pi/2 - 2 * pi * share, length.out = max(8, ceiling(360 * share)))
            grid.polygon(x = unit(c(cx, cx + radius * cos(a), cx), "cm"), y = unit(c(H - yy(yp), H - yy(yp) + radius *
                sin(a), H - yy(yp)), "cm"), gp = gpar(fill = "#647F90", col = NA))
            txt(fmt(d$N[i]), cx, 4.8, 6, maxwidth = 0.78)
        }
        tree <- FALSE
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
            txt(gsub("Density change", "Density\nchange", p$title, fixed = TRUE), 0.2,
                mean(c(p$top, p$bottom)) - 0.12, 7, p$col, just = "left", maxwidth = left - 0.45)
            if (nzchar(p$sub))
                txt(p$sub, 0.2, mean(c(p$top, p$bottom)) + 0.25, 6, just = "left")
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
        txt("Most\nfavourable\nprofile", note_x, 16.48, 7, just = "left", maxwidth = W - note_x - 0.3)
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
                if (j > 1) {
                  mark <- mismatch[profile == ord[i] & outcome == r$targets[j]]
                  if (nrow(mark) && isTRUE(mark$discordant))
                    triangle(xx[i] - cell/2 + 0.115, cy[j] - cell/2 + 0.115,
                      if (mark$pattern == "Higher SES / higher risk") "#D45C69" else "#8F65B7")
                }
                if (record)
                  geometry$squares[[length(geometry$squares) + 1L]] <<- list(profile = ord[i], outcome = r$targets[j],
                    x = xx[i], y = cy[j], w = cell, h = cell, estimate = value, optimum = isTRUE(a$lowest))
            }
            a <- r$best[outcome == r$targets[j]]
            txt(a$profile, note_x, cy[j], 7, face = "bold", just = "left")
        }
        ly <- 26.1
        bar <- function(pal, x, y, width, ticks, title) {
            txt(title, x, y - 0.35, 7, just = "left")
            grid.rect(x = unit(x + (seq_along(pal) - 0.5) * width/length(pal), "cm"),
                y = unit(H - yy(y), "cm"), width = unit(width/length(pal) + 0.0005, "cm"),
                height = unit(0.23, "cm"), gp = gpar(fill = pal, col = NA))
            for (i in seq_along(ticks)) txt(ticks[i], x + (i - 1) * width/(length(ticks) - 1), y + 0.31, 6.2)
        }
        bar(ses_pal, 2.4, ly, 3.1, sgn(c(ses_range[1], mean(ses_range), ses_range[2])), if (risk)
            "Adjusted SES (SD)"
        else expression(italic(beta)[SES]))
        bar(hr_pal, 6.9, ly, 3.1, if (risk)
            c("Lower", "", "Higher")
        else sprintf("%.2f", c(0.5, 1, 2)), if (risk)
            "10-year risk (%)"
        else paste0("HR (ref. ", r$reference, ")"))
        rect(11.3, ly, 0.33, 0.33, "white", "#17232C", 1.25)
        diamond(11.405, ly - 0.105, 0.055)
        txt("Most favourable\npoint estimate", 11.68, ly, 6.2, just = "left")
        rect(15.15, ly, 0.33, 0.33, "#E8EBEE", "#D5DBDF", 0.24)
        txt("Insufficient\ndata", 15.53, ly, 6.2, just = "left")
        keyx <- 18.65
        keyy <- ly
        line(c(keyx, keyx), c(keyy - 0.38, keyy + 0.38), "#64727D", 0.48)
        line(keyx + c(-0.1, 0.1), rep(keyy - 0.38, 2), "#64727D", 0.48)
        line(keyx + c(-0.1, 0.1), rep(keyy + 0.38, 2), "#64727D", 0.48)
        rect(keyx, keyy, 0.36, 0.4, "#BECBD3", NA)
        line(keyx + c(-0.18, 0.18), rep(keyy, 2), ink, 0.72)
        txt(paste0("Median; box: Q1-Q3\nWhiskers: 1.5 ", times, " IQR"), 19.12, keyy, 6.2, just = "left", maxwidth = 5.1)
        triangle(2.8, 26.85, "#D45C69")
        txt("Higher SES / higher risk", 3.04, 26.85, 6.2, just = "left")
        triangle(8.9, 26.85, "#8F65B7")
        txt("Lower SES / lower risk", 9.14, 26.85, 6.2, just = "left")
    }
}
mismatch_draw <- function(cache) {
  r <- cache
  H <- 28.0
  grid::grid.newpage()
  original <- profile_draw
  colours <- c('Higher SES / lower risk'='#39A6B1',
               'Higher SES / higher risk'='#D45C69',
               'Lower SES / lower risk'='#8F65B7',
               'Lower SES / higher risk'='#AAB3BB')
  labels <- c('All-cause mortality','CVD mortality','CAD','Stroke','HF',
              'AF','T2D','Depression','Anxiety','Dementia')
  ink <- '#28343D'
  maximum <- max(r$description$N)
  minus <- intToUtf8(8722)
  fmt_x <- function(x) gsub('-',minus,format(x,trim=TRUE,scientific=FALSE),fixed=TRUE)
  
  bubble <- function(target,title) {
    d <- copy(r$points[outcome==target & available==TRUE])
    stopifnot(nrow(d)==21L,uniqueN(d$profile)==21L)
    sx <- d$ses_reference[1]; sy <- d$risk_reference[1]
    xr <- c(-.73,.66)
    yr <- range(d$estimate)
    span <- diff(yr)
    yr <- yr + c(-.16,.20)*span
    yr[1] <- max(0,yr[1])
    mark <- unique(c(d[discordant==TRUE,profile],d[which.max(ses),profile],d[which.min(estimate),profile]))
    d[,label:=fifelse(profile%in%mark,profile,'')]
    p <- ggplot(d,aes(ses,estimate)) +
      annotate('rect',xmin=sx,xmax=xr[2],ymin=sy,ymax=yr[2],
               fill=colours[['Higher SES / higher risk']],alpha=.075) +
      annotate('rect',xmin=xr[1],xmax=sx,ymin=yr[1],ymax=sy,
               fill=colours[['Lower SES / lower risk']],alpha=.075) +
      geom_vline(xintercept=sx,colour='#85939D',linewidth=.17,linetype=2) +
      geom_hline(yintercept=sy,colour='#85939D',linewidth=.17,linetype=2) +
      geom_point(aes(size=N,fill=pattern),shape=21,colour='white',stroke=.15,alpha=.94) +
      geom_text_repel(aes(label=label),family=font,size=6/ggplot2::.pt,
        colour=ink,box.padding=.12,point.padding=.12,min.segment.length=0,
        segment.colour='#A0ACB5',segment.size=.15,max.overlaps=Inf,
        force=1.5,force_pull=.5,max.iter=20000,max.time=Inf,
        seed=20261027,show.legend=FALSE) +
      scale_fill_manual(values=colours,drop=FALSE) +
      scale_size_area(max_size=3.6,limits=c(0,maximum)) +
      scale_x_continuous(breaks=c(-.5,0,.5),labels=fmt_x,expand=expansion(mult=0)) +
      scale_y_continuous(breaks=function(x)pretty(x,n=3),
        labels=function(x)formatC(x,digits=1,format='f'),expand=expansion(mult=0)) +
      coord_cartesian(xlim=xr,ylim=yr,clip='on') +
      labs(title=title,x=NULL,y=NULL) +
      theme_classic(base_size=6,base_family=font) +
      theme(text=element_text(colour=ink),
        plot.title=element_text(size=7,face='plain',hjust=.5,margin=margin(b=3)),
        axis.text=element_text(size=6,colour=ink),
        axis.line=element_line(linewidth=.18,colour='#657580'),
        axis.ticks=element_line(linewidth=.16,colour='#657580'),
        axis.ticks.length=grid::unit(.065,'cm'),
        axis.text.x=element_text(margin=margin(t=2)),
        axis.text.y=element_text(margin=margin(r=2)),
        legend.position='none',aspect.ratio=1,
        plot.margin=margin(2,2,2,2),plot.background=element_blank())
    ggplotGrob(p)
  }
  plots <- lapply(seq_along(r$targets),function(i)bubble(r$targets[i],labels[i]))
  widths <- do.call(grid::unit.pmax,lapply(plots,`[[`,'widths'))
  heights <- do.call(grid::unit.pmax,lapply(plots,`[[`,'heights'))
  for(i in seq_along(plots)) {
    plots[[i]]$widths <- widths
    plots[[i]]$heights <- heights
  }
  at <- function(x,y,w,h,hh) grid::viewport(x=grid::unit(x,'cm'),y=grid::unit(hh-y,'cm'),
    width=grid::unit(w,'cm'),height=grid::unit(h,'cm'),just=c('left','top'),clip='off')
  text <- function(value,x,y,hh,size=6.5,just='centre',face='plain',rot=0)
    grid::grid.text(value,x=grid::unit(x,'cm'),y=grid::unit(hh-y,'cm'),just=just,rot=rot,
      gp=grid::gpar(fontfamily='Arial',fontsize=size,col=ink,fontface=face,lineheight=1.05))
  grid::pushViewport(at(.3,.3,23.35,27.4,H))
  original(r)
  grid::popViewport()
  grid::pushViewport(at(23.73,.3,10,27.4,H))
  bh <- 27.4
  cw <- (10-.15-.12)/2
  ch <- 5.05
  gap <- .08
  for(i in seq_along(plots)) {
    x <- .15+((i-1)%%2)*(cw+.12)
    y <- floor((i-1)/2)*(ch+gap)
    grid::pushViewport(at(x,y,cw,ch,bh))
    grid::grid.draw(plots[[i]])
    grid::popViewport()
  }
  bottom <- 5*ch+4*gap
  text('SES (SD)',5,bottom+.20,bh,7)
  text('10-year risk (%)',.015,bottom/2,bh,7,rot=90)
  ly <- bottom+.73
  for(i in 1:4) {
    x <- .35+(i-1)*2.4
    grid::grid.circle(grid::unit(x,'cm'),grid::unit(bh-ly,'cm'),r=grid::unit(.065,'cm'),
      gp=grid::gpar(fill=colours[i],col=NA))
    text(gsub(' / ','\n',names(colours)[i],fixed=TRUE),x+.15,ly,bh,6,just='left')
  }
  ly <- bottom+1.31
  text('N',.8,ly,bh,7)
  for(i in 1:3) {
    x <- 1.8+(i-1)*2.65
    value <- c(5000,25000,75000)[i]
    grid::grid.circle(grid::unit(x,'cm'),grid::unit(bh-ly,'cm'),r=grid::unit(.18*sqrt(value/maximum),'cm'),
      gp=grid::gpar(fill='#8C9CA8',col='white',lwd=.4))
    text(c('5,000','25,000','75,000')[i],x+.25,ly,bh,6,just='left')
  }
  grid::popViewport()
  text('a',.16,.18,H,10,just='left',face='bold')
  text('b',23.75,.18,H,10,just='left',face='bold')
}

r <- plot_data$result
draw_ppt(file.path(out,'Fig3.pptx'),function()mismatch_draw(r),34.03*cm,28.0*cm)
ragg::agg_png(file.path(out,'Fig3.png'),width=34.03,height=28.0,units='cm',res=600)
mismatch_draw(r);dev.off()
grDevices::cairo_pdf(file.path(out,'Fig3.pdf'),width=34.03*cm,height=28.0*cm,family=font)
mismatch_draw(r);dev.off()
finish('Fig3')
