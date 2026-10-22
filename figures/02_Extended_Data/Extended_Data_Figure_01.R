# 1. Extended Data Figure 1
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1])))), "style.R"))

ed1_plot <- function(r) {
    cols <- c(UKB = "#5F99D2", UKHLS = "#B76585", HRS = "#39A6B1")
    W <- 24
    H <- 17.8
    txt <- function(s, x, y, size = 6, col = "#222222", face = "plain", just = "centre") grid::grid.text(s, grid::unit(x, 
        "cm"), grid::unit(H - y, "cm"), just = just, gp = grid::gpar(fontfamily = font, fontsize = size, col = col, 
        fontface = face))
    line <- function(x, y, arrow = FALSE) grid::grid.lines(grid::unit(x, "cm"), grid::unit(H - y, "cm"), gp = grid::gpar(col = "#7A7A7A", 
        lwd = 0.4, linejoin = "round"), arrow = if (arrow) 
        grid::arrow(length = grid::unit(0.075, "cm"), type = "closed")
    else NULL)
    rect <- function(x, y, w, h, fill, border, lwd = 0.4) {
        grid::grid.rect(grid::unit(x, "cm"), grid::unit(H - y - h/2, "cm"), grid::unit(w, "cm"), grid::unit(h, 
            "cm"), gp = grid::gpar(fill = fill, col = border, lwd = lwd))
    }
    note <- function(s, x, y, size = 6, col = "#666666") {
        g <- grid::textGrob(s, gp = grid::gpar(fontfamily = font, fontsize = size))
        w <- grid::convertWidth(grid::grobWidth(g), "cm", valueOnly = TRUE) + 0.12
        rect(x, y - 0.13, w, 0.26, "white", NA)
        txt(s, x, y, size, col)
    }
    box <- function(x, y, w, h, col, lines, sizes = NULL, final = FALSE) {
        rect(x, y, w, h, if (final) 
            scales::alpha(col, 0.09)
        else "white", scales::alpha(col, 0.65))
        if (final) 
            rect(x, y, w, 0.065, col, NA)
        gap <- min(0.31, (h - 0.25)/max(length(lines) - 1, 1))
        yy <- y + h/2 + (seq_along(lines) - (length(lines) + 1)/2) * gap
        if (is.null(sizes)) 
            sizes <- rep(6, length(lines))
        for (k in seq_along(lines)) txt(lines[k], x, yy[k], sizes[k], if (grepl("^N =", lines[k])) 
            col
        else "#222222")
    }
    num <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)
    N <- function(n) paste0("N = ", num(n))
    exclude <- function(x, y, w, s) {
        rect(x, y, w, 0.76, "#F7F7F7", "#D7D7D7")
        for (i in seq_along(s)) txt(s[i], x, y + 0.38 + (i - (length(s) + 1)/2) * 0.27, 6, col = "#555555")
    }
    result <- function(x, y, w, h, col, title, n, migrant, non, md = NULL) {
        rect(x, y, w, h, scales::alpha(col, 0.09), scales::alpha(col, 0.65))
        rect(x, y, w, 0.065, col, NA)
        txt(title, x, y + 0.27, 7)
        txt(N(n), x, y + 0.62, 7, col)
        line(c(x - w/2 + 0.18, x + w/2 - 0.18), rep(y + 0.87, 2))
        txt("Migrants", x - w/2 + 0.2, y + 1.1, 6, just = "left")
        txt(num(migrant), x + w/2 - 0.2, y + 1.1, 6, just = "right")
        txt("Non-migrants", x - w/2 + 0.2, y + 1.4, 6, just = "left")
        txt(num(non), x + w/2 - 0.2, y + 1.4, 6, just = "right")
        if (!is.null(md)) 
            txt(paste0("Migrants with MD > 0: ", num(md)), x, y + 1.73, 6)
    }
    function() {
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(width = grid::unit(W, "cm"), height = grid::unit(H, "cm")))
        grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
        line(c(13.65, 13.65), c(0.35, 16.8))
        d <- r$display$UKB
        col <- cols["UKB"]
        txt("a", 0.5, 0.4, 10, face = "bold", just = "left")
        txt("UK Biobank", 1.1, 0.42, 7, col, face = "bold", just = "left")
        box(4, 1, 6.7, 1.2, col, c("Available baseline participants", N(d$n[1])), c(6, 7))
        box(4, 2.85, 6.7, 1.2, col, c("White participants", N(d$n[2])), c(6, 7))
        box(4, 4.7, 6.7, 1.3, col, c("Valid Birthplace and Recruitment", "residence coordinates", N(d$n[3])), c(6, 
            6, 7))
        line(c(4, 4), c(2.2, 2.77), TRUE)
        line(c(4, 8), rep(2.53, 2), TRUE)
        exclude(10.45, 2.15, 4.8, c("Non-White or ethnicity missing", paste0("Excluded: ", num(d$n[1] - d$n[2]))))
        line(c(4, 4), c(4.05, 4.62), TRUE)
        line(c(4, 8), rep(4.38, 2), TRUE)
        exclude(10.45, 4, 4.8, c("Geographic coordinates missing", paste0("Excluded: ", num(d$n[2] - d$n[3]))))
        line(c(4, 4, 3.6, 3.6), c(6, 6.3, 6.3, 7.3), TRUE)
        line(c(4, 10.2, 10.2), c(6.3, 6.3, 7.3), TRUE)
        note("Genetic analyses", 3.6, 6.75, 7, col)
        note("Phenotypic analyses", 10.2, 6.75, 7, col)
        box(3.6, 7.4, 5.6, 1.25, col, c("Genotypes, exposure and genetic-model data", "Age, age squared, sex, centre, batch", 
            "and PC1-PC40 complete"), rep(6, 3))
        line(c(3.6, 3.6), c(8.65, 9.22), TRUE)
        box(3.6, 9.3, 5.6, 3, col, c("MD GWAS", N(d$gwas$n[1]), "", paste0(dtdi_lab, " GWAS"), N(d$gwas$n[2]), 
            ""), c(7, 7, 6, 7, 7, 6), final = TRUE)
        line(c(0.98, 6.22), rep(10.78, 2))
        note(paste0("Core data missing: ", num(d$n[3] - d$n[4])), 10.2, 7.12)
        box(10.2, 7.4, 5.6, 1.25, col, c("Complete core covariates", "Age, sex, centre, batch, PC1-PC5,", "BMI, smoking and alcohol", 
            N(d$n[4])), c(6, 6, 6, 7))
        line(c(10.2, 10.2), c(8.65, 9.37), TRUE)
        note(paste0("Linkage / TDI / date incomplete: ", num(d$n[4] - d$n[5])), 10.2, 9.04)
        box(10.2, 9.45, 5.6, 1.25, col, c("Linked geographic and recruitment data", "Birthplace TDI and recruitment date", 
            N(d$n[5])), c(6, 6, 7))
        line(c(10.2, 10.2), c(10.7, 11.42), TRUE)
        note(paste0("Migration measures missing: ", num(d$n[5] - d$n[6])), 10.2, 11.08)
        box(10.2, 11.5, 5.6, 1.25, col, c("Complete migration measures", paste0("MD and ", dtdi_lab, " with model data"), 
            N(d$n[6])), c(6, 6, 7))
        line(c(10.2, 10.2), c(12.75, 13.2))
        line(c(3.6, 10.2), rep(13.2, 2))
        line(c(3.6, 3.6), c(13.2, 13.92), TRUE)
        line(c(10.2, 10.2), c(13.2, 13.92), TRUE)
        note("Outcome-specific exclusions", 3.6, 13.56)
        note(paste0("SES components missing: ", num(d$n[6] - d$ses)), 10.2, 13.56)
        box(3.6, 14, 5.6, 2.6, col, c("Health analyses", "10 outcomes", paste0("N = ", paste(num(d$health), collapse = "-")), 
            "Prevalent cases excluded for each disease", "Valid prospective follow-up required"), c(7, 6, 7, 6, 
            6), final = TRUE)
        box(10.2, 14, 5.6, 2.6, col, c("SES analysis", N(d$ses), "Education years, household income", "and CAMSIS complete"), 
            c(7, 7, 6, 6), final = TRUE)
        for (i in 1:2) {
            cohort <- c("UKHLS", "HRS")[i]
            col <- cols[cohort]
            d <- r$display[[cohort]]
            shift <- (i - 1) * 8.7
            cx <- 17.1
            txt(letters[i + 1], 14.25, 0.4 + shift, 10, face = "bold", just = "left")
            txt(cohort, 14.85, 0.42 + shift, 7, col, face = "bold", just = "left")
            top <- c(1, 2.62, 4.33) + shift
            heights <- c(0.95, 1, 0.95)
            for (j in 1:3) {
                label <- c(strsplit(d$labels[j], "\n", fixed = TRUE)[[1]], N(d$n[j]))
                box(cx, top[j], 5.2, heights[j], col, label, c(rep(6, length(label) - 1), 7))
                if (j < 3) 
                  line(c(cx, cx), c(top[j] + heights[j], top[j + 1] - 0.08), TRUE)
            }
            y <- c(2.13, 3.91) + shift
            for (j in 1:2) {
                line(c(cx, 20), rep(y[j], 2), TRUE)
                exclude(21.8, y[j] - 0.38, 3.4, c(if (j == 1) {
                  if (cohort == "UKHLS") "Age / ethnicity ineligible" else "Age ineligible or missing"
                } else "Migration history unavailable", paste0("Excluded: ", num(d$n[j] - d$n[j + 1]))))
            }
            bx <- c(16.6, 21.25)
            fork <- 5.55 + shift
            leaf <- 6.2 + shift
            line(c(cx, cx), c(top[3] + heights[3], fork))
            line(bx, rep(fork, 2))
            for (x in bx) line(c(x, x), c(fork, leaf - 0.08), TRUE)
            note(paste0("SES / covariates incomplete: ", num(d$n[3] - d$ses)), bx[1], 5.91 + shift)
            note(paste0("Follow-up / data incomplete: ", num(d$n[3] - d$mortality)), bx[2], 5.91 + shift)
            result(bx[1], leaf, 4.1, 1.96, col, "SES analysis", d$ses, d$ses_migrant, d$ses_non)
            result(bx[2], leaf, 4.1, 1.96, col, "All-cause mortality", d$mortality, d$mortality_migrant, d$mortality_non, 
                if (d$mortality_migrant != d$mortality_md) 
                  d$mortality_md
                else NULL)
        }
        grid::popViewport()
    }
}

draw_ppt(file.path(out, "ED1.pptx"), ed1_plot(plot_data$result), 24 * cm, 17.8 * cm)

finish("ED1")

