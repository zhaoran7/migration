# 1. Paths and module loading
if (!exists('.move_ready', inherits=FALSE)) {
  .move_ready <- TRUE
  suppressPackageStartupMessages(library(data.table))
  code <- Sys.getenv('MOVE_CODE',dirname(normalizePath(sub('^--file=', '',grep('^--file=',commandArgs(),value=TRUE)[1]))))
  root <- Sys.getenv('MOVE_ROOT', getwd())
  Sys.setenv(MOVE_ROOT=root)
  phe <- Sys.getenv('MOVE_PHE', file.path(root,'data/phe'))
  ukb <- Sys.getenv('MOVE_UKB',file.path(root,'data/ukb'))
  out <- Sys.getenv('MOVE_RESULTS',file.path(root,'results'))
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  core <- as.integer(Sys.getenv('MOVE_CORES','1'))
  testn <- 0L
  Sys.setenv(OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1)
  proj <- file.path(dirname(dirname(R.home())),'share/proj')
  if(!nzchar(Sys.getenv('PROJ_DATA'))&&file.exists(file.path(proj,'proj.db')))
    Sys.setenv(PROJ_DATA=proj,PROJ_LIB=proj)
  modules <- fread(file.path(code,'modules.tsv'))
  source(file.path(code,'00_common/02_definitions.R'))
  source(file.path(code,'00_common/03_source.R'))
  .loaded <- character()
  use <- function(names) for(n in names) if(!n%in%.loaded) {
    .loaded <<- c(.loaded,n)
    packages <- switch(n,ukb='haven',ukhls='haven',hrs='haven',
      assoc=c('survival','parallel'),region=c('coxme','survival'),
      external=c('survival','lme4'),omics=c('survival','parallel'),
      sensitivity=c('survival','parallel'),ses='lme4',character())
    for(p in packages) suppressPackageStartupMessages(library(p,character.only=TRUE))
    path <- modules$Script[match(n,modules$Module)]
    if(is.na(path)) stop('Unknown module: ',n)
    source(file.path(code,path),local=.GlobalEnv)
  }
}

rint <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    y <- rep(NA_real_, length(x))
    i <- is.finite(x)
    if (sum(i) >= 3) 
        y[i] <- qnorm((rank(x[i], ties.method = "average") - 0.5)/sum(i))
    y
}

z <- function(x) as.numeric(scale(suppressWarnings(as.numeric(x))))

ed_out <- function(x) {
  a <- c('CAD/MI'='CAD','Type 2 diabetes'='T2D','Anxiety disorder'='Anxiety')
  ifelse(x%in%names(a),a[x],x)
}

fmt <- function(x) {
    x <- toupper(trimws(as.character(x)))
    x <- gsub("\\.", "", x)
    x[x %in% c("", "NA", "NAN")] <- NA_character_
    x
}

hit <- function(x, p) Reduce(`|`, lapply(p, function(z) startsWith(fmt(x), z)))

first0 <- function(x) {
    x <- x[!is.na(x)]
    if (length(x)) 
        as.numeric(x[1])
    else NA_real_
}

last0 <- function(x) {
    x <- x[!is.na(x)]
    if (length(x)) 
        as.numeric(x[length(x)])
    else NA_real_
}

read_ids <- function(x) {
    f <- file.path(phe, paste0(x, ".ids"))
    if (file.exists(f)) 
        unique(as.character(fread(f, header = FALSE)[[2]]))
    else character()
}

stale <- function(file, inputs) !file.exists(file) || any(file.info(inputs)$mtime > file.info(file)$mtime, 
    na.rm = TRUE)

ukb_inputs <- function() c(file.path(phe, "move.rds"), if (file.exists(file.path(phe, "ttwa.rds"))) file.path(phe, 
    "ttwa.rds") else file.path(phe, "ukb.pheno"), vapply(c("death", "icd10.code0", "icd10.date0"), ukb_input, 
    character(1)))

ukb_input <- function(name) {
    local <- file.path(phe, "ukb", paste0(name, ".rds"))
    if (file.exists(local)) 
        local
    else file.path(ukb, paste0(name, ".rds"))
}

xlsx_clean <- function(file) {
    tmp <- tempfile("xlsx_")
    dir.create(tmp)
    unzip(file, exdir = tmp)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    f <- file.path(tmp, "xl/workbook.xml")
    d <- xml2::read_xml(f)
    r <- xml2::xml_find_all(d, "//*[local-name()='workbookView']")
    for (node in r) for (key in c("xWindow", "yWindow", "windowWidth", "windowHeight", "activeTab", "firstSheet", 
        "tabRatio")) if (identical(xml2::xml_attr(node, key), "")) 
        xml2::xml_attr(node, key) <- NULL
    xml2::write_xml(d, f)
    f <- file.path(tmp, "xl/styles.xml")
    d <- xml2::read_xml(f)
    r <- xml2::xml_find_all(d, "//*[@rgb]")
    color <- sub("^#", "", xml2::xml_attr(r, "rgb"))
    xml2::xml_attr(r, "rgb") <- ifelse(nchar(color) == 6, paste0("FF", color), color)
    xml2::write_xml(d, f)
    for (f in list.files(tmp, "\\.rels$", recursive = TRUE, full.names = TRUE)) {
        d <- xml2::read_xml(f)
        r <- xml2::xml_children(d)
        target <- xml2::xml_attr(r, "Target")
        mode <- xml2::xml_attr(r, "TargetMode")
        paths <- ifelse(startsWith(target, "/"), file.path(tmp, substring(target, 2)), file.path(dirname(dirname(f)), 
            target))
        bad <- which(!file.exists(paths) & (is.na(mode) | mode != "External"))
        if (!length(bad)) 
            next
        owner <- file.path(dirname(dirname(f)), sub("\\.rels$", "", basename(f)))
        if (file.exists(owner)) {
            sheet <- xml2::read_xml(owner)
            for (id in xml2::xml_attr(r[bad], "Id")) {
                nodes <- xml2::xml_find_all(sheet, paste0("//*[@r:id='", id, "']"), c(r = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"))
                for (node in nodes) if (xml2::xml_name(node) %in% c("drawing", "legacyDrawing")) 
                  xml2::xml_remove(node)
                else xml2::xml_attr(node, "r:id") <- NULL
            }
            xml2::write_xml(sheet, owner)
        }
        xml2::xml_remove(r[bad])
        xml2::write_xml(d, f)
    }
    f <- file.path(tmp, "[Content_Types].xml")
    d <- xml2::read_xml(f)
    r <- xml2::xml_find_all(d, "//*[local-name()='Override']")
    xml2::xml_remove(r[!file.exists(file.path(tmp, substring(xml2::xml_attr(r, "PartName"), 2)))])
    xml2::write_xml(d, f)
    zip::zipr(file, list.files(tmp, recursive = TRUE, all.files = TRUE, no.. = TRUE), root = tmp, mode = "mirror")
}

save_xlsx <- function(x, file) {
    book <- openxlsx::createWorkbook(creator = "MOVE")
    openxlsx::modifyBaseFont(book, fontName = "Arial", fontSize = 10, fontColour = "#000000")
    header <- openxlsx::createStyle(fontName = "Arial", fontSize = 10, textDecoration = "bold", fgFill = "#F2F2F2", 
        border = "bottom", borderColour = "#808080", wrapText = TRUE, valign = "center")
    for (sheet in names(x)) {
        d <- as.data.frame(x[[sheet]])
        stopifnot(!anyDuplicated(names(d)), ncol(d) > 0, nrow(d) < 1048576)
        openxlsx::addWorksheet(book, sheet, gridLines = FALSE)
        openxlsx::writeData(book, sheet, d, headerStyle = header, withFilter = nrow(d) > 0, keepNA = FALSE)
        openxlsx::freezePane(book, sheet, firstRow = TRUE)
        openxlsx::setRowHeights(book, sheet, 1, 30)
        width <- pmin(38, pmax(12, nchar(names(d)) + 2))
        for (j in seq_along(d)) {
            if (is.character(d[[j]])) 
                width[j] <- min(38, max(width[j], suppressWarnings(max(nchar(d[[j]]), na.rm = TRUE)) + 
                  2))
            if (nrow(d) && is.numeric(d[[j]])) {
                is_p <- grepl("(^|[ _])(p|q)( |$)|FDR", names(d)[j])
                integer0 <- all(is.na(d[[j]]) | d[[j]] == trunc(d[[j]])) && !is_p
                st <- openxlsx::createStyle(numFmt = if (is_p) 
                  "0.00E+00"
                else if (integer0) 
                  "0"
                else "0.0000")
                openxlsx::addStyle(book, sheet, st, rows = 2:(nrow(d) + 1), cols = j, stack = TRUE)
            }
            if (grepl("(^|[ _])(p|q)( |$)|^Beta", names(d)[j])) 
                openxlsx::addStyle(book, sheet, openxlsx::createStyle(textDecoration = c("bold", "italic")), 
                  rows = 1, cols = j, stack = TRUE)
        }
        openxlsx::setColWidths(book, sheet, seq_along(d), width)
    }
    stage <- tempfile(fileext = ".xlsx")
    on.exit(unlink(stage), add = TRUE)
    openxlsx::saveWorkbook(book, stage, overwrite = TRUE)
    xlsx_clean(stage)
    if (!file.copy(stage, file, overwrite = TRUE)) 
        stop("Close the workbook before replacing: ", file)
    invisible(file)
}

source_data <- function(figure, panels) {
    labels <- c(cohort = "Cohort", outcome = "Outcome", exposure = "Exposure", term = "Term", model = "Model", 
        type = "Model", group = "Group", level = "Level", strat = "Stratum", n = "N", N = "N", event = "Events", 
        events = "Events", estimate = "Estimate", fit = "Estimate", beta = "Beta", se = "SE", HR = "HR", 
        ratio = "RHR", lo = "95% CI lower", hi = "95% CI upper", lower = "95% CI lower", upper = "95% CI upper", 
        p = "p", q = "q", FDR = "q (BH)", p_current = "p current", p_reference = "p reference", p_overall = "p overall", 
        p_nonlinear = "p non-linear", p_interaction = "p interaction", interaction_FDR = "q interaction", 
        interaction = "Interaction", df = "Degrees of freedom", x = "X", y = "Y", wave = "Wave", count = "MT", 
        mover = "Migrant", md = "MD (km)", val = "Correlation", half = "Correlation type", h2 = "SNP heritability", 
        h2_se = "SE", context = "Location", code = "Area code", area = "Area", MHR = "MHR", conditional_lower = "95% conditional CI lower", 
        conditional_upper = "95% conditional CI upper", p_approx = "p approximate", q_approx = "q approximate", 
        panel_q = "q regional variance")
    format_table <- function(d) {
        d <- copy(as.data.table(d))
        drop <- intersect(names(d), c("eid", "IID", "FID", "pid", "pidp", "HHIDPN"))
        if (length(drop)) 
            stop("Source data must not contain participant identifiers: ", paste(drop, collapse = ", "))
        for (v in names(d)) {
            if (is.factor(d[[v]])) 
                set(d, j = v, value = as.character(d[[v]]))
            if (is.list(d[[v]])) 
                stop("Nested source column: ", v)
        }
        old <- names(d)
        new <- unname(labels[old])
        new[is.na(new)] <- tools::toTitleCase(gsub("_", " ", old[is.na(new)]))
        i <- grepl("^[pq]($|_)", old)
        new[i] <- paste0(substr(old[i], 1, 1), substring(new[i], 2))
        setnames(d, make.unique(new, sep = " "))
        for (v in names(d)) if (is.character(d[[v]])) 
            d[get(v) == "Per SD PRS", `:=`((v), "Per SD PGS")]
        if ("Outcome" %in% names(d)) 
            d[Outcome == "CAD/MI", `:=`(Outcome, "CAD")]
        d
    }
    tables <- lapply(panels, function(x) {
        if (is.data.frame(x)) 
            return(format_table(x))
        rbindlist(lapply(names(x), function(nm) {
            d <- format_table(x[[nm]])
            d[, `:=`(Analysis, nm)]
            setcolorder(d, c("Analysis", setdiff(names(d), "Analysis")))
            d
        }), fill = TRUE)
    })
    save_xlsx(tables, file.path(out, paste0(figure, "_source_data.xlsx")))
}

