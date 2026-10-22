# 1. Approved public point-map rasters
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))

use("geo")

point_layer <- function(d, file, dpi = 600) {
    ragg::agg_png(file, width = 5.69, height = 9.35, units = "cm", res = dpi, background = "transparent")
    on.exit(dev.off())
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(xscale = c(-45000, 710000), yscale = c(-20000, 1220000), clip = "on"))
    d <- unique(d)
    angle <- seq(0, 2 * pi, length.out = 9)[-9]
    radius <- 0.0075 * 755000/5.69
    for (start in seq.int(1L, nrow(d), by = 5000L)) {
        p <- d[start:min(nrow(d), start + 4999L)]
        xx <- as.vector(t(outer(p$easting, radius * cos(angle), "+")))
        yy <- as.vector(t(outer(p$northing, radius * sin(angle), "+")))
        grid::grid.path(grid::unit(xx, "native"), grid::unit(yy, "native"), id.lengths = rep(8L, nrow(p)), rule = "winding", 
            gp = grid::gpar(fill = "#01A6AC", col = NA))
    }
}

ed2_export <- function(destination = Sys.getenv("MOVE_SOURCE", file.path(root, "source"))) {
    dir.create(destination, recursive = TRUE, showWarnings = FALSE)
    x <- ed2_data()
    layers <- c(birth = "ED2_birth.png", home = "ED2_home.png")
    for (n in names(layers)) point_layer(x[[n]], file.path(destination, layers[n]))
    x[c("birth", "home")] <- NULL
    x$layers <- layers
    source_write(list(result = x), file.path(destination, "ED2_source_data.xlsx"))
}

if (sys.nframe() == 0L) ed2_export()

