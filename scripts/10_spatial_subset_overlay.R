# Interactive CXCL13 and CCL21A subset overlays
#
# Loupe Browser exports contain only barcodes. Each subset is joined directly to
# the corresponding slide projection, whose required columns are Barcode,
# X Coordinate, and Y Coordinate. The projection coordinates are expected to
# match the tissue image; optional alignment values default to no transformation.

# Import libraries and data ------------------------------------------------------
source(file.path("config", "analysis_config.R"))
source(file.path("R", "reproducibility.R"))
require_packages(c("dplyr", "ggplot2", "png", "jsonlite", "patchwork"))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(png)
  library(jsonlite)
  library(patchwork)
})

overlay_inputs <- c(
  analysis_config$paths$ctrl_projection,
  analysis_config$paths$oxd4_projection,
  analysis_config$paths$ctrl_hires_png,
  analysis_config$paths$oxd4_hires_png,
  analysis_config$paths$ctrl_scale_json,
  analysis_config$paths$oxd4_scale_json,
  analysis_config$paths$wt_ctrl_barcodes,
  analysis_config$paths$mt_ctrl_barcodes,
  analysis_config$paths$wt_oxd4_barcodes,
  analysis_config$paths$mt_oxd4_barcodes,
  analysis_config$paths$wt_ctrl_ccl21a_barcodes,
  analysis_config$paths$mt_ctrl_ccl21a_barcodes,
  analysis_config$paths$wt_oxd4_ccl21a_barcodes,
  analysis_config$paths$mt_oxd4_ccl21a_barcodes
)
require_inputs(overlay_inputs)

slide_files <- list(
  CTRL = list(
    projection = read_projection_csv(analysis_config$paths$ctrl_projection),
    image = png::readPNG(analysis_config$paths$ctrl_hires_png),
    scale_factors = jsonlite::fromJSON(analysis_config$paths$ctrl_scale_json)
  ),
  OXD4 = list(
    projection = read_projection_csv(analysis_config$paths$oxd4_projection),
    image = png::readPNG(analysis_config$paths$oxd4_hires_png),
    scale_factors = jsonlite::fromJSON(analysis_config$paths$oxd4_scale_json)
  )
)

subset_manifest <- data.frame(
  gene = c(rep("Cxcl13", 4), rep("Ccl21a", 4)),
  group = rep(c("WT_CTRL", "MT_CTRL", "WT_OXD4", "MT_OXD4"), 2),
  slide = rep(c("CTRL", "CTRL", "OXD4", "OXD4"), 2),
  path = c(
    analysis_config$paths$wt_ctrl_barcodes,
    analysis_config$paths$mt_ctrl_barcodes,
    analysis_config$paths$wt_oxd4_barcodes,
    analysis_config$paths$mt_oxd4_barcodes,
    analysis_config$paths$wt_ctrl_ccl21a_barcodes,
    analysis_config$paths$mt_ctrl_ccl21a_barcodes,
    analysis_config$paths$wt_oxd4_ccl21a_barcodes,
    analysis_config$paths$mt_oxd4_ccl21a_barcodes
  ),
  stringsAsFactors = FALSE
)

# OPTIONAL MANUAL ALIGNMENT ------------------------------------------------------
# Leave these values unchanged when projection coordinates match the image.
# Positive x_shift moves right; positive y_shift moves down. Adjust only after
# comparing the preview to known tissue landmarks.
overlay_alignment <- list(
  WT_CTRL = list(swap_xy = FALSE, flip_x = FALSE, flip_y = FALSE, x_shift = 0, y_shift = 0, x_scale = 1, y_scale = 1),
  MT_CTRL = list(swap_xy = FALSE, flip_x = FALSE, flip_y = FALSE, x_shift = 0, y_shift = 0, x_scale = 1, y_scale = 1),
  WT_OXD4 = list(swap_xy = FALSE, flip_x = FALSE, flip_y = FALSE, x_shift = 0, y_shift = 0, x_scale = 1, y_scale = 1),
  MT_OXD4 = list(swap_xy = FALSE, flip_x = FALSE, flip_y = FALSE, x_shift = 0, y_shift = 0, x_scale = 1, y_scale = 1)
)

transform_coordinates <- function(x, image, settings) {
  image_width <- dim(image)[2]
  image_height <- dim(image)[1]
  x$x_plot <- x[["X Coordinate"]]
  x$y_plot <- x[["Y Coordinate"]]
  if (settings$swap_xy) {
    temporary <- x$x_plot
    x$x_plot <- x$y_plot
    x$y_plot <- temporary
  }
  x$x_plot <- x$x_plot * settings$x_scale
  x$y_plot <- x$y_plot * settings$y_scale
  if (settings$flip_x) x$x_plot <- image_width - x$x_plot
  if (settings$flip_y) x$y_plot <- image_height - x$y_plot
  x$x_plot <- x$x_plot + settings$x_shift
  x$y_plot <- x$y_plot + settings$y_shift
  x$inside_image <- x$x_plot >= 0 & x$x_plot <= image_width & x$y_plot >= 0 & x$y_plot <= image_height
  x
}

make_overlay <- function(coordinates, image, title) {
  image_width <- dim(image)[2]
  image_height <- dim(image)[1]
  plot_data <- coordinates[coordinates$matched & coordinates$inside_image, , drop = FALSE]
  ggplot() +
    annotation_raster(image, xmin = 0, xmax = image_width, ymin = 0, ymax = image_height) +
    geom_point(data = plot_data, aes(x_plot, y_plot), color = "#D55E00", size = 0.35, alpha = 0.8) +
    coord_fixed(xlim = c(0, image_width), ylim = c(image_height, 0), expand = FALSE) +
    theme_void() +
    labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

# Actual analysis ---------------------------------------------------------------
overlay_coordinates <- list(Cxcl13 = list(), Ccl21a = list())
overlay_plots <- list(Cxcl13 = list(), Ccl21a = list())

for (row_index in seq_len(nrow(subset_manifest))) {
  item <- subset_manifest[row_index, ]
  slide <- slide_files[[item$slide]]
  barcodes <- read_barcode_only_csv(item$path)
  coordinates <- match_projection_coordinates(barcodes, slide$projection, paste(item$gene, item$group))
  coordinates$matched <- stats::complete.cases(coordinates[, c("X Coordinate", "Y Coordinate")])
  coordinates$gene <- item$gene
  coordinates$group <- item$group
  coordinates$slide <- item$slide
  coordinates <- transform_coordinates(coordinates, slide$image, overlay_alignment[[item$group]])
  message(item$gene, " ", item$group, ": ", sum(coordinates$inside_image & coordinates$matched),
          " matched barcode(s) inside the image.")
  overlay_coordinates[[item$gene]][[item$group]] <- coordinates
  overlay_plots[[item$gene]][[item$group]] <- make_overlay(
    coordinates, slide$image, paste(item$gene, item$group)
  )
}

cxcl13_overlay_panel <- wrap_plots(overlay_plots$Cxcl13, ncol = 2)
ccl21a_overlay_panel <- wrap_plots(overlay_plots$Ccl21a, ncol = 2)
print(cxcl13_overlay_panel)
print(ccl21a_overlay_panel)

