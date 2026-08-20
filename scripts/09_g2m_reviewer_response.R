# G2M analysis for the reviewer response
#
# WT CTRL LN1 and WT OXD4 LN2 are scored separately. The H5 barcodes are joined
# to the corresponding whole-slide projection for interactive tissue overlays.
# Spatial bins can contain mixtures of cells, so these scores are spatial
# expression summaries rather than proof of a single-cell state.

# Import libraries and data ------------------------------------------------------
source(file.path("config", "analysis_config.R"))
source(file.path("R", "reproducibility.R"))
require_packages(c("Seurat", "ggplot2", "png", "jsonlite", "patchwork"))
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(png); library(jsonlite); library(patchwork) })
set.seed(analysis_config$seed)
require_inputs(c(
  analysis_config$paths$wt_ctrl_h5,
  analysis_config$paths$wt_oxd4_h5,
  analysis_config$paths$ctrl_projection,
  analysis_config$paths$oxd4_projection,
  analysis_config$paths$ctrl_hires_png,
  analysis_config$paths$oxd4_hires_png,
  analysis_config$paths$ctrl_scale_json,
  analysis_config$paths$oxd4_scale_json
))

g2m_inputs <- list(
  WT_CTRL = list(h5 = analysis_config$paths$wt_ctrl_h5,
                 projection = analysis_config$paths$ctrl_projection,
                 image = analysis_config$paths$ctrl_hires_png,
                 scale_factors = analysis_config$paths$ctrl_scale_json),
  WT_OXD4 = list(h5 = analysis_config$paths$wt_oxd4_h5,
                 projection = analysis_config$paths$oxd4_projection,
                 image = analysis_config$paths$oxd4_hires_png,
                 scale_factors = analysis_config$paths$oxd4_scale_json)
)

# OPTIONAL MANUAL ALIGNMENT ------------------------------------------------------
# Projection coordinates should match the images. Leave unchanged unless a
# landmark-based preview demonstrates that adjustment is required.
g2m_alignment <- list(
  WT_CTRL = list(swap_xy = FALSE, flip_x = FALSE, flip_y = FALSE, x_shift = 0, y_shift = 0, x_scale = 1, y_scale = 1),
  WT_OXD4 = list(swap_xy = FALSE, flip_x = FALSE, flip_y = FALSE, x_shift = 0, y_shift = 0, x_scale = 1, y_scale = 1)
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

make_score_overlay <- function(data, image, feature, title) {
  image_width <- dim(image)[2]
  image_height <- dim(image)[1]
  plot_data <- data[data$matched & data$inside_image & is.finite(data[[feature]]), , drop = FALSE]
  ggplot() +
    annotation_raster(image, xmin = 0, xmax = image_width, ymin = 0, ymax = image_height) +
    geom_point(data = plot_data, aes(x = x_plot, y = y_plot, color = .data[[feature]]), size = 0.35, alpha = 0.8) +
    scale_color_viridis_c() +
    coord_fixed(xlim = c(0, image_width), ylim = c(image_height, 0), expand = FALSE) +
    theme_void() + labs(title = title, color = feature) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

# Actual analysis ---------------------------------------------------------------
g2m_objects <- list()
g2m_coordinates <- list()
g2m_plots <- list()

for (group_name in names(g2m_inputs)) {
  input <- g2m_inputs[[group_name]]
  scale_factors <- jsonlite::fromJSON(input$scale_factors)
  message(group_name, " scale-factor keys: ", paste(names(scale_factors), collapse = ", "))
  counts <- select_gene_expression(Read10X_h5(input$h5), input$h5)
  object <- CreateSeuratObject(counts, project = group_name)
  DefaultAssay(object) <- "RNA"
  object <- NormalizeData(object, scale.factor = analysis_config$normalization_scale_factor, verbose = FALSE)

  s_genes_mouse <- unique(stats::na.omit(CaseMatch(search = cc.genes.updated.2019$s.genes, match = rownames(object))))
  g2m_genes_mouse <- unique(stats::na.omit(CaseMatch(search = cc.genes.updated.2019$g2m.genes, match = rownames(object))))
  message(group_name, ": matched ", length(s_genes_mouse), " S genes and ", length(g2m_genes_mouse), " G2M genes.")
  if (length(s_genes_mouse) < 10L || length(g2m_genes_mouse) < 10L) {
    stop("Too few mouse cell-cycle genes matched in ", group_name, ".", call. = FALSE)
  }
  object <- CellCycleScoring(object, s.features = s_genes_mouse,
                             g2m.features = g2m_genes_mouse, set.ident = FALSE)
  print(table(group = group_name, phase = object$Phase, useNA = "ifany"))

  projection <- read_projection_csv(input$projection)
  coordinates <- match_projection_coordinates(colnames(object), projection, group_name)
  coordinates$matched <- stats::complete.cases(coordinates[, c("X Coordinate", "Y Coordinate")])
  expression <- FetchData(object, vars = intersect(c("Mki67", "Top2a", "G2M.Score"), c(rownames(object), colnames(object[[]]))))
  expression$Barcode <- rownames(expression)
  coordinates <- merge(coordinates, expression, by = "Barcode", all.x = TRUE, sort = FALSE)
  image <- png::readPNG(input$image)
  coordinates <- transform_coordinates(coordinates, image, g2m_alignment[[group_name]])

  features <- intersect(c("Mki67", "Top2a", "G2M.Score"), colnames(coordinates))
  g2m_plots[[group_name]] <- lapply(features, function(feature) {
    make_score_overlay(coordinates, image, feature, paste(group_name, feature))
  })
  names(g2m_plots[[group_name]]) <- features
  g2m_objects[[group_name]] <- object
  g2m_coordinates[[group_name]] <- coordinates
}

g2m_reviewer_panel <- wrap_plots(c(g2m_plots$WT_CTRL, g2m_plots$WT_OXD4), ncol = 3)
print(g2m_reviewer_panel)
