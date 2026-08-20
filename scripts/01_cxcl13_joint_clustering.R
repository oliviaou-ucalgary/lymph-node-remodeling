# CXCL13 joint clustering
#
# Mouse lymph-node cells/spots were selected manually in Loupe Browser using
# Cxcl13. The exported barcode CSV files are inputs to this projection; this
# script does not recreate the manual Loupe selection.

# Import libraries and data ------------------------------------------------------
source(file.path("config", "analysis_config.R"))
source(file.path("R", "reproducibility.R"))
require_packages(c("Seurat", "dplyr", "ggplot2", "patchwork"))
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})
set.seed(analysis_config$seed)
selection_gene <- "Cxcl13"

input_paths <- c(
  analysis_config$paths$wt_ctrl_h5,
  analysis_config$paths$mt_ctrl_h5,
  analysis_config$paths$wt_oxd4_h5,
  analysis_config$paths$mt_oxd4_h5,
  analysis_config$paths$wt_ctrl_barcodes,
  analysis_config$paths$mt_ctrl_barcodes,
  analysis_config$paths$wt_oxd4_barcodes,
  analysis_config$paths$mt_oxd4_barcodes
)
require_inputs(input_paths)

read_ln_object <- function(path, project) {
  counts <- select_gene_expression(Read10X_h5(path), path)
  CreateSeuratObject(counts, project = project)
}
ln_objects <- list(
  WT_CTRL = read_ln_object(analysis_config$paths$wt_ctrl_h5, "WT_CTRL"),
  MT_CTRL = read_ln_object(analysis_config$paths$mt_ctrl_h5, "MT_CTRL"),
  WT_OXD4 = read_ln_object(analysis_config$paths$wt_oxd4_h5, "WT_OXD4"),
  MT_OXD4 = read_ln_object(analysis_config$paths$mt_oxd4_h5, "MT_OXD4")
)
if (!all(vapply(ln_objects, function(x) selection_gene %in% rownames(x), logical(1)))) {
  stop("Selection gene is absent from one or more inputs: ", selection_gene, call. = FALSE)
}

# Actual analysis ---------------------------------------------------------------
cxcl13_groups <- list(
  WT_CTRL = subset(ln_objects$WT_CTRL, cells = require_barcode_overlap(ln_objects$WT_CTRL, read_barcode_only_csv(analysis_config$paths$wt_ctrl_barcodes), "WT CTRL CXCL13 selection")),
  MT_CTRL = subset(ln_objects$MT_CTRL, cells = require_barcode_overlap(ln_objects$MT_CTRL, read_barcode_only_csv(analysis_config$paths$mt_ctrl_barcodes), "MT CTRL CXCL13 selection")),
  WT_OXD4 = subset(ln_objects$WT_OXD4, cells = require_barcode_overlap(ln_objects$WT_OXD4, read_barcode_only_csv(analysis_config$paths$wt_oxd4_barcodes), "WT OXD4 CXCL13 selection")),
  MT_OXD4 = subset(ln_objects$MT_OXD4, cells = require_barcode_overlap(ln_objects$MT_OXD4, read_barcode_only_csv(analysis_config$paths$mt_oxd4_barcodes), "MT OXD4 CXCL13 selection"))
)

for (group_name in names(cxcl13_groups)) cxcl13_groups[[group_name]]$group <- group_name

cxcl13_joint <- merge(
  cxcl13_groups$WT_CTRL,
  y = cxcl13_groups[c("MT_CTRL", "WT_OXD4", "MT_OXD4")],
  add.cell.ids = names(cxcl13_groups),
  project = "CXCL13_JOINT"
)
cxcl13_joint$genotype <- ifelse(grepl("^WT_", cxcl13_joint$group), "WT", "MT")
cxcl13_joint$condition <- ifelse(grepl("_CTRL$", cxcl13_joint$group), "CTRL", "OXD4")
DefaultAssay(cxcl13_joint) <- analysis_config$assay
cxcl13_joint <- NormalizeData(cxcl13_joint, scale.factor = analysis_config$normalization_scale_factor, verbose = FALSE)
cxcl13_joint <- FindVariableFeatures(cxcl13_joint, nfeatures = analysis_config$variable_features, verbose = FALSE)
cxcl13_joint <- ScaleData(cxcl13_joint, verbose = FALSE)
cxcl13_joint <- RunPCA(cxcl13_joint, npcs = max(analysis_config$pca_dims), verbose = FALSE)
cxcl13_joint <- FindNeighbors(cxcl13_joint, dims = analysis_config$pca_dims, verbose = FALSE)
cxcl13_joint <- FindClusters(cxcl13_joint, resolution = analysis_config$clustering_resolution, verbose = FALSE)
cxcl13_joint <- RunUMAP(cxcl13_joint, dims = analysis_config$pca_dims,
                        n.neighbors = analysis_config$umap_neighbors,
                        min.dist = analysis_config$umap_min_dist,
                        seed.use = analysis_config$seed, verbose = FALSE)

print(table(cxcl13_joint$group, cxcl13_joint$seurat_clusters))
cxcl13_umap <- DimPlot(cxcl13_joint, group.by = "group") |
  DimPlot(cxcl13_joint, group.by = "seurat_clusters", label = TRUE, repel = TRUE)
print(cxcl13_umap)
