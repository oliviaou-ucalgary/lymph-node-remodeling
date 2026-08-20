# CCL21A joint clustering
#
# Mouse lymph-node cells/spots were selected manually in Loupe Browser using
# Ccl21a. These CSV files are external projection inputs; R does not recreate
# the manual selection.

# Import libraries and data ------------------------------------------------------
source(file.path("config", "analysis_config.R"))
source(file.path("R", "reproducibility.R"))
require_packages(c("Seurat", "dplyr", "ggplot2", "patchwork"))
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2); library(patchwork) })
set.seed(analysis_config$seed)
selection_gene <- "Ccl21a"

ccl21a_paths <- c(
  analysis_config$paths$wt_ctrl_h5, analysis_config$paths$mt_ctrl_h5,
  analysis_config$paths$wt_oxd4_h5, analysis_config$paths$mt_oxd4_h5,
  analysis_config$paths$wt_ctrl_ccl21a_barcodes,
  analysis_config$paths$mt_ctrl_ccl21a_barcodes,
  analysis_config$paths$wt_oxd4_ccl21a_barcodes,
  analysis_config$paths$mt_oxd4_ccl21a_barcodes
)
require_inputs(ccl21a_paths)

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
ccl21a_groups <- list(
  WT_CTRL = subset(ln_objects$WT_CTRL, cells = require_barcode_overlap(ln_objects$WT_CTRL, read_barcode_only_csv(analysis_config$paths$wt_ctrl_ccl21a_barcodes), "WT CTRL CCL21A selection")),
  MT_CTRL = subset(ln_objects$MT_CTRL, cells = require_barcode_overlap(ln_objects$MT_CTRL, read_barcode_only_csv(analysis_config$paths$mt_ctrl_ccl21a_barcodes), "MT CTRL CCL21A selection")),
  WT_OXD4 = subset(ln_objects$WT_OXD4, cells = require_barcode_overlap(ln_objects$WT_OXD4, read_barcode_only_csv(analysis_config$paths$wt_oxd4_ccl21a_barcodes), "WT OXD4 CCL21A selection")),
  MT_OXD4 = subset(ln_objects$MT_OXD4, cells = require_barcode_overlap(ln_objects$MT_OXD4, read_barcode_only_csv(analysis_config$paths$mt_oxd4_ccl21a_barcodes), "MT OXD4 CCL21A selection"))
)
for (group_name in names(ccl21a_groups)) ccl21a_groups[[group_name]]$group <- group_name

ccl21a_joint <- merge(
  ccl21a_groups$WT_CTRL,
  y = ccl21a_groups[c("MT_CTRL", "WT_OXD4", "MT_OXD4")],
  add.cell.ids = names(ccl21a_groups), project = "CCL21A_JOINT"
)
ccl21a_joint$genotype <- ifelse(grepl("^WT_", ccl21a_joint$group), "WT", "MT")
ccl21a_joint$condition <- ifelse(grepl("_CTRL$", ccl21a_joint$group), "CTRL", "OXD4")
DefaultAssay(ccl21a_joint) <- analysis_config$assay
ccl21a_joint <- NormalizeData(ccl21a_joint, scale.factor = analysis_config$normalization_scale_factor, verbose = FALSE)
ccl21a_joint <- FindVariableFeatures(ccl21a_joint, nfeatures = analysis_config$variable_features, verbose = FALSE)
ccl21a_joint <- ScaleData(ccl21a_joint, verbose = FALSE)
ccl21a_joint <- RunPCA(ccl21a_joint, npcs = max(analysis_config$pca_dims), verbose = FALSE)
ccl21a_joint <- FindNeighbors(ccl21a_joint, dims = analysis_config$pca_dims, verbose = FALSE)
ccl21a_joint <- FindClusters(ccl21a_joint, resolution = analysis_config$clustering_resolution, verbose = FALSE)
ccl21a_joint <- RunUMAP(ccl21a_joint, dims = analysis_config$pca_dims,
                        n.neighbors = analysis_config$umap_neighbors,
                        min.dist = analysis_config$umap_min_dist,
                        seed.use = analysis_config$seed, verbose = FALSE)

print(table(ccl21a_joint$group, ccl21a_joint$seurat_clusters))
ccl21a_umap <- DimPlot(ccl21a_joint, group.by = "group") |
  DimPlot(ccl21a_joint, group.by = "seurat_clusters", label = TRUE, repel = TRUE)
print(ccl21a_umap)
