# Optional CXCL13 gene-expression plots
# Add genes only after opening the joint object in RStudio.

# Import libraries and data ------------------------------------------------------
source(file.path("config", "analysis_config.R"))
source(file.path("R", "reproducibility.R"))
require_packages(c("Seurat", "ggplot2", "patchwork"))
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })
if (!exists("cxcl13_joint") || !inherits(cxcl13_joint, "Seurat")) stop("Run scripts/01_cxcl13_joint_clustering.R first.", call. = FALSE)
validate_seurat(cxcl13_joint, assay = "RNA", metadata = c("group", "seurat_clusters"))

# Genes previously used included: Cxcl13, Ccl21a, and Ccl19.
genes_to_plot <- c(
  # Add mouse gene symbols here.
)

# Actual analysis ---------------------------------------------------------------
if (!length(genes_to_plot)) {
  message("No optional genes specified; edit genes_to_plot to display plots.")
} else {
  missing_genes <- setdiff(genes_to_plot, rownames(cxcl13_joint))
  if (length(missing_genes)) stop("Genes not found: ", paste(missing_genes, collapse = ", "), call. = FALSE)
  cxcl13_feature_plot <- FeaturePlot(cxcl13_joint, features = genes_to_plot, split.by = "group")
  cxcl13_violin_plot <- VlnPlot(cxcl13_joint, features = genes_to_plot, group.by = "seurat_clusters", split.by = "group", pt.size = 0)
  print(cxcl13_feature_plot)
  print(cxcl13_violin_plot)
}

