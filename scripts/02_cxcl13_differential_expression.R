# CXCL13 cluster-wise differential expression
#
# This preserves the original exploratory cell/spot-level Wilcoxon comparison.
# Cells/spots are not biological replicates; do not interpret the p-values as a
# replicate-aware genotype effect.

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
if (!exists("cxcl13_joint") || !inherits(cxcl13_joint, "Seurat")) {
  stop("Run scripts/01_cxcl13_joint_clustering.R first in this RStudio session.", call. = FALSE)
}
validate_seurat(cxcl13_joint, assay = "RNA", metadata = c("group", "seurat_clusters"))

# Actual analysis ---------------------------------------------------------------
cxcl13_de_object <- cxcl13_joint
cxcl13_de_object[["RNA"]] <- JoinLayers(cxcl13_de_object[["RNA"]])
Idents(cxcl13_de_object) <- "seurat_clusters"
cxcl13_de_results <- list()
cxcl13_volcano_plots <- list()

for (cluster_id in levels(Idents(cxcl13_de_object))) {
  group_counts <- table(cxcl13_de_object$group[cxcl13_de_object$seurat_clusters == cluster_id])
  comparison_counts <- group_counts[c("WT_CTRL", "MT_CTRL")]
  comparison_counts[is.na(comparison_counts)] <- 0
  if (any(comparison_counts < analysis_config$min_cells_per_group)) next

  result <- FindMarkers(
    cxcl13_de_object,
    subset.ident = cluster_id,
    group.by = "group",
    ident.1 = "WT_CTRL",
    ident.2 = "MT_CTRL",
    logfc.threshold = analysis_config$logfc_threshold,
    min.pct = analysis_config$min_pct,
    test.use = "wilcox"
  )
  result$gene <- rownames(result)
  result$cluster <- cluster_id
  cxcl13_de_results[[cluster_id]] <- result

  result$category <- "Not significant"
  result$category[result$p_val_adj < analysis_config$adjusted_p_cutoff & result$avg_log2FC > analysis_config$lfc_cutoff] <- "WT higher"
  result$category[result$p_val_adj < analysis_config$adjusted_p_cutoff & result$avg_log2FC < -analysis_config$lfc_cutoff] <- "MT higher"
  cxcl13_volcano_plots[[cluster_id]] <- ggplot(result, aes(avg_log2FC, -log10(pmax(p_val_adj, .Machine$double.xmin)), color = category)) +
    geom_point(alpha = 0.6, size = 1) +
    theme_classic() +
    labs(title = paste("Cluster", cluster_id, "WT CTRL vs MT CTRL"), x = "Average log2 fold change", y = "-log10 adjusted p-value")
}

print(cxcl13_de_results)
if (length(cxcl13_volcano_plots)) print(wrap_plots(cxcl13_volcano_plots))

