# CCL21A cluster-wise differential expression
# Preserves the original exploratory cell/spot-level Wilcoxon comparison.

# Import libraries and data ------------------------------------------------------
source(file.path("config", "analysis_config.R"))
source(file.path("R", "reproducibility.R"))
require_packages(c("Seurat", "dplyr", "ggplot2", "patchwork"))
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2); library(patchwork) })
if (!exists("ccl21a_joint") || !inherits(ccl21a_joint, "Seurat")) stop("Run scripts/05_ccl21a_joint_clustering.R first.", call. = FALSE)
validate_seurat(ccl21a_joint, assay = "RNA", metadata = c("group", "seurat_clusters"))

# Actual analysis ---------------------------------------------------------------
ccl21a_de_object <- ccl21a_joint
ccl21a_de_object[["RNA"]] <- JoinLayers(ccl21a_de_object[["RNA"]])
Idents(ccl21a_de_object) <- "seurat_clusters"
ccl21a_de_results <- list()
ccl21a_volcano_plots <- list()

for (cluster_id in levels(Idents(ccl21a_de_object))) {
  group_counts <- table(ccl21a_de_object$group[ccl21a_de_object$seurat_clusters == cluster_id])
  comparison_counts <- group_counts[c("WT_CTRL", "MT_CTRL")]
  comparison_counts[is.na(comparison_counts)] <- 0
  if (any(comparison_counts < analysis_config$min_cells_per_group)) next
  result <- FindMarkers(
    ccl21a_de_object, subset.ident = cluster_id, group.by = "group",
    ident.1 = "WT_CTRL", ident.2 = "MT_CTRL",
    logfc.threshold = analysis_config$logfc_threshold,
    min.pct = analysis_config$min_pct, test.use = "wilcox"
  )
  result$gene <- rownames(result)
  result$cluster <- cluster_id
  ccl21a_de_results[[cluster_id]] <- result
  result$category <- "Not significant"
  result$category[result$p_val_adj < analysis_config$adjusted_p_cutoff & result$avg_log2FC > analysis_config$lfc_cutoff] <- "WT higher"
  result$category[result$p_val_adj < analysis_config$adjusted_p_cutoff & result$avg_log2FC < -analysis_config$lfc_cutoff] <- "MT higher"
  ccl21a_volcano_plots[[cluster_id]] <- ggplot(result, aes(avg_log2FC, -log10(pmax(p_val_adj, .Machine$double.xmin)), color = category)) +
    geom_point(alpha = 0.6, size = 1) + theme_classic() +
    labs(title = paste("Cluster", cluster_id, "WT CTRL vs MT CTRL"), x = "Average log2 fold change", y = "-log10 adjusted p-value")
}
print(ccl21a_de_results)
if (length(ccl21a_volcano_plots)) print(wrap_plots(ccl21a_volcano_plots))

