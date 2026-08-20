# CCL21A cluster composition visualization
# This is descriptive because biological-replicate uncertainty is unavailable.

# Import libraries and data ------------------------------------------------------
source(file.path("config", "analysis_config.R"))
source(file.path("R", "reproducibility.R"))
require_packages(c("Seurat", "dplyr", "ggplot2", "scales"))
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2); library(scales) })
if (!exists("ccl21a_joint") || !inherits(ccl21a_joint, "Seurat")) stop("Run scripts/05_ccl21a_joint_clustering.R first.", call. = FALSE)
validate_seurat(ccl21a_joint, metadata = c("group", "seurat_clusters"))

# Actual analysis ---------------------------------------------------------------
ccl21a_composition <- ccl21a_joint[[]] |>
  count(group, seurat_clusters, name = "n") |>
  group_by(group) |>
  mutate(proportion = n / sum(n)) |>
  ungroup()
print(ccl21a_composition)
ccl21a_composition_plot <- ggplot(ccl21a_composition, aes(group, proportion, fill = seurat_clusters)) +
  geom_col(color = "black") + scale_y_continuous(labels = percent) + theme_classic() +
  labs(x = NULL, y = "Proportion", fill = "Cluster", title = "CCL21A cluster composition")
print(ccl21a_composition_plot)

