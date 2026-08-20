# CXCL13 cluster composition visualization
# This is descriptive because biological-replicate uncertainty is unavailable.

# Import libraries and data ------------------------------------------------------
source(file.path("config", "analysis_config.R"))
source(file.path("R", "reproducibility.R"))
require_packages(c("Seurat", "dplyr", "ggplot2", "scales"))
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2); library(scales) })
if (!exists("cxcl13_joint") || !inherits(cxcl13_joint, "Seurat")) stop("Run scripts/01_cxcl13_joint_clustering.R first.", call. = FALSE)
validate_seurat(cxcl13_joint, metadata = c("group", "seurat_clusters"))

# Actual analysis ---------------------------------------------------------------
cxcl13_composition <- cxcl13_joint[[]] |>
  count(group, seurat_clusters, name = "n") |>
  group_by(group) |>
  mutate(proportion = n / sum(n)) |>
  ungroup()
print(cxcl13_composition)
cxcl13_composition_plot <- ggplot(cxcl13_composition, aes(group, proportion, fill = seurat_clusters)) +
  geom_col(color = "black") +
  scale_y_continuous(labels = percent) +
  theme_classic() +
  labs(x = NULL, y = "Proportion", fill = "Cluster", title = "CXCL13 cluster composition")
print(cxcl13_composition_plot)

