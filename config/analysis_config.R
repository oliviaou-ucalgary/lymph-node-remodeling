# Non-secret, version-controlled defaults. Copy settings that differ locally to
# config/local.R; that file is ignored by Git.
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

analysis_config <- list(
  seed = 1234L,
  input_dir = file.path(project_root, "data"),
  output_dir = file.path(project_root, "outputs"),
  figure_dpi = 300L,
  assay = "RNA",
  normalized_layer = "data",
  variable_features = 2000L,
  normalization_scale_factor = 10000,
  pca_dims = 1:15,
  pca_dims_30 = 1:30,
  clustering_resolution = 0.6,
  clustering_resolution_30 = 0.4,
  umap_neighbors = 15L,
  umap_min_dist = 0.05,
  min_cells_per_group = 20L,
  min_pct = 0.1,
  logfc_threshold = 0,
  adjusted_p_cutoff = 0.05,
  lfc_cutoff = 0.25,
  top_genes = 10L,
  heatmap_dpi = 300L,
  # Loupe Browser selections are external inputs. R validates and projects the
  # selected barcodes but does not recreate the manual selection.
  paths = list(
    wt_ctrl_h5 = file.path(project_root, "data", "ctrl", "wt ctrl ln1.h5"),
    mt_ctrl_h5 = file.path(project_root, "data", "ctrl", "mt ctrl ln1.h5"),
    wt_oxd4_h5 = file.path(project_root, "data", "oxd4", "wt oxd4 ln2.h5"),
    mt_oxd4_h5 = file.path(project_root, "data", "oxd4", "mt oxd4 ln3.h5"),
    wt_ctrl_barcodes = file.path(project_root, "data", "ctrl", "cxcl13", "wt ctrl ln1 cxcl13 coordinate.csv"),
    mt_ctrl_barcodes = file.path(project_root, "data", "ctrl", "cxcl13", "mt ctrl ln1 cxcl13 coordinate.csv"),
    wt_oxd4_barcodes = file.path(project_root, "data", "oxd4", "cxcl13", "wt oxd4 ln2 cxcl13 coordinate.csv"),
    mt_oxd4_barcodes = file.path(project_root, "data", "oxd4", "cxcl13", "mt oxd4 ln3 cxcl13 coordinate.csv"),
    wt_ctrl_ccl21a_barcodes = file.path(project_root, "data", "ctrl", "ccl21", "wt ctrl ln1 ccl21 coordinate.csv"),
    mt_ctrl_ccl21a_barcodes = file.path(project_root, "data", "ctrl", "ccl21", "mt ctrl ln1 ccl21 coordinate.csv"),
    wt_oxd4_ccl21a_barcodes = file.path(project_root, "data", "oxd4", "ccl21", "wt oxd4 ln2 ccl21 coordinate.csv"),
    mt_oxd4_ccl21a_barcodes = file.path(project_root, "data", "oxd4", "ccl21", "mt oxd4 ln3 ccl21 coordinate.csv"),
    ctrl_projection = file.path(project_root, "data", "ctrl", "projection.csv"),
    oxd4_projection = file.path(project_root, "data", "oxd4", "projection.csv"),
    joint_rds = file.path(project_root, "outputs", "cxcl13_joint", "cxcl13_joint4_umap.rds"),
    scrna_rds = file.path(project_root, "data", "derived", "joint_WT_Lysm_CD45_filtered_umap.rds"),
    ctrl_spatial_dir = file.path(project_root, "data", "ctrl", "binned_outputs", "square_008um"),
    oxd4_spatial_dir = file.path(project_root, "data", "oxd4", "binned_outputs", "square_008um"),
    ctrl_barcode_map = file.path(project_root, "data", "ctrl", "barcode_mappings.parquet"),
    oxd4_barcode_map = file.path(project_root, "data", "oxd4", "barcode_mappings.parquet"),
    ctrl_hires_png = file.path(project_root, "data", "ctrl", "tissue_hires_image.png"),
    oxd4_hires_png = file.path(project_root, "data", "oxd4", "tissue_hires_image.png"),
    ctrl_scale_json = file.path(project_root, "data", "ctrl", "scalefactors_json.json"),
    oxd4_scale_json = file.path(project_root, "data", "oxd4", "scalefactors_json.json")
  )
)

local_config <- file.path(project_root, "config", "local.R")
if (file.exists(local_config)) source(local_config, local = TRUE)
set.seed(analysis_config$seed)
