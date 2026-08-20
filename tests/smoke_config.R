source(file.path("config", "analysis_config.R"))
source(file.path("R", "reproducibility.R"))

stopifnot(
  analysis_config$seed == 1234L,
  grepl("/outputs$", analysis_config$output_dir),
  is.function(require_inputs),
  is.function(require_layers),
  is.function(require_barcode_overlap),
  is.function(read_barcode_only_csv),
  is.function(read_projection_csv),
  is.function(match_projection_coordinates),
  is.function(write_provenance)
)

description <- read.dcf("DESCRIPTION")
stopifnot(nrow(description) == 1L)

required_path_keys <- c(
  "wt_ctrl_h5", "mt_ctrl_h5", "wt_oxd4_h5", "mt_oxd4_h5",
  "ctrl_projection", "oxd4_projection",
  "wt_ctrl_barcodes", "mt_ctrl_barcodes", "wt_oxd4_barcodes", "mt_oxd4_barcodes",
  "wt_ctrl_ccl21a_barcodes", "mt_ctrl_ccl21a_barcodes",
  "wt_oxd4_ccl21a_barcodes", "mt_oxd4_ccl21a_barcodes",
  "ctrl_hires_png", "oxd4_hires_png", "ctrl_scale_json", "oxd4_scale_json"
)
stopifnot(all(required_path_keys %in% names(analysis_config$paths)))
message("Configuration, helpers, and DESCRIPTION smoke checks passed.")
