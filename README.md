# Paper Seurat analysis

This repository contains interactive R/Seurat v5 workflows for `Cxcl13`- and `Ccl21a`-selected cells/spots from a mouse lymph-node model, plus a G2M spatial analysis used for a reviewer response. Raw and controlled data are intentionally not included.

## Public-release status

The code has been made portable, but the workflows cannot be run from this checkout alone: required 10x matrices, spatial images/coordinates, barcode selections, and derived Seurat objects must be obtained through the study's approved data-access route. Before release, replace the placeholder author/contact details, add the manuscript citation and selected license, and confirm every input is de-identified and legally shareable.

## Setup

Use R 4.3 or later and Seurat 5. Install dependencies listed in `DESCRIPTION`. For an exact environment, initialize `renv` on the analysis machine and commit the generated `renv.lock`; no lockfile is supplied because the original R environment was unavailable during repository preparation.

All scripts expect the repository root as the working directory. Version-controlled defaults are in `config/analysis_config.R`. Machine-specific overrides belong in `config/local.R`, which is ignored. Never put credentials or patient identifiers in configuration.

Expected private inputs are organized under `data/`:

- `ctrl/` and `oxd4/`: LN-specific H5 matrices, whole-slide `projection.csv`, tissue image, scale factors, and barcode mappings.
- `ctrl/cxcl13`, `ctrl/ccl21`, `oxd4/cxcl13`, and `oxd4/ccl21`: Loupe Browser barcode-only CSVs. Loupe selection is manual and is not reproduced in R.
- `derived/`: upstream Seurat RDS inputs.

## Interactive workflows and execution order

Open the repository as an RStudio project or set the working directory to the repository root. Source scripts in the same RStudio session because downstream scripts consume in-memory objects and do not read or write checkpoints.

CXCL13 workflow:

1. `scripts/01_cxcl13_joint_clustering.R`
2. `scripts/02_cxcl13_differential_expression.R`
3. `scripts/03_cxcl13_cluster_composition.R`
4. `scripts/04_cxcl13_gene_plots.R`

CCL21A workflow (same steps with distinct Loupe inputs):

1. `scripts/05_ccl21a_joint_clustering.R`
2. `scripts/06_ccl21a_differential_expression.R`
3. `scripts/07_ccl21a_cluster_composition.R`
4. `scripts/08_ccl21a_gene_plots.R`

Reviewer-response workflow:

1. `scripts/09_g2m_reviewer_response.R`

Interactive tissue overlays:

1. `scripts/10_spatial_subset_overlay.R`

The overlay script joins each barcode-only Loupe export directly to `Barcode` in the corresponding slide's `projection.csv`, using the required `X Coordinate` and `Y Coordinate` columns. Alignment defaults to no transformation and is controlled by a clearly labeled optional adjustment block.

The new `scripts/` files are the active publication workflow. Files under the older singular `script/` directory are retained only as a recoverable legacy reference and should not be used as the execution guide.

The active scripts display objects, tables, and plots in RStudio. They do not call `write.csv()`, `ggsave()`, `png()`, `saveRDS()`, or provenance-writing functions. Optional gene plotting vectors are empty by default; comments record genes used previously.



## Interpretation limitations

- Current cell/spot-level Wilcoxon tests are exploratory. Cells or spots are not biological replicates. Condition-level inference requires donor/animal identifiers and replicate-aware pseudobulk or mixed-model analysis.
- The joint workflow merges and normalizes samples; it is not batch integration. Confirm common chemistry and assess batch effects before interpreting clusters.
- CXCL13- and CCL21A-selected subsets originate from external Loupe Browser barcode selections. Preserve the exported CSVs and selection notes as input provenance.
- Spatial bins can contain mixtures of cells; cell-state and cell-cycle labels require appropriate caution.
- Composition plots are descriptive unless uncertainty is estimated across biological replicates.
- Manual image transforms in the overlay workflow must be visually verified against tissue landmarks.


## Reproducibility checklist

From the repository root, parse all scripts and source the relevant sequence interactively. Confirm expected cell/spot counts, assays, layers, images, metadata, barcode overlap, and plot appearance. Record platform/chemistry, sample identities, biological replicates, batch structure, QC rules, and exact Loupe selection provenance in the manuscript or data-access documentation.
