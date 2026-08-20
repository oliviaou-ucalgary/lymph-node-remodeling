require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required package(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

require_inputs <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing required input(s):\n- ", paste(missing, collapse = "\n- "), call. = FALSE)
  invisible(normalizePath(paths, winslash = "/", mustWork = TRUE))
}

ensure_output_dir <- function(path) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create output directory: ", path, call. = FALSE)
  }
  invisible(path)
}

select_gene_expression <- function(x, source) {
  if (!is.list(x)) return(x)
  if ("Gene Expression" %in% names(x)) return(x[["Gene Expression"]])
  stop("Multimodal input lacks a 'Gene Expression' matrix: ", source, call. = FALSE)
}

validate_seurat <- function(object, assay = "RNA", metadata = character()) {
  if (!inherits(object, "Seurat")) stop("Expected a Seurat object.", call. = FALSE)
  if (!assay %in% SeuratObject::Assays(object)) stop("Missing assay: ", assay, call. = FALSE)
  missing_meta <- setdiff(metadata, colnames(object[[]]))
  if (length(missing_meta)) stop("Missing metadata: ", paste(missing_meta, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

require_layers <- function(object, assay, layers) {
  validate_seurat(object, assay = assay)
  available <- SeuratObject::Layers(object[[assay]])
  missing <- setdiff(layers, available)
  if (length(missing)) stop("Missing layer(s) in assay ", assay, ": ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

require_reductions <- function(object, reductions) {
  missing <- setdiff(reductions, names(object@reductions))
  if (length(missing)) stop("Missing reduction(s): ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

require_barcode_overlap <- function(object, barcodes, label, minimum = 1L) {
  matched <- intersect(colnames(object), unique(as.character(barcodes)))
  if (length(matched) < minimum) {
    stop(label, " matched ", length(matched), " barcode(s); expected at least ", minimum,
         ". Check barcode suffixes and input provenance.", call. = FALSE)
  }
  matched
}

read_barcode_only_csv <- function(path) {
  require_inputs(path)
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(x) != 1L) stop("Expected exactly one barcode column in: ", path, call. = FALSE)
  barcodes <- trimws(as.character(x[[1]]))
  barcodes <- barcodes[!is.na(barcodes) & nzchar(barcodes)]
  if (!length(barcodes)) stop("No usable barcodes in: ", path, call. = FALSE)
  if (anyDuplicated(barcodes)) message("Removing duplicated barcodes from: ", path)
  unique(barcodes)
}

read_projection_csv <- function(path) {
  require_inputs(path)
  projection <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("Barcode", "X Coordinate", "Y Coordinate")
  missing <- setdiff(required, colnames(projection))
  if (length(missing)) stop("Projection file is missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  projection <- projection[, required, drop = FALSE]
  projection$Barcode <- trimws(as.character(projection$Barcode))
  projection[["X Coordinate"]] <- suppressWarnings(as.numeric(projection[["X Coordinate"]]))
  projection[["Y Coordinate"]] <- suppressWarnings(as.numeric(projection[["Y Coordinate"]]))
  if (anyNA(projection$Barcode) || any(!nzchar(projection$Barcode))) stop("Projection contains missing barcodes: ", path, call. = FALSE)
  if (anyNA(projection[["X Coordinate"]]) || anyNA(projection[["Y Coordinate"]])) stop("Projection contains missing or nonnumeric coordinates: ", path, call. = FALSE)
  unique_positions <- unique(projection)
  ambiguous <- unique_positions$Barcode[duplicated(unique_positions$Barcode) | duplicated(unique_positions$Barcode, fromLast = TRUE)]
  if (length(ambiguous)) stop("Projection maps barcode(s) to multiple positions, including: ", paste(head(unique(ambiguous), 10), collapse = ", "), call. = FALSE)
  unique_positions[!duplicated(unique_positions$Barcode), , drop = FALSE]
}

match_projection_coordinates <- function(barcodes, projection, label) {
  matched <- merge(data.frame(Barcode = unique(barcodes), stringsAsFactors = FALSE), projection,
                   by = "Barcode", all.x = TRUE, sort = FALSE)
  n_matched <- sum(stats::complete.cases(matched[, c("X Coordinate", "Y Coordinate")]))
  message(label, ": matched ", n_matched, " of ", nrow(matched), " barcode(s).")
  if (!n_matched) stop("No projection coordinates matched for ", label, call. = FALSE)
  matched
}

write_csv_with_id <- function(x, path, id_name = "feature") {
  ensure_output_dir(dirname(path))
  out <- as.data.frame(x, check.names = FALSE)
  if (!is.null(rownames(out))) out <- cbind(stats::setNames(data.frame(rownames(out)), id_name), out)
  utils::write.csv(out, path, row.names = FALSE, na = "")
  invisible(path)
}

save_png <- function(plot, path, width, height, dpi = analysis_config$figure_dpi) {
  ensure_output_dir(dirname(path))
  ggplot2::ggsave(path, plot = plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  invisible(path)
}

write_provenance <- function(output_dir, parameters = analysis_config, object = NULL) {
  ensure_output_dir(output_dir)
  capture.output(str(parameters), file = file.path(output_dir, "analysis_parameters.txt"))
  capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))
  loaded_namespaces <- loadedNamespaces()
  namespace_versions <- vapply(loaded_namespaces, function(x) as.character(utils::packageVersion(x)), character(1))
  versions <- data.frame(
    item = c("timestamp_utc", "seed", "R", "Seurat", "SeuratObject"),
    value = c(format(Sys.time(), tz = "UTC", usetz = TRUE), parameters$seed,
              R.version.string,
              if (requireNamespace("Seurat", quietly = TRUE)) as.character(utils::packageVersion("Seurat")) else NA,
              if (requireNamespace("SeuratObject", quietly = TRUE)) as.character(utils::packageVersion("SeuratObject")) else NA)
  )
  utils::write.csv(versions, file.path(output_dir, "software_versions.csv"), row.names = FALSE)
  utils::write.csv(data.frame(package = loaded_namespaces, version = unname(namespace_versions)),
                   file.path(output_dir, "loaded_package_versions.csv"), row.names = FALSE)
  capture.output(Sys.info(), file = file.path(output_dir, "system_info.txt"))
  if (!is.null(object)) {
    audit <- data.frame(
      class = paste(class(object), collapse = ";"), features = nrow(object), cells = ncol(object),
      default_assay = SeuratObject::DefaultAssay(object),
      assays = paste(SeuratObject::Assays(object), collapse = ";"),
      reductions = paste(names(object@reductions), collapse = ";"),
      images = paste(names(object@images), collapse = ";"),
      identities = paste(levels(SeuratObject::Idents(object)), collapse = ";"),
      layers = paste(unlist(lapply(SeuratObject::Assays(object), function(a) {
        paste0(a, ":", paste(SeuratObject::Layers(object[[a]]), collapse = "|"))
      })), collapse = ";")
    )
    utils::write.csv(audit, file.path(output_dir, "seurat_object_audit.csv"), row.names = FALSE)
    utils::write.csv(data.frame(metadata_column = colnames(object[[]])), file.path(output_dir, "metadata_columns.csv"), row.names = FALSE)
  }
  invisible(output_dir)
}
