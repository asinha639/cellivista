#' Convert Matrix/Barcode/Gene Files to a Seurat RDS
#'
#' @param matrix_path Path to count matrix (.mtx or .mtx.gz).
#' @param cell_path Path to cell/barcode file.
#' @param gene_path Path to gene/features file.
#' @param output_path Path to save Seurat .rds output.
#' @param project_name Seurat project name.
#' @param sample_label_mode Mode to create SampleLabel metadata.
#' @param single_sample_label Default single/fallback sample label.
#' @param sample_label_column Column name or 1-based index for SampleLabel when
#'   `sample_label_mode = "from_cell_file_column"`.
#' @param min_cells Minimum cells per feature for CreateSeuratObject.
#' @param min_features Minimum features per cell for CreateSeuratObject.
#'
#' @return Seurat object with required SampleLabel metadata.
#' @export
convert_matrix_to_rds <- function(
  matrix_path,
  cell_path,
  gene_path,
  output_path,
  project_name = "SeuratProject",
  sample_label_mode = c("barcode_suffix", "single_sample", "from_cell_file_column"),
  single_sample_label = "Sample1",
  sample_label_column = NULL,
  min_cells = 3,
  min_features = 200
) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package is required.")
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required.")
  }
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("readr package is required.")
  }

  sample_label_mode <- match.arg(sample_label_mode)
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  for (p in c(matrix_path, cell_path, gene_path)) {
    if (is.null(p) || !nzchar(p) || !file.exists(p)) {
      stop("Input file does not exist: ", p %||% "<missing>")
    }
  }

  single_sample_label <- trimws(as.character(single_sample_label %||% "Sample1"))
  if (!nzchar(single_sample_label)) single_sample_label <- "Sample1"

  detect_delim <- function(path) {
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
    on.exit(close(con), add = TRUE)
    first_line <- readLines(con, n = 1, warn = FALSE)
    if (length(first_line) == 0) return("\t")
    line <- first_line[[1]]
    n_comma <- lengths(regmatches(line, gregexpr(",", line, fixed = TRUE)))
    n_tab <- lengths(regmatches(line, gregexpr("\t", line, fixed = TRUE)))
    if (n_comma > n_tab) "," else "\t"
  }

  read_table_no_header <- function(path) {
    delim <- detect_delim(path)
    readr::read_delim(
      file = path,
      delim = delim,
      col_names = FALSE,
      col_types = readr::cols(.default = readr::col_character()),
      progress = FALSE,
      show_col_types = FALSE,
      trim_ws = TRUE
    )
  }

  read_table_with_header <- function(path) {
    delim <- detect_delim(path)
    readr::read_delim(
      file = path,
      delim = delim,
      col_names = TRUE,
      col_types = readr::cols(.default = readr::col_character()),
      progress = FALSE,
      show_col_types = FALSE,
      trim_ws = TRUE
    )
  }

  counts <- if (grepl("\\.gz$", matrix_path, ignore.case = TRUE)) {
    con <- gzfile(matrix_path, open = "rt")
    on.exit(close(con), add = TRUE)
    Matrix::readMM(con)
  } else {
    Matrix::readMM(matrix_path)
  }
  if (!inherits(counts, "dgCMatrix")) {
    counts <- methods::as(counts, "dgCMatrix")
  }

  cell_df <- read_table_no_header(cell_path)
  gene_df <- read_table_no_header(gene_path)

  if (ncol(cell_df) < 1 || nrow(cell_df) < 1) {
    stop("Cell/barcode file must contain at least one column and one row.")
  }
  if (ncol(gene_df) < 1 || nrow(gene_df) < 1) {
    stop("Gene/features file must contain at least one column and one row.")
  }

  cell_barcodes <- as.character(cell_df[[1]])
  cell_barcodes <- trimws(cell_barcodes)
  if (any(!nzchar(cell_barcodes))) {
    stop("Cell/barcode file contains blank barcode values in the first column.")
  }

  gene_name_primary <- as.character(gene_df[[1]])
  gene_name_secondary <- if (ncol(gene_df) >= 2) as.character(gene_df[[2]]) else rep(NA_character_, length(gene_name_primary))
  gene_names <- ifelse(!is.na(gene_name_secondary) & nzchar(trimws(gene_name_secondary)), gene_name_secondary, gene_name_primary)
  gene_names <- trimws(gene_names)
  if (any(!nzchar(gene_names))) {
    stop("Gene/features file contains blank gene names after processing.")
  }
  gene_names <- make.unique(gene_names)

  if (nrow(counts) != length(gene_names)) {
    stop("Matrix row count (", nrow(counts), ") does not match gene count (", length(gene_names), ").")
  }
  if (ncol(counts) != length(cell_barcodes)) {
    stop("Matrix column count (", ncol(counts), ") does not match cell count (", length(cell_barcodes), ").")
  }

  rownames(counts) <- gene_names
  colnames(counts) <- cell_barcodes

  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts,
    project = project_name,
    min.cells = as.numeric(min_cells),
    min.features = as.numeric(min_features)
  )

  build_sample_labels <- function() {
    active_barcodes <- colnames(seurat_obj)

    if (identical(sample_label_mode, "single_sample")) {
      return(rep(single_sample_label, length(active_barcodes)))
    }

    if (identical(sample_label_mode, "barcode_suffix")) {
      suffix <- sub(".*-([0-9]+)$", "\\1", active_barcodes)
      matched <- grepl(".*-[0-9]+$", active_barcodes)
      if (!any(matched)) {
        return(rep(single_sample_label, length(active_barcodes)))
      }
      out <- rep(single_sample_label, length(active_barcodes))
      out[matched] <- paste0("Sample", suffix[matched])
      return(out)
    }

    if (is.null(sample_label_column) || !nzchar(trimws(as.character(sample_label_column)))) {
      stop("sample_label_column is required when sample_label_mode = 'from_cell_file_column'.")
    }

    labels_by_barcode <- NULL
    sample_label_column_chr <- trimws(as.character(sample_label_column))
    sample_col_index <- suppressWarnings(as.integer(sample_label_column_chr))
    if (!is.na(sample_col_index) && sample_col_index >= 1) {
      if (sample_col_index > ncol(cell_df)) {
        stop("sample_label_column index is out of bounds for the cell/barcode file.")
      }
      labels_by_barcode <- as.character(cell_df[[sample_col_index]])
    } else {
      cell_df_with_header <- read_table_with_header(cell_path)
      if (!(sample_label_column_chr %in% names(cell_df_with_header))) {
        stop("sample_label_column '", sample_label_column_chr, "' was not found in cell/barcode file headers.")
      }
      header_barcodes <- as.character(cell_df_with_header[[1]])
      header_labels <- as.character(cell_df_with_header[[sample_label_column_chr]])
      header_map <- stats::setNames(header_labels, header_barcodes)
      mapped_labels <- unname(header_map[active_barcodes])
      labels <- as.character(mapped_labels)
      labels[is.na(labels) | !nzchar(trimws(labels))] <- single_sample_label
      return(labels)
    }

    barcode_map <- stats::setNames(labels_by_barcode, cell_barcodes)
    labels <- as.character(unname(barcode_map[active_barcodes]))
    if (any(is.na(labels) | !nzchar(trimws(labels)))) {
      warning("Some SampleLabel values were missing/blank and were replaced with fallback label.")
      labels[is.na(labels) | !nzchar(trimws(labels))] <- single_sample_label
    }
    labels
  }

  sample_labels <- build_sample_labels()
  sample_labels <- as.character(sample_labels)
  sample_labels[is.na(sample_labels) | !nzchar(trimws(sample_labels))] <- single_sample_label
  seurat_obj$SampleLabel <- sample_labels

  meta <- seurat_obj[[]]
  if (!"SampleLabel" %in% names(meta)) {
    stop("Failed to create SampleLabel metadata column.")
  }
  sample_meta <- as.character(meta$SampleLabel)
  if (any(is.na(sample_meta) | !nzchar(trimws(sample_meta)))) {
    stop("SampleLabel contains NA or blank values after processing.")
  }

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(seurat_obj, output_path)
  seurat_obj
}
