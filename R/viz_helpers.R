#' Detect the Cluster Column in Seurat Metadata
#'
#' Looks for the cluster assignments written by the Downstream Analysis step:
#' `seurat_clusters` if present, otherwise the first `*_snn_res.*` column
#' written by `Seurat::FindClusters()`.
#'
#' @param seurat_obj A Seurat object.
#'
#' @return The name of the cluster column, or NULL if none is found (callers
#'   then fall back to the active identities).
#' @export
#'
#' @examples
#' \dontrun{
#' detect_cluster_column(seurat_obj)
#' }
detect_cluster_column <- function(seurat_obj) {
  if (is.null(seurat_obj) || !inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  meta_cols <- colnames(seurat_obj@meta.data)
  if ("seurat_clusters" %in% meta_cols) {
    return("seurat_clusters")
  }
  res_cols <- grep("_snn_res\\.", meta_cols, value = TRUE)
  if (length(res_cols) > 0) {
    return(res_cols[[1]])
  }
  NULL
}

#' Detect the Sample Column in Seurat Metadata
#'
#' Returns `SampleLabel`, the sample column used throughout Cellivista, when it
#' exists and contains at least one non-empty value.
#'
#' @param seurat_obj A Seurat object.
#'
#' @return `"SampleLabel"`, or NULL if the object has no usable sample column.
#' @export
#'
#' @examples
#' \dontrun{
#' detect_sample_column(seurat_obj)
#' }
detect_sample_column <- function(seurat_obj) {
  if (is.null(seurat_obj) || !inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  if (!"SampleLabel" %in% colnames(seurat_obj@meta.data)) {
    return(NULL)
  }
  labels <- as.character(seurat_obj@meta.data$SampleLabel)
  if (all(is.na(labels) | !nzchar(trimws(labels)))) {
    return(NULL)
  }
  "SampleLabel"
}

#' List Categorical Metadata Columns for Grouping Cells
#'
#' Returns metadata columns that can group cells in a plot: factor, character,
#' or logical columns with at least one value and at most `max_levels`
#' distinct values. Numeric QC columns and per-cell identifiers are left out.
#'
#' @param seurat_obj A Seurat object.
#' @param max_levels Maximum number of distinct values a column may have.
#'   Default is 200.
#'
#' @return A character vector of column names (possibly empty).
#' @export
#'
#' @examples
#' \dontrun{
#' list_group_columns(seurat_obj)
#' }
list_group_columns <- function(seurat_obj, max_levels = 200) {
  if (is.null(seurat_obj) || !inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  meta <- seurat_obj@meta.data
  keep <- vapply(colnames(meta), function(col) {
    values <- meta[[col]]
    if (!(is.factor(values) || is.character(values) || is.logical(values))) {
      return(FALSE)
    }
    n_values <- length(unique(stats::na.omit(as.character(values))))
    n_values >= 1 && n_values <= max_levels
  }, logical(1))
  colnames(meta)[keep]
}

#' Choose the Default Grouping for Visualization Plots
#'
#' Prefers the cell-type annotation written by [annotate_clusters()]
#' (`cell_annotation`), then the cluster column from [detect_cluster_column()],
#' then the active identity.
#'
#' @param seurat_obj A Seurat object.
#'
#' @return A metadata column name, or `"ident"` for the active identity.
#' @export
#'
#' @examples
#' \dontrun{
#' default_group_column(seurat_obj)
#' }
default_group_column <- function(seurat_obj) {
  if (is.null(seurat_obj) || !inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  if ("cell_annotation" %in% colnames(seurat_obj@meta.data)) {
    labels <- as.character(seurat_obj@meta.data$cell_annotation)
    if (!all(is.na(labels) | !nzchar(trimws(labels)))) {
      return("cell_annotation")
    }
  }
  cluster_col <- detect_cluster_column(seurat_obj)
  if (is.null(cluster_col)) "ident" else cluster_col
}

# Resolves a grouping choice ("ident", a metadata column, or NULL for the
# default) into per-cell group values and a display label.
resolve_viz_groups <- function(seurat_obj, group_by = NULL) {
  if (is.null(group_by)) {
    group_by <- default_group_column(seurat_obj)
  }
  if (identical(group_by, "ident")) {
    values <- Seurat::Idents(seurat_obj)
    label <- "Identity"
  } else if (group_by %in% colnames(seurat_obj@meta.data)) {
    values <- seurat_obj@meta.data[[group_by]]
    label <- if (identical(group_by, "seurat_clusters")) "Cluster" else group_by
  } else {
    stop("Grouping column not found in metadata: ", group_by)
  }
  group_levels <- if (is.factor(values)) {
    levels(droplevels(values))
  } else {
    sort(unique(stats::na.omit(as.character(values))))
  }
  values <- factor(as.character(values), levels = group_levels)
  names(values) <- colnames(seurat_obj)
  list(column = group_by, values = values, levels = group_levels, label = label)
}

#' Save a Visualization Plot to Disk
#'
#' Writes a ggplot object as png, pdf, or svg. png and svg use the same
#' `ggplot2::ggsave()` devices as the other Cellivista plot exports; pdf uses
#' base `grDevices::pdf` so it does not require cairo.
#'
#' @param plot_obj A ggplot object.
#' @param path Output path. Its extension is replaced to match `format`.
#' @param format One of "png", "pdf", or "svg". Default is "png".
#' @param width Width in inches. Default is 10.
#' @param height Height in inches. Default is 8.
#'
#' @return The path of the written file, invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#' save_viz_plot(p, "Plots/umap_palette", format = "pdf")
#' }
save_viz_plot <- function(plot_obj, path, format = c("png", "pdf", "svg"),
                          width = 10, height = 8) {
  format <- match.arg(format)
  width <- as.numeric(width)
  height <- as.numeric(height)
  if (length(width) != 1 || is.na(width) || width <= 0 ||
      length(height) != 1 || is.na(height) || height <= 0) {
    stop("width and height must be positive numbers.")
  }

  path <- paste0(tools::file_path_sans_ext(path), ".", format)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  device <- if (identical(format, "pdf")) grDevices::pdf else format

  ggplot2::ggsave(
    filename = path,
    plot = plot_obj,
    device = device,
    width = width,
    height = height
  )
  invisible(path)
}
