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
