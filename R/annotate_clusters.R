#' Annotate Clusters with User-Defined Cell Types and Visualize
#'
#' @param seurat_obj A clustered Seurat object (with `seurat_clusters`).
#' @param annotations A character vector of annotation labels, or a CSV path.
#' @param output_dir Directory to save plots. Default is "./Plots/Annotation".
#' @param rds_path Path to save annotated object. Default is "RDS_Files/annotated_seurat.obj.rds".
#' @param save_rds Logical; whether to save the annotated Seurat object. Default is TRUE.
#' @param save_plots Logical; whether to save annotation plots. Default is TRUE.
#' @param plot_options Optional list of plot control parameters for each annotation plot.
#'
#' @return Annotated Seurat object.
#' @export
annotate_clusters <- function(seurat_obj,
                              annotations,
                              output_dir = "./Plots/Annotation",
                              rds_path = "RDS_Files/annotated_seurat.obj.rds",
                              save_rds = TRUE,
                              save_plots = TRUE,
                              plot_options = NULL) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required but not installed.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required but not installed.")
  }
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("readr is required but not installed.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr is required but not installed.")
  }
  
  if (is.null(seurat_obj)) {
    stop("seurat_obj is NULL.")
  }
  if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    stop("seurat_obj must contain 'seurat_clusters' in metadata.")
  }
  if (!"SampleLabel" %in% colnames(seurat_obj@meta.data)) {
    stop("seurat_obj must contain 'SampleLabel' in metadata.")
  }
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  save_plot <- function(plot_obj, opts = NULL, default_filename, default_title, default_format = "jpeg", default_width = 10, default_height = 8) {
    opts <- opts %||% list()
    format <- tolower(opts$format %||% default_format)
    device <- switch(format,
                     jpg = "jpeg",
                     jpeg = "jpeg",
                     png = "png",
                     svg = "svg",
                     format)
    filename <- paste0(
      tools::file_path_sans_ext(opts$filename %||% default_filename),
      ".",
      format
    )
    plot_title <- opts$title %||% default_title
    plot_obj <- plot_obj + ggplot2::ggtitle(if (nzchar(plot_title)) plot_title else NULL)
    ggplot2::ggsave(
      filename = file.path(output_dir, filename),
      plot = plot_obj,
      device = device,
      width = as.numeric(opts$width %||% default_width),
      height = as.numeric(opts$height %||% default_height)
    )
    invisible(file.path(output_dir, filename))
  }
  
  cluster_levels <- sort(unique(as.character(seurat_obj@meta.data$seurat_clusters)))
  
  # Build annotation table
  if (is.character(annotations) && length(annotations) == 1 && file.exists(annotations)) {
    cluster_annot <- readr::read_csv(annotations, show_col_types = FALSE)
    
    required_cols <- c("seurat_clusters", "cell_annotation")
    missing_cols <- setdiff(required_cols, colnames(cluster_annot))
    if (length(missing_cols) > 0) {
      stop("Annotation CSV must contain columns: seurat_clusters, cell_annotation")
    }
    
    cluster_annot <- cluster_annot[, required_cols, drop = FALSE]
    cluster_annot$seurat_clusters <- as.character(cluster_annot$seurat_clusters)
    cluster_annot$cell_annotation <- as.character(cluster_annot$cell_annotation)
    
  } else if (is.character(annotations) && length(annotations) == length(cluster_levels)) {
    cluster_annot <- data.frame(
      seurat_clusters = cluster_levels,
      cell_annotation = as.character(annotations),
      stringsAsFactors = FALSE
    )
  } else {
    stop("annotations must be either a valid CSV path or a character vector matching the number of clusters.")
  }
  
  # Validate coverage
  missing_clusters <- setdiff(cluster_levels, cluster_annot$seurat_clusters)
  extra_clusters <- setdiff(cluster_annot$seurat_clusters, cluster_levels)
  
  if (length(missing_clusters) > 0) {
    stop("Missing annotations for cluster(s): ", paste(missing_clusters, collapse = ", "))
  }
  if (length(extra_clusters) > 0) {
    warning("Annotation table contains unused cluster(s): ", paste(extra_clusters, collapse = ", "))
  }
  
  cluster_annot <- cluster_annot[match(cluster_levels, cluster_annot$seurat_clusters), , drop = FALSE]
  
  # Map annotations safely without row reordering
  annotation_map <- stats::setNames(cluster_annot$cell_annotation, cluster_annot$seurat_clusters)
  seurat_obj$cell_annotation <- unname(annotation_map[as.character(seurat_obj$seurat_clusters)])
  
  if (any(is.na(seurat_obj$cell_annotation))) {
    stop("Failed to assign annotations to all cells.")
  }
  
  seurat_obj$cell_annotation <- factor(
    seurat_obj$cell_annotation,
    levels = unique(cluster_annot$cell_annotation)
  )
  Seurat::Idents(seurat_obj) <- seurat_obj$cell_annotation
  
  umap_annotated_opts <- plot_options$umap_annotated %||% list()
  umap_split_opts <- plot_options$umap_split_by_sample %||% list()
  pie_opts <- plot_options$celltype_piechart %||% list()
  save_umap_annotated <- isTRUE(save_plots) && (
    is.null(plot_options) || is.null(plot_options$umap_annotated) || isTRUE(umap_annotated_opts$save)
  )
  save_umap_split <- isTRUE(save_plots) && (
    is.null(plot_options) || is.null(plot_options$umap_split_by_sample) || isTRUE(umap_split_opts$save)
  )
  save_pie <- isTRUE(save_plots) && (
    is.null(plot_options) || is.null(plot_options$celltype_piechart) || isTRUE(pie_opts$save)
  )
  
  umap_annotated <- Seurat::DimPlot(
    seurat_obj,
    reduction = "umap",
    label = FALSE,
    group.by = "cell_annotation"
  )
  if (save_umap_annotated) {
    save_plot(
      plot_obj = umap_annotated,
      opts = umap_annotated_opts,
      default_filename = "umap_annotated",
      default_title = "UMAP: Annotated Cell Types",
      default_width = 10,
      default_height = 8
    )
  }
  
  umap_split <- Seurat::DimPlot(
    seurat_obj,
    reduction = "umap",
    split.by = "SampleLabel",
    group.by = "cell_annotation",
    label = FALSE
  )
  
  if (save_umap_split) {
    save_plot(
      plot_obj = umap_split,
      opts = umap_split_opts,
      default_filename = "umap_split_by_sample",
      default_title = "UMAP: Annotated Cell Types by Sample",
      default_width = 18,
      default_height = 8
    )
  }
  
  pie_data <- as.data.frame(table(seurat_obj$cell_annotation), stringsAsFactors = FALSE)
  colnames(pie_data) <- c("CellType", "Count")
  
  pie_plot <- ggplot2::ggplot(
    pie_data,
    ggplot2::aes(x = "", y = Count, fill = CellType)
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::theme_minimal() +
    ggplot2::geom_text(
      ggplot2::aes(label = Count),
      position = ggplot2::position_stack(vjust = 0.5)
    ) +
    ggplot2::labs(x = "", y = "") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
  
  if (save_pie) {
    save_plot(
      plot_obj = pie_plot,
      opts = pie_opts,
      default_filename = "celltype_piechart",
      default_title = "Cell Type Composition",
      default_width = 10,
      default_height = 10
    )
  }
  
  if (isTRUE(save_rds)) {
    dir.create(dirname(rds_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(seurat_obj, rds_path)
  }
  
  return(seurat_obj)
}

