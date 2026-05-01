#' Run Cluster-wise DGEA for All Clusters
#'
#' Performs DGEA and optionally saves violin, dot, and feature plots for top markers.
#'
#' @param seurat_obj A Seurat object.
#' @param output_dir Directory to save output CSVs and plots.
#' @param logfc_threshold Log fold change threshold.
#' @param min_pct Minimum percentage expression.
#' @param plot_options List of plot control parameters:
#'        - save_violin, save_dot, save_feature (logical)
#'        - format (jpeg, png, svg), width, height
#'        - title_prefix (optional)
#'
#' @return A data frame of marker genes.
#' @export
run_clusterwise_dgea <- function(seurat_obj,
                                 output_dir = "./DGEA",
                                 logfc_threshold = 1,
                                 min_pct = 0.25,
                                 plot_options = NULL) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required but not installed.")
  }
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("readr is required but not installed.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr is required but not installed.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required but not installed.")
  }
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  sanitize_file_component <- function(x) gsub("[^A-Za-z0-9_\\-]", "_", x)
  
  if (is.null(seurat_obj)) {
    stop("seurat_obj is NULL.")
  }
  if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    stop("seurat_clusters not found. Please run downstream clustering before cluster-wise DGEA.")
  }
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (!"SCT" %in% Seurat::Assays(seurat_obj)) {
    stop("SCT assay not found. Please run SCTransform/integration before cluster-wise DGEA.")
  }
  Seurat::DefaultAssay(seurat_obj) <- "SCT"
  seurat_obj <- Seurat::PrepSCTFindMarkers(object = seurat_obj)
  Seurat::Idents(seurat_obj) <- seurat_obj$seurat_clusters
  
  all_markers <- Seurat::FindAllMarkers(
    object = seurat_obj,
    logfc.threshold = logfc_threshold,
    min.pct = min_pct,
    only.pos = FALSE
  )
  
  if (!is.data.frame(all_markers) || nrow(all_markers) == 0) {
    stop("No marker genes detected from SCT assay.")
  }
  
  readr::write_csv(all_markers, file.path(output_dir, "all_marker.csv"))
  
  top15_marker <- all_markers %>%
    dplyr::filter(.data$avg_log2FC > 0) %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::arrange(dplyr::desc(.data$avg_log2FC), .by_group = TRUE) %>%
    dplyr::slice_head(n = 15) %>%
    dplyr::ungroup() %>%
    dplyr::select(.data$avg_log2FC, .data$cluster, .data$gene)
  
  readr::write_csv(
    top15_marker %>% dplyr::select(.data$cluster, .data$gene),
    file.path(output_dir, "top15_marker.csv")
  )
  
  if (is.null(plot_options)) {
    return(all_markers)
  }
  
  plot_dir <- file.path(output_dir, "Plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  save_plot <- function(plot_obj, kind_opts, default_filename, default_title, filename_suffix = NULL, default_width = 10, default_height = 8) {
    kind_opts <- kind_opts %||% list()
    plot_options_fallback <- plot_options %||% list()
    title_prefix <- kind_opts$title_prefix %||% plot_options_fallback$title_prefix %||% "Cluster"
    format <- tolower(kind_opts$format %||% plot_options_fallback$format %||% "jpeg")
    device <- switch(format,
                     jpg = "jpeg",
                     jpeg = "jpeg",
                     png = "png",
                     svg = "svg",
                     format)
    width <- as.numeric(kind_opts$width %||% plot_options_fallback$width %||% default_width)
    height <- as.numeric(kind_opts$height %||% plot_options_fallback$height %||% default_height)
    filename_prefix <- kind_opts$filename_prefix %||% plot_options_fallback$filename_prefix %||% default_filename
    filename_prefix <- tools::file_path_sans_ext(filename_prefix)
    plot_title <- if (nzchar(title_prefix)) paste(title_prefix, default_title) else default_title
    plot_obj <- plot_obj + ggplot2::ggtitle(plot_title)
    file_name <- paste0(filename_prefix, filename_suffix %||% "", ".", format)
    ggplot2::ggsave(
      filename = file.path(plot_dir, file_name),
      plot = plot_obj,
      device = device,
      width = width,
      height = height
    )
    invisible(file.path(plot_dir, file_name))
  }
  
  clusters <- unique(top15_marker$cluster)
  
  for (clust in clusters) {
    markers_for_cluster <- top15_marker %>%
      dplyr::filter(.data$cluster == clust) %>%
      dplyr::pull(.data$gene)
    
    markers_for_cluster <- unique(markers_for_cluster)
    markers_for_cluster <- markers_for_cluster[markers_for_cluster %in% rownames(seurat_obj)]
    
    if (length(markers_for_cluster) == 0) next
    
    violin_opts <- plot_options$violin %||% list()
    dot_opts <- plot_options$dot %||% list()
    feature_opts <- plot_options$feature %||% list()
    
    violin_save <- isTRUE(violin_opts$save) || isTRUE(plot_options$save_violin)
    dot_save <- isTRUE(dot_opts$save) || isTRUE(plot_options$save_dot)
    feature_save <- isTRUE(feature_opts$save) || isTRUE(plot_options$save_feature)
    
    if (violin_save) {
      p <- Seurat::VlnPlot(
        seurat_obj,
        features = markers_for_cluster,
        stack = TRUE,
        flip = TRUE
      ) +
        ggplot2::theme(legend.position = "none")
      
      save_plot(
        plot_obj = p,
        kind_opts = violin_opts,
        default_filename = "violin_cluster",
        default_title = paste(clust, "- Violin"),
        filename_suffix = clust,
        default_width = 10,
        default_height = 8
      )
    }
    
    if (dot_save) {
      p <- Seurat::DotPlot(
        seurat_obj,
        features = markers_for_cluster
      ) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      
      save_plot(
        plot_obj = p,
        kind_opts = dot_opts,
        default_filename = "dotplot_cluster",
        default_title = paste(clust, "- DotPlot"),
        filename_suffix = clust,
        default_width = 10,
        default_height = 8
      )
    }
    
    if (feature_save) {
      for (gene in markers_for_cluster) {
        p <- Seurat::FeaturePlot(seurat_obj, features = gene) +
          ggplot2::ggtitle(NULL)
        
        save_plot(
          plot_obj = p,
          kind_opts = feature_opts,
          default_filename = "featureplot_cluster",
          default_title = paste(clust, "-", gene),
          filename_suffix = paste0(clust, "_", sanitize_file_component(gene)),
          default_width = 10,
          default_height = 8
        )
      }
    }
  }
  
  return(all_markers)
}

