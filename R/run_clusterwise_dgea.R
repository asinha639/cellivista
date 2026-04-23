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

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  Seurat::DefaultAssay(seurat_obj) <- "RNA"

  all_markers <- Seurat::FindAllMarkers(
    seurat_obj,
    logfc.threshold = logfc_threshold,
    min.pct = min_pct,
    only.pos = FALSE
  )
  readr::write_csv(all_markers, file.path(output_dir, "all_marker.csv"))

  top15_marker <- all_markers %>%
    dplyr::filter(avg_log2FC > 0) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice(1:15) %>%
    dplyr::select(avg_log2FC, cluster, gene)

  readr::write_csv(top15_marker %>% dplyr::select(cluster, gene),
                   file.path(output_dir, "top15_marker.csv"))

  # Skip plot generation if plot_options is NULL or no plots are selected
  if (is.null(plot_options)) return(all_markers)

  plot_dir <- file.path(output_dir, "Plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  clusters <- unique(top15_marker$cluster)

  for (clust in clusters) {
    markers_for_cluster <- top15_marker %>% dplyr::filter(cluster == clust) %>% dplyr::pull(gene)
    if (length(markers_for_cluster) == 0) next

    title_prefix <- plot_options$title_prefix %||% "Cluster"
    ext <- plot_options$format %||% "jpeg"
    width <- plot_options$width %||% 10
    height <- plot_options$height %||% 8

    # Violin plot
    if (isTRUE(plot_options$save_violin)) {
      p <- Seurat::VlnPlot(seurat_obj, features = markers_for_cluster, stack = TRUE, flip = TRUE) +
        ggplot2::ggtitle(paste(title_prefix, clust, "- Violin")) +
        ggplot2::theme(legend.position = "none")

      ggsave(file.path(plot_dir, paste0("violin_cluster", clust, ".", ext)), plot = p,
             device = ext, width = width, height = height)
    }

    # Dot plot
    if (isTRUE(plot_options$save_dot)) {
      p <- Seurat::DotPlot(seurat_obj, features = markers_for_cluster) +
        ggplot2::ggtitle(paste(title_prefix, clust, "- DotPlot")) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

      ggsave(file.path(plot_dir, paste0("dotplot_cluster", clust, ".", ext)), plot = p,
             device = ext, width = width, height = height)
    }

    # Feature plots
    if (isTRUE(plot_options$save_feature)) {
      for (gene in markers_for_cluster) {
        p <- Seurat::FeaturePlot(seurat_obj, features = gene) +
          ggplot2::ggtitle(paste(title_prefix, clust, "-", gene))

        ggsave(file.path(plot_dir, paste0("featureplot_cluster", clust, "_", gene, ".", ext)),
               plot = p, device = ext, width = width, height = height)
      }
    }
  }

  return(all_markers)
}

