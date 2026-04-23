#' Performs scaling, dimensionality reduction, clustering, and generates summary plots.
#'
#' @param seurat_obj A Seurat object (e.g., integrated dataset).
#' @param dims Dimensions to use for PCA/UMAP. Default is 1:15.
#' @param resolution Clustering resolution. Default is 0.1.
#' @param output_dir Directory to save plots. Default is "./Plots".
#' @param save_rds Logical; whether to save the Seurat object. Default is TRUE.
#' @param rds_path Path to save the processed object. Default is "RDS_Files/integrated_postclustering.rds".
#'
#' @return A Seurat object with downstream analysis applied.
#' @export
#'
#' @examples
#' integrated <- run_downstream_analysis(integrated)


  #' Run PCA, UMAP, Clustering, and Downstream Visualization
  #'
  #' Performs scaling, dimensionality reduction, clustering, and generates summary plots.
  #'
  #' @param seurat_obj A Seurat object (e.g., integrated dataset).
  #' @param dims Dimensions to use for PCA/UMAP. Default is 1:15.
  #' @param resolution Clustering resolution. Default is 0.1.
  #' @param output_dir Directory to save plots. Default is "./Plots".
  #' @param save_rds Logical; whether to save the Seurat object. Default is TRUE.
  #' @param rds_path Path to save the processed object. Default is "RDS_Files/downstream_analyzed.rds".
  #' @param plot_options Optional list of plot settings for saving (filenames, formats, sizes).
  #'
  #' @return A Seurat object with downstream analysis applied.
  #' @export
  #'
  #' @examples
  #' integrated <- run_downstream_analysis(integrated)
  run_downstream_analysis <- function(seurat_obj,
                                      dims = 1:15,
                                      resolution = 0.1,
                                      output_dir = "./Plots",
                                      save_rds = TRUE,
                                      rds_path = "RDS_Files/downstream_analyzed.rds",
                                      plot_options = NULL) {

    if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required.")
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")

    `%||%` <- function(a, b) if (!is.null(a)) a else b

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # Dimensionality Reduction and Clustering
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
    seurat_obj <- RunUMAP(seurat_obj, dims = dims)
    seurat_obj <- FindNeighbors(seurat_obj, dims = dims)
    seurat_obj <- FindClusters(seurat_obj, resolution = resolution)

    # Plot: UMAP grouped
    title <- opts$title %||% "UMAP: Grouped by Samples"
    umap_grouped <- DimPlot(seurat_obj, reduction = "umap", group.by = "SampleLabel") +
      ggplot2::ggtitle(title)

    if (is.null(plot_options) || is.null(plot_options$umap_grouped) || !isTRUE(plot_options$umap_grouped$save)) {
      ggsave(file.path(output_dir, "integrated_umap_grouped.jpeg"), umap_grouped, width = 12, height = 12)
    } else {
      opts <- plot_options$umap_grouped
      ggsave(
        filename = file.path(output_dir, paste0(opts$filename %||% "integrated_umap_grouped", ".", opts$format %||% "jpeg")),
        plot = umap_grouped,
        device = opts$format %||% "jpeg",
        width = opts$width %||% 12,
        height = opts$height %||% 12
      )
    }

    # Plot: UMAP split
    title <- opts$title %||% "UMAP: Split by Samples"
    umap_split <- DimPlot(seurat_obj, reduction = "umap", split.by = "SampleLabel") +
      ggplot2::ggtitle(title)

    if (is.null(plot_options) || is.null(plot_options$umap_split) || !isTRUE(plot_options$umap_split$save)) {
      ggsave(file.path(output_dir, "integrated_umap_split.jpeg"), umap_split, width = 18)
    } else {
      opts <- plot_options$umap_split
      ggsave(
        filename = file.path(output_dir, paste0(opts$filename %||% "integrated_umap_split", ".", opts$format %||% "jpeg")),
        plot = umap_split,
        device = opts$format %||% "jpeg",
        width = opts$width %||% 18,
        height = opts$height %||% 12
      )
    }

    # Plot: Cell count pie chart
    data <- data.frame(table(seurat_obj$SampleLabel))
    colnames(data) <- c("SampleLabel", "Count")

    title <- opts$title %||% "Cell Count by Sample"
    pie_plot <- ggplot2::ggplot(data, ggplot2::aes(x = "", y = Count, fill = SampleLabel)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::coord_polar(theta = "y") +
      ggplot2::labs(x = "", y = "", title = title) +
      ggplot2::theme_minimal() +
      ggplot2::geom_text(ggplot2::aes(label = Count), position = ggplot2::position_stack(vjust = 0.5)) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank())

    if (is.null(plot_options) || is.null(plot_options$pie_chart) || !isTRUE(plot_options$pie_chart$save)) {
      ggsave(file.path(output_dir, "cellcount_piechart.jpg"), pie_plot, width = 10)
    } else {
      opts <- plot_options$pie_chart
      ggsave(
        filename = file.path(output_dir, paste0(opts$filename %||% "cellcount_piechart", ".", opts$format %||% "jpeg")),
        plot = pie_plot,
        device = opts$format %||% "jpeg",
        width = opts$width %||% 10,
        height = opts$height %||% 10
      )
    }

    # Save Seurat object
    if (save_rds) {
      dir.create(dirname(rds_path), recursive = TRUE, showWarnings = FALSE)
      saveRDS(seurat_obj, rds_path)
    }

    return(seurat_obj)
  }

