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
#' @param save_plots Logical; whether downstream plots should be saved. Default is TRUE.
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
                                    save_plots = TRUE,
                                    plot_options = NULL) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required.")
  }
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  if (is.null(seurat_obj)) {
    stop("seurat_obj is NULL.")
  }
  if (!"SampleLabel" %in% colnames(seurat_obj@meta.data)) {
    stop("seurat_obj must contain a 'SampleLabel' column in metadata.")
  }
  
  active_assay <- SeuratObject::DefaultAssay(seurat_obj)
  if (is.null(active_assay) || !nzchar(active_assay)) {
    stop("seurat_obj must have an active assay before downstream analysis.")
  }
  data_layers <- SeuratObject::Layers(seurat_obj, assay = active_assay, search = "data")
  if (length(data_layers) == 0) {
    seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
  }
  if (length(Seurat::VariableFeatures(seurat_obj, assay = active_assay)) == 0) {
    seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, verbose = FALSE)
  }
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  dims <- as.integer(dims)
  dims <- dims[!is.na(dims) & dims > 0]
  if (length(dims) == 0) {
    stop("dims must contain at least one positive integer.")
  }
  
  # Core downstream workflow
  seurat_obj <- Seurat::ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- Seurat::RunPCA(seurat_obj, verbose = FALSE)
  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = dims, verbose = FALSE)
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = dims, verbose = FALSE)
  seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = resolution, verbose = FALSE)
  
  # Helper to save plots consistently
  save_plot <- function(plot_obj, opts = NULL, default_filename, default_title,
                        default_format = "jpeg", default_width = 10, default_height = 8) {
    opts <- opts %||% list()
    
    plot_obj <- plot_obj + ggplot2::ggtitle(opts$title %||% default_title)
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
    
    ggplot2::ggsave(
      filename = file.path(output_dir, filename),
      plot = plot_obj,
      device = device,
      width = as.numeric(opts$width %||% default_width),
      height = as.numeric(opts$height %||% default_height)
    )
    
    invisible(file.path(output_dir, filename))
  }
  
  if (isTRUE(save_plots)) {
    # Grouped UMAP
    umap_grouped <- Seurat::DimPlot(
      seurat_obj,
      reduction = "umap",
      group.by = "SampleLabel"
    )
    
    if (is.null(plot_options)) {
      save_plot(
        plot_obj = umap_grouped,
        default_filename = "integrated_umap_grouped",
        default_title = "UMAP: Grouped by Samples",
        default_format = "jpeg",
        default_width = 12,
        default_height = 12
      )
    } else if (!is.null(plot_options$umap_grouped) && isTRUE(plot_options$umap_grouped$save)) {
      save_plot(
        plot_obj = umap_grouped,
        opts = plot_options$umap_grouped,
        default_filename = "integrated_umap_grouped",
        default_title = "UMAP: Grouped by Samples",
        default_format = "jpeg",
        default_width = 12,
        default_height = 12
      )
    }
    
    # Split UMAP
    umap_split <- Seurat::DimPlot(
      seurat_obj,
      reduction = "umap",
      split.by = "SampleLabel"
    )
    
    if (is.null(plot_options)) {
      save_plot(
        plot_obj = umap_split,
        default_filename = "integrated_umap_split",
        default_title = "UMAP: Split by Samples",
        default_format = "jpeg",
        default_width = 18,
        default_height = 12
      )
    } else if (!is.null(plot_options$umap_split) && isTRUE(plot_options$umap_split$save)) {
      save_plot(
        plot_obj = umap_split,
        opts = plot_options$umap_split,
        default_filename = "integrated_umap_split",
        default_title = "UMAP: Split by Samples",
        default_format = "jpeg",
        default_width = 18,
        default_height = 12
      )
    }
    
    # Cell count pie chart
    pie_data <- as.data.frame(table(seurat_obj$SampleLabel), stringsAsFactors = FALSE)
    colnames(pie_data) <- c("SampleLabel", "Count")
    
    pie_plot <- ggplot2::ggplot(
      pie_data,
      ggplot2::aes(x = "", y = Count, fill = SampleLabel)
    ) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::coord_polar(theta = "y") +
      ggplot2::labs(x = "", y = "") +
      ggplot2::theme_minimal() +
      ggplot2::geom_text(
        ggplot2::aes(label = Count),
        position = ggplot2::position_stack(vjust = 0.5)
      ) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank())
    
    if (is.null(plot_options)) {
      save_plot(
        plot_obj = pie_plot,
        default_filename = "cellcount_piechart",
        default_title = "Cell Count by Sample",
        default_format = "jpeg",
        default_width = 10,
        default_height = 10
      )
    } else if (!is.null(plot_options$pie_chart) && isTRUE(plot_options$pie_chart$save)) {
      save_plot(
        plot_obj = pie_plot,
        opts = plot_options$pie_chart,
        default_filename = "cellcount_piechart",
        default_title = "Cell Count by Sample",
        default_format = "jpeg",
        default_width = 10,
        default_height = 10
      )
    }
  }
  
  if (save_rds) {
    dir.create(dirname(rds_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(seurat_obj, rds_path)
  }
  
  return(seurat_obj)
}
