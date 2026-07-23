#' SCTransform-Based Integration of Multiple Seurat Samples
#'
#' Performs normalization, feature selection, anchor finding, and data integration across multiple Seurat samples.
#'
#' @param seurat_obj A Seurat object with multiple samples labeled in the `SampleLabel` metadata column.
#' @param output_path Path to save the integrated Seurat object (RDS). Default is "RDS_Files/integrated_seurat.obj.rds".
#' @param nfeatures Number of integration features to select. Default is 3000.
#' @param future_max_size Maximum memory allocation for future parallelization. Default is 2GB.
#' @param save_rds Logical; whether to save the integrated object. Default is TRUE.
#' @param output_dir Optional directory to save integration plots. Default is NULL.
#' @param plot_options Optional list of integration plot settings.
#'
#' @return A Seurat object with integrated expression data.
#' @export
#'
#' @examples
#' integrated_obj <- integrate_seurat_samples(seurat_obj_doubletfiltered)
integrate_seurat_samples <- function(seurat_obj,
                                     output_path = "RDS_Files/integrated_seurat.obj.rds",
                                     nfeatures = 3000,
                                     dims = NULL,
                                     k.anchor = NULL,
                                     k.filter = NULL,
                                     k.score = NULL,
                                     k.weight = NULL,
                                     future_max_size = 2 * 1024^3,
                                     save_rds = TRUE,
                                     output_dir = NULL,
                                     plot_options = NULL) {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat package is required.")
  if (!requireNamespace("future", quietly = TRUE)) stop("future package is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 package is required.")
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  if (is.null(seurat_obj) || !inherits(seurat_obj, "Seurat")) {
    stop("integrate_seurat_samples() requires a valid Seurat object.")
  }
  meta <- seurat_obj[[]]
  if (!"SampleLabel" %in% names(meta)) {
    stop("integrate_seurat_samples() requires a 'SampleLabel' column in metadata.")
  }
  if (is.null(SeuratObject::DefaultAssay(seurat_obj))) {
    stop("integrate_seurat_samples() requires a default assay to be set.")
  }
  if (!"RNA" %in% SeuratObject::Assays(seurat_obj)) {
    stop("integrate_seurat_samples() requires an RNA assay.")
  }
  
  # Setup future memory limit
  options(future.globals.maxSize = future_max_size)

  # Seurat v5 layer-aware split before the existing per-sample workflow.
  seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$SampleLabel)
  
  # Split object into a list by SampleLabel
  seurat_list <- Seurat::SplitObject(seurat_obj, split.by = "SampleLabel")
  
  # Normalize each sample using SCTransform
  seurat_list <- lapply(X = seurat_list, FUN = function(x) {
    Seurat::SCTransform(x, verbose = FALSE)
  })

  sample_sizes <- vapply(seurat_list, ncol, integer(1))
  smallest_sample_size <- min(sample_sizes)
  if (smallest_sample_size < 3L) {
    stop("integrate_seurat_samples() requires at least 3 cells in each sample.")
  }

  max_supported_dims <- max(1L, smallest_sample_size - 1L)
  max_supported_dims <- min(max_supported_dims, min(vapply(seurat_list, nrow, integer(1))) - 1L)
  max_supported_dims <- max(1L, as.integer(max_supported_dims))
  default_dims <- seq_len(min(30L, max_supported_dims))

  dims <- if (is.null(dims)) {
    default_dims
  } else {
    dims <- as.integer(dims)
    dims <- dims[!is.na(dims) & dims > 0]
    unique(dims)
  }
  dims <- dims[dims <= max_supported_dims]
  if (length(dims) == 0) {
    stop("No valid PCA dimensions remain after adjusting for the smallest sample size.")
  }

  resolve_k <- function(value, default_value, max_value) {
    if (is.null(value)) {
      value <- default_value
    }
    value <- as.integer(value)
    value <- value[!is.na(value) & value > 0]
    if (length(value) == 0) {
      stop("Integration neighborhood settings must be positive integers.")
    }
    value <- value[[1L]]
    max_value <- as.integer(max_value)
    if (is.na(max_value) || max_value < 1L) {
      stop("integrate_seurat_samples() requires at least 3 cells in each sample.")
    }
    min(value, max_value)
  }

  max_anchor_neighbors <- max(2L, smallest_sample_size - 1L)
  max_filter_neighbors <- max(2L, smallest_sample_size - 1L)
  max_score_neighbors <- max(2L, smallest_sample_size - 1L)
  max_weight_neighbors <- max(2L, smallest_sample_size - 1L)

  k.anchor <- resolve_k(k.anchor, 5L, max_anchor_neighbors)
  k.filter <- resolve_k(k.filter, 200L, max_filter_neighbors)
  k.score <- resolve_k(k.score, 30L, max_score_neighbors)
  k.weight <- resolve_k(k.weight, 100L, max_weight_neighbors)
  
  # Feature selection
  features <- Seurat::SelectIntegrationFeatures(seurat_list, nfeatures = nfeatures)
  
  # Prepare for SCT integration
  seurat_list <- Seurat::PrepSCTIntegration(seurat_list, anchor.features = features)
  
  # Find integration anchors
  anchors <- Seurat::FindIntegrationAnchors(object.list = seurat_list,
                                            normalization.method = "SCT",
                                            anchor.features = features,
                                            dims = dims,
                                            nn.method = "rann",
                                            k.anchor = k.anchor,
                                            k.filter = k.filter,
                                            k.score = k.score)
  
  # Integrate
  integrated <- Seurat::IntegrateData(
    anchorset = anchors,
    normalization.method = "SCT",
    k.weight = k.weight,
    dims = dims
  )

  # Generate optional integration plots
  if (!is.null(plot_options) && length(plot_options) > 0) {
    output_dir <- output_dir %||% dirname(output_path)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    SeuratObject::DefaultAssay(integrated) <- "integrated"
    integrated <- Seurat::ScaleData(integrated, verbose = FALSE)
    integrated <- Seurat::RunPCA(integrated, verbose = FALSE)

    available_dims <- tryCatch({
      ncol(Seurat::Embeddings(integrated, reduction = "pca"))
    }, error = function(e) 0L)
    dims_use <- seq_len(min(15, as.integer(available_dims)))
    if (length(dims_use) == 0) {
      stop("Failed to compute PCA dimensions required for integration UMAP plots.")
    }
    integrated <- Seurat::RunUMAP(integrated, dims = dims_use, verbose = FALSE)

    save_plot <- function(plot_obj, opts = NULL, default_filename, default_title,
                          default_format = "jpeg", default_width = 10, default_height = 8) {
      opts <- opts %||% list()
      plot_obj <- plot_obj + ggplot2::ggtitle(opts$title %||% default_title)
      format <- tolower(opts$format %||% default_format)
      device <- switch(
        format,
        jpg = "jpeg",
        jpeg = "jpeg",
        png = "png",
        svg = "svg",
        format
      )
      filename <- paste0(
        tools::file_path_sans_ext(opts$filename %||% default_filename),
        ".",
        format
      )
      path <- file.path(output_dir, filename)
      ggplot2::ggsave(
        filename = path,
        plot = plot_obj,
        device = device,
        width = as.numeric(opts$width %||% default_width),
        height = as.numeric(opts$height %||% default_height)
      )
      invisible(path)
    }

    if (!is.null(plot_options$umap_grouped) && isTRUE(plot_options$umap_grouped$save)) {
      umap_grouped <- Seurat::DimPlot(
        integrated,
        reduction = "umap",
        group.by = "SampleLabel"
      )
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

    if (!is.null(plot_options$umap_split) && isTRUE(plot_options$umap_split$save)) {
      umap_split <- Seurat::DimPlot(
        integrated,
        reduction = "umap",
        split.by = "SampleLabel"
      )
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

    if (!is.null(plot_options$pie_chart) && isTRUE(plot_options$pie_chart$save)) {
      pie_data <- as.data.frame(table(integrated$SampleLabel), stringsAsFactors = FALSE)
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
      save_plot(
        plot_obj = pie_plot,
        opts = plot_options$pie_chart,
        default_filename = "integration_cellcount_piechart",
        default_title = "Cell Count by Sample",
        default_format = "jpeg",
        default_width = 10,
        default_height = 10
      )
    }
  }
  
  # Save
  if (save_rds) {
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(integrated, output_path)
  }
  
  return(integrated)
}
