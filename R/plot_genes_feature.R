#' Plot Gene Expression for Display in Shiny and Optional Saving
#'
#' @param seurat_obj A Seurat object.
#' @param genes Character vector of gene names to plot.
#' @param cell_type Optional. Exact name of the cell type in `cell_annotation` to subset. Default is NULL (no filtering).
#' @param split_by Optional. Metadata column name to split FeaturePlot. Default is NULL (no splitting).
#' @param output_dir Directory to save plots. Default is "./Plots/FeaturePlots".
#' @param save_plots Logical; whether to save plots. Default is TRUE.
#' @param plot_options Optional list controlling saved feature plot exports:
#'        filename_prefix, title_prefix, format, width, height.
#'
#' @return A named list of ggplot objects, one per gene (for use in Shiny display).
#' @export
#'
#' @examples
#' plot_genes_feature(seurat_obj, genes = c("GeneA", "GeneB"))
plot_genes_feature <- function(seurat_obj,
                               genes,
                               cell_type = NULL,
                               split_by = NULL,
                               output_dir = "./Plots/FeaturePlots",
                               save_plots = TRUE,
                               plot_options = NULL) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required but not installed.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required but not installed.")
  }
  
  if (is.null(seurat_obj)) {
    stop("seurat_obj is NULL.")
  }
  
  if (length(genes) == 1 && is.character(genes) && grepl(",", genes, fixed = TRUE)) {
    genes <- trimws(strsplit(genes, ",", fixed = TRUE)[[1]])
  }
  genes <- unique(trimws(as.character(genes)))
  genes <- genes[nzchar(genes)]
  
  if (length(genes) == 0) {
    stop("No valid genes were provided.")
  }
  
  if (!is.null(split_by) && !split_by %in% colnames(seurat_obj@meta.data)) {
    stop("split_by not found in seurat_obj@meta.data: ", split_by)
  }
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (!is.null(cell_type)) {
    if (!"cell_annotation" %in% colnames(seurat_obj@meta.data)) {
      stop("The Seurat object does not have a 'cell_annotation' column.")
    }
    if (!cell_type %in% unique(as.character(seurat_obj$cell_annotation))) {
      stop("Cell type not found in 'cell_annotation': ", cell_type)
    }
    
    seurat_obj <- Seurat::subset(seurat_obj, subset = cell_annotation == cell_type)
    
    output_dir <- file.path(output_dir, gsub("[^A-Za-z0-9_\\-]", "_", cell_type))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  available_genes <- rownames(seurat_obj)
  valid_genes <- intersect(genes, available_genes)
  missing_genes <- setdiff(genes, available_genes)
  
  if (length(valid_genes) == 0) {
    stop("None of the requested genes were found in the Seurat object.")
  }
  if (length(missing_genes) > 0) {
    warning("Skipping missing gene(s): ", paste(missing_genes, collapse = ", "))
  }
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  sanitize_file_component <- function(x) gsub("[^A-Za-z0-9_\\-]", "_", x)
  use_custom_export <- !is.null(plot_options)
  plot_options <- plot_options %||% list()
  filename_prefix <- plot_options$filename_prefix %||% "featureplot_"
  title_prefix <- plot_options$title_prefix %||% ""
  format <- tolower(plot_options$format %||% "jpeg")
  device <- switch(format,
                   jpg = "jpeg",
                   jpeg = "jpeg",
                   png = "png",
                   svg = "svg",
                   format)
  width <- as.numeric(plot_options$width %||% 7)
  height <- as.numeric(plot_options$height %||% 6)
  filename_prefix <- tools::file_path_sans_ext(filename_prefix)
  
  plot_list <- list()
  
  for (gene in valid_genes) {
    p <- Seurat::FeaturePlot(
      seurat_obj,
      features = gene,
      split.by = split_by
    ) + ggplot2::ggtitle(if (nzchar(title_prefix)) paste(title_prefix, gene) else gene)
    
    if (isTRUE(save_plots)) {
      if (use_custom_export) {
        ggplot2::ggsave(
          filename = file.path(output_dir, paste0(filename_prefix, sanitize_file_component(gene), ".", format)),
          plot = p,
          device = device,
          width = width,
          height = height
        )
      } else {
        ggplot2::ggsave(
          filename = file.path(output_dir, paste0(gene, "_featureplot.jpeg")),
          plot = p,
          width = 7,
          height = 6
        )
      }
    }
    
    plot_list[[gene]] <- p
  }
  
  return(plot_list)
}
