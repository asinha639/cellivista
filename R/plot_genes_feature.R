#' Plot Gene Expression for Display in Shiny and Optional Saving
#'
#' @param seurat_obj A Seurat object.
#' @param genes Character vector of gene names to plot.
#' @param cell_type Optional. Exact name of the cell type in `cell_annotation` to subset. Default is NULL (no filtering).
#' @param split_by Optional. Metadata column name to split FeaturePlot. Default is NULL (no splitting).
#' @param output_dir Directory to save plots. Default is "./Plots/FeaturePlots".
#' @param save_plots Logical; whether to save plots. Default is TRUE.
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
                               save_plots = TRUE) {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat package is required but not installed.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 package is required but not installed.")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Subset Seurat object if cell_type is specified
  if (!is.null(cell_type)) {
    if (!"cell_annotation" %in% colnames(seurat_obj@meta.data)) {
      stop("The Seurat object does not have a 'cell_annotation' column.")
    }
    if (!cell_type %in% unique(seurat_obj$cell_annotation)) {
      stop(paste("Cell type", cell_type, "not found in 'cell_annotation' column."))
    }
    seurat_obj <- subset(seurat_obj, subset = cell_annotation == cell_type)
    message("Subsetting to cell type: ", cell_type)
    
    output_dir <- file.path(output_dir, gsub(" ", "_", cell_type))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  plot_list <- list()
  
  for (gene in genes) {
    if (!gene %in% rownames(seurat_obj)) {
      warning(paste("Gene", gene, "not found in Seurat object. Skipping."))
      next
    }
    
    p <- Seurat::FeaturePlot(seurat_obj, features = gene, split.by = split_by) +
      ggplot2::ggtitle(gene)
    
    if (save_plots) {
      ggsave(
        filename = file.path(output_dir, paste0(gene, "_featureplot.jpeg")),
        plot = p,
        width = 7, height = 6
      )
    }
    
    plot_list[[gene]] <- p
  }
  
  return(plot_list) 
}
