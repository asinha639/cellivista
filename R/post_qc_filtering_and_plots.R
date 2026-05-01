#' Filter Seurat Object Based on QC Metrics and Generate Post-QC Plots
#'
#' Applies user-defined QC filters to a Seurat object, generates post-QC plots, and optionally saves the object as an RDS file.
#'
#' @param seurat_obj A Seurat object with QC metrics already calculated.
#' @param min_features Minimum number of genes per cell. Default is 200.
#' @param max_features Maximum number of genes per cell. Default is 7500.
#' @param max_mt_percent Maximum mitochondrial content. Default is 25.
#' @param output_dir Directory path to save output plots. Default is "./Plots".
#' @param sample_split_var Column in metadata to split plots. Default is "SampleLabel".
#' @param save_rds Logical; whether to save the filtered Seurat object. Default is TRUE.
#' @param rds_path Path to save the RDS file. Default is "RDS_Files/seurat.obj_qc.rds".
#'
#' @return Filtered Seurat object.
#' @export
#'
#' @examples
#' seurat_obj <- post_qc_filtering_and_plots(seurat_obj)
post_qc_filtering_and_plots <- function(seurat_obj,
                                        min_features = 200,
                                        max_features = 7500,
                                        max_mt_percent = 25,
                                        output_dir = "./Plots",
                                        sample_split_var = "SampleLabel",
                                        save_rds = TRUE,
                                        rds_path = "RDS_Files/seurat.obj_qc.rds") {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required but not installed.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required but not installed.")
  }
  
  if (is.null(seurat_obj)) {
    stop("seurat_obj is NULL.")
  }
  
  required_cols <- c("nFeature_RNA", "nCount_RNA", "percent.mt", sample_split_var)
  missing_cols <- setdiff(required_cols, colnames(seurat_obj@meta.data))
  if (length(missing_cols) > 0) {
    stop(
      "seurat_obj is missing required metadata columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (!is.numeric(min_features) || !is.numeric(max_features) || !is.numeric(max_mt_percent)) {
    stop("min_features, max_features, and max_mt_percent must be numeric.")
  }
  if (min_features >= max_features) {
    stop("min_features must be less than max_features.")
  }
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  seurat_obj <- subset(
    x = seurat_obj,
    subset = nFeature_RNA > min_features &
      nFeature_RNA < max_features &
      percent.mt < max_mt_percent
  )
  
  if (ncol(seurat_obj) == 0) {
    stop("Filtering removed all cells. Please relax QC thresholds.")
  }
  
  post_qc_vplot <- Seurat::VlnPlot(
    object = seurat_obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    split.by = sample_split_var,
    ncol = 3
  ) & ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    legend.position = "right"
  )
  
  ggplot2::ggsave(
    filename = file.path(output_dir, "post_qc_vplot.jpg"),
    plot = post_qc_vplot,
    width = 18,
    height = 12
  )
  
  meta_df <- seurat_obj@meta.data
  
  density_plot <- ggplot2::ggplot(
    meta_df,
    ggplot2::aes(
      x = .data$nCount_RNA,
      color = .data[[sample_split_var]],
      fill = .data[[sample_split_var]]
    )
  ) +
    ggplot2::geom_density(alpha = 0.2) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = "Post-QC nCount_RNA Density",
      x = "nCount_RNA",
      y = "Density",
      color = sample_split_var,
      fill = sample_split_var
    )
  
  ggplot2::ggsave(
    filename = file.path(output_dir, "post_qc_densityplot.jpg"),
    plot = density_plot,
    width = 18,
    height = 12
  )
  
  pie_data <- as.data.frame(table(meta_df[[sample_split_var]]), stringsAsFactors = FALSE)
  colnames(pie_data) <- c("SampleLabel", "Count")
  
  cellcount_piechart <- ggplot2::ggplot(
    pie_data,
    ggplot2::aes(x = "", y = Count, fill = SampleLabel)
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::labs(title = "Cell Count by SampleLabel", x = "", y = "") +
    ggplot2::theme_minimal() +
    ggplot2::geom_text(
      ggplot2::aes(label = Count),
      position = ggplot2::position_stack(vjust = 0.5)
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
  
  ggplot2::ggsave(
    filename = file.path(output_dir, "cellcount_piechart_post.jpg"),
    plot = cellcount_piechart,
    width = 10,
    height = 10
  )
  
  if (isTRUE(save_rds)) {
    dir.create(dirname(rds_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(seurat_obj, rds_path)
  }
  
  return(seurat_obj)
}