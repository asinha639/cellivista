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
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required but not installed.")
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Filter Seurat object based on QC metrics
  seurat_obj <- subset(seurat_obj,
                       subset = nFeature_RNA > min_features &
                         nFeature_RNA < max_features &
                         percent.mt < max_mt_percent)
  
  # Post-QC violin plot
  post_qc_vplot <- Seurat::VlnPlot(
    seurat_obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    split.by = sample_split_var,
    ncol = 3
  ) & ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    legend.position = "right"
  )
  ggplot2::ggsave(filename = file.path(output_dir, "post_qc_vplot.jpg"),
                  plot = post_qc_vplot, width = 18, height = 12)
  
  # Density plot for nCount_RNA
  density_plot <- seurat_obj@meta.data %>%
    ggplot2::ggplot(ggplot2::aes(color = .data[[sample_split_var]], 
                                 x = nCount_RNA, 
                                 fill = .data[[sample_split_var]])) +
    ggplot2::geom_density(alpha = 0.2) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_log10()
  ggplot2::ggsave(file.path(output_dir, "post_qc_densityplot.jpg"),
                  plot = density_plot, width = 18, height = 12)
  
  # Pie chart showing cell count distribution by sample
  pie_data <- data.frame(table(seurat_obj@meta.data[[sample_split_var]]))
  colnames(pie_data) <- c("SampleLabel", "Count")
  
  cellcount_piechart <- ggplot2::ggplot(pie_data, ggplot2::aes(x = "", y = Count, fill = SampleLabel)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::labs(title = "Cell Count by SampleLabel", x = "", y = "") +
    ggplot2::theme_minimal() +
    ggplot2::geom_text(ggplot2::aes(label = Count),
                       position = ggplot2::position_stack(vjust = 0.5)) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
  
  ggplot2::ggsave(file.path(output_dir, "cellcount_piechart_post.jpg"), cellcount_piechart)
  
  # Optional: Save Seurat object
  if (save_rds) {
    rds_dir <- dirname(rds_path)
    if (!dir.exists(rds_dir)) dir.create(rds_dir, recursive = TRUE)
    saveRDS(seurat_obj, rds_path)
  }
  
  return(seurat_obj)
}
