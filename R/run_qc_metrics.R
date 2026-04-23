#' Perform Quality Control and Generate QC Plots
#'
#' Calculates mitochondrial read percentages, creates standard violin and scatter plots,
#' and optionally saves them to a specified directory.
#'
#' @param seurat_obj A Seurat object.
#' @param mt_pattern Regex pattern to identify mitochondrial genes. Default is "^MT-".
#' @param output_dir Directory path to save output plots. Default is "./Plots".
#' @param sample_split_var Column in metadata to split violin plots. Default is "SampleLabel".
#' @param ylimit Violin plot Y-axis limits for percent.mt. Default is c(0, 80).
#' @param save_plots Logical. Whether to save the generated plots. Default is TRUE.
#'
#' @return A Seurat object with updated QC metrics in metadata.
#' @export
#'
#' @examples
#' seurat_obj <- run_qc_metrics(seurat_obj, mt_pattern = "^mt-", output_dir = "QC_Plots/", save_plots = FALSE)
run_qc_metrics <- function(seurat_obj,
                           mt_pattern = "^MT-",
                           output_dir = "./Plots",
                           sample_split_var = "SampleLabel",
                           ylimit = c(0, 80),
                           save_plots = TRUE) {

  if (!requireNamespace("Seurat", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Required packages Seurat and ggplot2 are not installed.")
  }

  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Detect mitochondrial genes and calculate percent.mt
  seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = mt_pattern)

  # Violin plots for QC metrics
  qc_violin <- Seurat::VlnPlot(
    object = seurat_obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    split.by = sample_split_var,
    ncol = 3
  ) & ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    legend.position = "right"
  )

  if (save_plots) {
    ggplot2::ggsave(
      filename = file.path(output_dir, "pre_qc_violin.jpg"),
      plot = print(qc_violin),
      width = 18, height = 12
    )
  }

  # Violin plot with Y-limit for percent.mt
  qc_mt <- Seurat::VlnPlot(
    object = seurat_obj,
    features = "percent.mt",
    split.by = sample_split_var
  ) + ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  ) + ggplot2::ylim(ylimit[1], ylimit[2])

  if (save_plots) {
    ggplot2::ggsave(
      filename = file.path(output_dir, "pre_qc_percent.mt.jpg"),
      plot = print(qc_mt),
      width = 10, height = 6
    )
  }

  # Feature scatter plots
  plot1 <- Seurat::FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- Seurat::FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

  combined_plot <- patchwork::wrap_plots(plot1, plot2)

  if (save_plots) {
    ggplot2::ggsave(
      filename = file.path(output_dir, "pre_qc_FeatureScatterplot.jpg"),
      plot = print(combined_plot),
      width = 18, height = 12
    )
  }

  return(seurat_obj)
}
