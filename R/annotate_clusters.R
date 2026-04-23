#' Annotate Clusters with User-Defined Cell Types and Visualize
#'
#' @param seurat_obj A clustered Seurat object (with `seurat_clusters`).
#' @param annotations A character vector of annotation labels, or a CSV path.
#' @param output_dir Directory to save plots. Default is "./Plots/Annotation".
#' @param rds_path Path to save annotated object. Default is "RDS_Files/annotated_seurat.obj.rds".
#' @param save_rds Logical; whether to save the annotated Seurat object. Default is TRUE.
#'
#' @return Annotated Seurat object.
#' @export
annotate_clusters <- function(seurat_obj,
                              annotations,
                              output_dir = "./Plots/Annotation",
                              rds_path = "RDS_Files/annotated_seurat.obj.rds",
                              save_rds = TRUE) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Handle CSV or vector
  if (is.character(annotations) && file.exists(annotations)) {
    cluster_annot <- readr::read_csv(annotations, col_names = c("seurat_clusters", "cell_annotation"))
    cluster_annot$seurat_clusters <- as.factor(cluster_annot$seurat_clusters)
  } else if (is.character(annotations) && length(annotations) == length(unique(seurat_obj$seurat_clusters))) {
    cluster_annot <- data.frame(
      seurat_clusters = as.factor(0:(length(annotations) - 1)),
      cell_annotation = annotations
    )
  } else {
    stop("Invalid annotations: must be a path to CSV or a character vector matching number of clusters.")
  }

  # Annotate
  temp_rowname <- rownames(seurat_obj@meta.data)
  seurat_obj@meta.data <- dplyr::left_join(seurat_obj@meta.data, cluster_annot, by = "seurat_clusters")
  seurat_obj@meta.data$cell_annotation <- as.character(seurat_obj@meta.data$cell_annotation)
  rownames(seurat_obj@meta.data) <- temp_rowname
  seurat_obj$cell_annotation <- factor(seurat_obj$cell_annotation, levels = unique(cluster_annot$cell_annotation))
  Seurat::Idents(seurat_obj) <- seurat_obj$cell_annotation

  # UMAP plots
  umap_annotated <- Seurat::DimPlot(seurat_obj, reduction = "umap", label = F, group.by = "cell_annotation") +
    ggplot2::ggtitle("UMAP: Annotated Cell Types")

  ggsave(file.path(output_dir, "umap_annotated.jpeg"), umap_annotated, width = 10, height = 8)

  umap_split <- Seurat::DimPlot(seurat_obj, reduction = "umap", split.by = "SampleLabel", group.by = "cell_annotation",
                                label = F) + ggplot2::ggtitle(NULL)  # Remove default title
  ggsave(file.path(output_dir, "umap_split_by_sample.jpeg"), umap_split, width = 18, height = 8)

  # Pie chart of cell type distribution
  pie_data <- data.frame(table(seurat_obj$cell_annotation))
  colnames(pie_data) <- c("CellType", "Count")

  pie_plot <- ggplot2::ggplot(pie_data, ggplot2::aes(x = "", y = Count, fill = CellType)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::theme_minimal() +
    ggplot2::geom_text(ggplot2::aes(label = Count), position = ggplot2::position_stack(vjust = 0.5)) +
    ggplot2::labs(title = "Cell Type Composition", x = "", y = "") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())

  ggsave(file.path(output_dir, "celltype_piechart.jpeg"), pie_plot)

  # Save annotated Seurat object
  if (save_rds) {
    dir.create(dirname(rds_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(seurat_obj, rds_path)
  }

  return(seurat_obj)
}

