#' Doublet Removal Using DoubletFinder per SampleLabel
#'
#' Performs doublet detection and removal for each sample in a Seurat object, using DoubletFinder.
#'
#' @param seurat_obj A Seurat object with a `SampleLabel` column in metadata.
#' @param output_dir Directory to save plots. Default is "./Plots/DoubletFinder".
#' @param save_rds Logical; whether to save the final filtered Seurat object. Default is TRUE.
#' @param rds_path Path to save filtered Seurat object. Default is "seurat.obj_doubletfiltered.rds".
#' @param doublet_rate Expected doublet rate (e.g., 0.075 for 7.5%). Default is 0.075.
#'
#' @return A Seurat object with doublets removed.
#' @export
remove_doublets <- function(seurat_obj,
                            output_dir = "./Plots/DoubletFinder",
                            save_rds = TRUE,
                            rds_path = "RDS_Files/seurat.obj_doubletfiltered.rds",
                            doublet_rate = 0.075) {

  if (!requireNamespace("DoubletFinder", quietly = TRUE)) stop("DoubletFinder is required.")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  library(ggplot2)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  sample_labels <- unique(seurat_obj$SampleLabel)
  singlet_list <- list()
  count_table <- list()

  for (sample in sample_labels) {
    message("Processing sample: ", sample)

    sample_seurat <- subset(seurat_obj, subset = SampleLabel == sample)
    sample_seurat <- NormalizeData(sample_seurat)
    sample_seurat <- FindVariableFeatures(sample_seurat, selection.method = "vst", nfeatures = 3000)
    sample_seurat <- ScaleData(sample_seurat)
    sample_seurat <- RunPCA(sample_seurat)
    sample_seurat <- RunUMAP(sample_seurat, dims = 1:10)

    sweep.res.list <- DoubletFinder::paramSweep(sample_seurat, PCs = 1:15, sct = FALSE)
    sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- DoubletFinder::find.pK(sweep.stats)

    pK <- bcmvn %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::pull(pK) %>% as.character() %>% as.numeric()

    annotations <- sample_seurat@meta.data$seurat_clusters
    homotypic.prop <- DoubletFinder::modelHomotypic(annotations)
    nExp_poi <- round(doublet_rate * nrow(sample_seurat@meta.data))
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

    sample_seurat <- DoubletFinder::doubletFinder(sample_seurat,
                                                  PCs = 1:15, pN = 0.25, pK = pK, nExp = nExp_poi,
                                                  reuse.pANN = FALSE, sct = FALSE)

    df_col <- colnames(sample_seurat@meta.data)[grepl("DF.classifications", colnames(sample_seurat@meta.data))]

    # Save UMAP
    df_plot <- DimPlot(sample_seurat, reduction = "umap", group.by = df_col) + ggtitle(paste(sample, "DoubletFinder"))
    ggsave(file.path(output_dir, paste0("doubletfinder_", sample, ".jpeg")), df_plot, width = 12, height = 12)

    # Filter singlets
    singlet_cells <- rownames(sample_seurat@meta.data)[sample_seurat@meta.data[[df_col]] == "Singlet"]
    singlet_list[[sample]] <- subset(sample_seurat, cells = singlet_cells)

    # Count stats
    count_table[[sample]] <- table(sample_seurat@meta.data[[df_col]])
  }

  # Convert list of tables to a single tidy dataframe
  df_count <- purrr::map2_dfr(
    count_table,
    names(count_table),
    ~ data.frame(Sample = .y, Category = names(.x), Count = as.integer(.x))
  )

  # Calculate percentages
  df_long <- df_count %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(Percentage = Count / sum(Count)) %>%
    dplyr::ungroup()

  # Plot
  bar_plot <- ggplot(df_long, aes(x = Sample, y = Count, fill = Category)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = paste0(Count, " (", round(Percentage * 100, 1), "%)")),
              position = position_dodge(width = 0.9), vjust = -0.5, size = 3,
              fontface = "bold", color = "gray35") +
    labs(title = "Singlet & Doublet Count per Sample", x = "Sample", y = "Count", fill = "Category") +
    theme_minimal()

  # Save the plot
  ggsave(file.path(output_dir, "singlet.doublet_count_barplot.jpeg"), bar_plot, width = 8, height = 6)

  # Combine all singlet cell names across samples
  all_singlet_cells <- unlist(lapply(singlet_list, function(s) colnames(s)))

  # Subset original Seurat object to only retain singlets
  singlet_merged <- subset(seurat_obj, cells = all_singlet_cells)

  # Optional save
  if (save_rds) {
    saveRDS(singlet_merged, rds_path)
  }

  return(singlet_merged)
}
