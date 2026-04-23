#' Run Pairwise DGEA Between Clusters or Sample Labels
#'
#' @param seurat_obj A Seurat object.
#' @param group_by Metadata column to group by (e.g., "SampleLabel", "seurat_clusters").
#' @param ident1 First group value (e.g., "Control" or 0).
#' @param ident2 Second group value (e.g., "Trt1" or 2).
#' @param output_path Path to save CSV results. Default is "pairwise_dgea.csv".
#' @param logfc_threshold Log fold change threshold. Default is 0.25.
#' @param min_pct Minimum percentage expression. Default is 0.1.
#'
#' @return A data frame with DGEA results.
#' @export
run_pairwise_dgea <- function(seurat_obj,
                              group_by = "SampleLabel",
                              ident1,
                              ident2,
                              output_path = "pairwise_dgea.csv",
                              logfc_threshold = 0.25,
                              min_pct = 0.1) {
  
  Seurat::DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- Seurat::SetIdent(seurat_obj, value = group_by)
  
  dgea_result <- Seurat::FindMarkers(
    seurat_obj,
    ident.1 = ident1,
    ident.2 = ident2,
    logfc.threshold = logfc_threshold,
    min.pct = min_pct
  )
  
  dgea_result$gene <- rownames(dgea_result)
  readr::write_csv(dgea_result, output_path)
  
  return(dgea_result)
}
