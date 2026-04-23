#' SCTransform-Based Integration of Multiple Seurat Samples
#'
#' Performs normalization, feature selection, anchor finding, and data integration across multiple Seurat samples.
#'
#' @param seurat_obj A Seurat object with multiple samples labeled in the `SampleLabel` metadata column.
#' @param output_path Path to save the integrated Seurat object (RDS). Default is "RDS_Files/integrated_seurat.obj.rds".
#' @param nfeatures Number of integration features to select. Default is 3000.
#' @param future_max_size Maximum memory allocation for future parallelization. Default is 2GB.
#' @param save_rds Logical; whether to save the integrated object. Default is TRUE.
#'
#' @return A Seurat object with integrated expression data.
#' @export
#'
#' @examples
#' integrated_obj <- integrate_seurat_samples(seurat_obj_doubletfiltered)
integrate_seurat_samples <- function(seurat_obj,
                                     output_path = "RDS_Files/integrated_seurat.obj.rds",
                                     nfeatures = 3000,
                                     future_max_size = 2 * 1024^3,
                                     save_rds = TRUE) {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat package is required.")
  if (!requireNamespace("future", quietly = TRUE)) stop("future package is required.")
  
  # Setup future memory limit
  options(future.globals.maxSize = future_max_size)
  
  # Split object into a list by SampleLabel
  seurat_list <- Seurat::SplitObject(seurat_obj, split.by = "SampleLabel")
  
  # Normalize each sample using SCTransform
  seurat_list <- lapply(X = seurat_list, FUN = function(x) {
    Seurat::SCTransform(x, verbose = FALSE)
  })
  
  # Feature selection
  features <- Seurat::SelectIntegrationFeatures(seurat_list, nfeatures = nfeatures)
  
  # Prepare for SCT integration
  seurat_list <- Seurat::PrepSCTIntegration(seurat_list, anchor.features = features)
  
  # Find integration anchors
  anchors <- Seurat::FindIntegrationAnchors(object.list = seurat_list,
                                            normalization.method = "SCT",
                                            anchor.features = features)
  
  # Integrate
  integrated <- Seurat::IntegrateData(anchorset = anchors, normalization.method = "SCT")
  
  # Save
  if (save_rds) {
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(integrated, output_path)
  }
  
  return(integrated)
}
