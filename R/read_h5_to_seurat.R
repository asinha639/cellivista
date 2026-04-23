#' Read 10x HDF5 File and Create Seurat Object with Sample Label Integration
#'
#' This function reads a 10x Genomics-formatted .h5 file and returns a Seurat object.
#' It also extracts sample labels from cell barcodes and allows relabeling.
#'
#' @param filepath Path to the .h5 file.
#' @param project_name Optional project name for the Seurat object.
#' @param min_features Minimum number of features per cell. Default is 200.
#' @param min_cells Minimum number of cells per feature. Default is 3.
#' @param add_metadata Optional data.frame of metadata to append to Seurat object.
#' @param sample_label_map Named character vector to relabel sample labels. Default is NULL.
#'
#' @return A Seurat object.
#' @export
#'
#' @examples
#' sample_map <- c("Sample1" = "Control", "Sample2" = "BPV")
#' seurat_obj <- read_h5_to_seurat("data/sample_data.h5", project_name = "PBMC", sample_label_map = sample_map)
read_h5_to_seurat <- function(filepath,
                              project_name = "SeuratProject",
                              min_features = 200,
                              min_cells = 3,
                              add_metadata = NULL,
                              sample_label_map = NULL) {

  # Load required packages
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required but not installed.")
  }
  # if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
  #   stop("SeuratDisk package is required but not installed.")
  # }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required but not installed.")
  }

  library(dplyr)

  # Read the HDF5 file
  message("Reading HDF5 file: ", filepath)
  counts <- Seurat::Read10X_h5(filepath)

  # Create Seurat object
  seurat_obj <- Seurat::CreateSeuratObject(counts = counts,
                                           project = project_name,
                                           min.cells = min_cells,
                                           min.features = min_features)

  # Optional: Add metadata
  if (!is.null(add_metadata)) {
    if (!all(rownames(add_metadata) %in% colnames(seurat_obj))) {
      warning("Some rownames in metadata do not match Seurat cell names. Skipping metadata merge.")
    } else {
      seurat_obj <- Seurat::AddMetaData(seurat_obj, metadata = add_metadata)
    }
  }

  if (!is.null(sample_label_map)) {
    if (!is.character(sample_label_map) || is.null(names(sample_label_map))) {
      stop("sample_label_map must be a named character vector.")
    }
    seurat_obj@meta.data <- seurat_obj@meta.data %>%
      mutate(SampleLabel = recode(SampleLabel, !!!sample_label_map))
  }

  # Parse SampleLabel from cell names (e.g., ending in -1, -2, etc.)
  seurat_obj$SampleLabel <- sub(".*-(\\d+)$", "Sample\\1", rownames(seurat_obj@meta.data))

  # Optionally relabel SampleLabel using user-supplied named vector
  if (!is.null(sample_label_map)) {
    seurat_obj@meta.data <- seurat_obj@meta.data %>%
      mutate(SampleLabel = recode(SampleLabel, !!!sample_label_map))
  }

  return(seurat_obj)
}
