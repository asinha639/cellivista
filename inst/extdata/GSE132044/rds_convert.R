# Load libraries
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(stringr)

# ---------------------------
# 1. Define file paths
# ---------------------------
raw_dir <- "Raw Data/"

mtx_file  <- file.path(raw_dir, "GSE132044_pbmc_hg38_count_matrix.mtx")
cell_file <- file.path(raw_dir, "GSE132044_pbmc_hg38_cell.tsv")
gene_file <- file.path(raw_dir, "GSE132044_pbmc_hg38_gene.tsv")

mtx_gz  <- paste0(mtx_file, ".gz")
cell_gz <- paste0(cell_file, ".gz")
gene_gz <- paste0(gene_file, ".gz")

# ---------------------------
# 2. Decompress files if needed
# ---------------------------
decompress_if_needed <- function(gz_file, output_file) {
  if (file.exists(output_file)) {
    message("Raw file already exists: ", output_file)
  } else if (file.exists(gz_file)) {
    message("Decompressing: ", gz_file)
    R.utils::gunzip(gz_file, destname = output_file, remove = FALSE, overwrite = TRUE)
  } else {
    stop("Missing both raw and compressed file: ", output_file)
  }
}

decompress_if_needed(mtx_gz, mtx_file)
decompress_if_needed(cell_gz, cell_file)
decompress_if_needed(gene_gz, gene_file)

# ---------------------------
# 3. Load data and create Seurat object
# ---------------------------
data <- ReadMtx(
  mtx = mtx_file,
  features = gene_file,
  cells = cell_file,
  feature.column = 1
)

seurat_obj <- CreateSeuratObject(data)

# ---------------------------
# 4. Extract sample metadata
# ---------------------------
meta <- seurat_obj@meta.data

# Extract PBMC sample labels
meta$SampleLabel <- str_extract(rownames(meta), "PBMC\\d+")

# Add to Seurat object
seurat_obj$SampleLabel <- meta$SampleLabel

# ---------------------------
# 5. Save files
# ---------------------------
SaveH5Seurat(seurat_obj, filename = "GSE132044_pbmc_hg38.h5Seurat")
saveRDS(seurat_obj, "GSE132044_pbmc_hg38.rds")

# ---------------------------
# 6. Quick checks
# ---------------------------
table(seurat_obj$SampleLabel, useNA = "ifany")
head(seurat_obj@meta.data)
unique(seurat_obj$SampleLabel)