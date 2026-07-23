make_toy_counts <- function(n_genes = 8,
                            n_cells = 6,
                            mt_genes = 2,
                            structured = FALSE) {
  stopifnot(n_genes >= 4, n_cells >= 4, mt_genes >= 1, mt_genes < n_genes)

  set.seed(101)
  counts <- matrix(
    rpois(n_genes * n_cells, lambda = 5),
    nrow = n_genes,
    ncol = n_cells
  )

  if (structured) {
    split_point <- floor(n_cells / 2)
    left_genes <- seq_len(max(1, floor(n_genes / 2)))
    right_genes <- seq(from = max(left_genes) + 1, to = n_genes)

    counts[left_genes, seq_len(split_point)] <- counts[left_genes, seq_len(split_point)] + 12
    counts[right_genes, seq.int(split_point + 1, n_cells)] <- counts[right_genes, seq.int(split_point + 1, n_cells)] + 12
  }

  rownames(counts) <- c(
    paste0("MT-", seq_len(mt_genes)),
    paste0("Gene", seq_len(n_genes - mt_genes))
  )
  colnames(counts) <- paste0("Cell", seq_len(n_cells))

  Matrix::Matrix(counts, sparse = TRUE)
}

make_toy_seurat <- function(counts = make_toy_counts(),
                           sample_labels = NULL,
                           clusters = NULL,
                           cell_annotations = NULL) {
  obj <- Seurat::CreateSeuratObject(
    counts = counts,
    project = "CellivistaTest",
    min.cells = 0,
    min.features = 0
  )

  if (is.null(sample_labels)) {
    sample_labels <- rep(c("Sample1", "Sample2"), length.out = ncol(obj))
  }
  obj$SampleLabel <- sample_labels

  if (!is.null(clusters)) {
    obj$seurat_clusters <- factor(clusters)
  }

  if (!is.null(cell_annotations)) {
    obj$cell_annotation <- factor(cell_annotations)
  }

  obj
}

make_umap_ready_seurat <- function(counts = make_toy_counts(n_genes = 10, n_cells = 8)) {
  obj <- make_toy_seurat(
    counts = counts,
    sample_labels = rep(c("Sample1", "Sample2"), each = ncol(counts) / 2),
    clusters = rep(c("0", "1"), each = ncol(counts) / 2),
    cell_annotations = rep(c("T cell", "B cell"), each = ncol(counts) / 2)
  )

  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = min(6, nrow(obj)), verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, npcs = min(4, nrow(obj) - 1, ncol(obj) - 1), verbose = FALSE)
  obj <- suppressWarnings(Seurat::RunUMAP(obj, dims = 1:2, n.neighbors = 3, verbose = FALSE))
  obj
}

make_marker_ready_seurat <- function() {
  counts <- matrix(2L, nrow = 12, ncol = 40)
  counts[1:6, 1:20] <- 60L
  counts[1:6, 21:40] <- 2L
  counts[7:12, 1:20] <- 2L
  counts[7:12, 21:40] <- 60L
  counts <- counts + matrix(
    rpois(length(counts), lambda = 2),
    nrow = nrow(counts),
    ncol = ncol(counts)
  )

  rownames(counts) <- paste0("Gene", seq_len(12))
  colnames(counts) <- paste0("Cell", seq_len(40))

  obj <- Seurat::CreateSeuratObject(
    counts = Matrix::Matrix(counts, sparse = TRUE),
    project = "MarkerTest",
    min.cells = 0,
    min.features = 0
  )
  obj$SampleLabel <- rep(c("Sample1", "Sample2"), each = 20)
  obj$seurat_clusters <- factor(rep(c("0", "1"), each = 20))
  Seurat::Idents(obj) <- obj$seurat_clusters
  obj <- Seurat::SCTransform(obj, verbose = FALSE)
  Seurat::Idents(obj) <- obj$seurat_clusters
  obj
}

write_matrix_market_fixture <- function(matrix_path, cell_path, gene_path) {
  counts <- Matrix::Matrix(
    c(
      5, 0, 1,
      0, 4, 0
    ),
    nrow = 2,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(counts) <- c("ENSG1", "ENSG2")
  colnames(counts) <- c("cellA-1", "cellB-1", "cellC-2")

  Matrix::writeMM(counts, matrix_path)
  write.table(
    c("cellA-1", "cellB-1", "cellC-2"),
    file = cell_path,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  write.table(
    data.frame(
      gene_id = c("ENSG1", "ENSG2"),
      gene_name = c("GeneA", "GeneB")
    ),
    file = gene_path,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t"
  )
}

write_minimal_10x_h5 <- function(filepath, counts, barcodes, gene_names) {
  stopifnot(length(barcodes) == ncol(counts), length(gene_names) == nrow(counts))

  if (file.exists(filepath)) {
    file.remove(filepath)
  }

  sparse_counts <- as(counts, "dgCMatrix")

  rhdf5::h5createFile(filepath)
  rhdf5::h5createGroup(filepath, "matrix")
  rhdf5::h5createGroup(filepath, "matrix/features")

  rhdf5::h5write(as.integer(sparse_counts@x), filepath, "matrix/data")
  rhdf5::h5write(as.integer(sparse_counts@i), filepath, "matrix/indices")
  rhdf5::h5write(as.integer(sparse_counts@p), filepath, "matrix/indptr")
  rhdf5::h5write(as.integer(dim(sparse_counts)), filepath, "matrix/shape")
  rhdf5::h5write(as.character(barcodes), filepath, "matrix/barcodes")
  rhdf5::h5write(as.character(gene_names), filepath, "matrix/features/name")
  rhdf5::h5write(as.character(gene_names), filepath, "matrix/features/id")
  rhdf5::h5write(rep("Gene Expression", length(gene_names)), filepath, "matrix/features/feature_type")
  rhdf5::h5write(rep("GRCh38", length(gene_names)), filepath, "matrix/features/genome")

  invisible(filepath)
}
