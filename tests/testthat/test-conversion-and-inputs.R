test_that("convert_matrix_to_rds reads a matrix market bundle and assigns sample labels", {
  matrix_path <- tempfile(fileext = ".mtx")
  cell_path <- tempfile(fileext = ".tsv")
  gene_path <- tempfile(fileext = ".tsv")
  output_path <- tempfile(fileext = ".rds")

  write_matrix_market_fixture(matrix_path, cell_path, gene_path)

  obj <- convert_matrix_to_rds(
    matrix_path = matrix_path,
    cell_path = cell_path,
    gene_path = gene_path,
    output_path = output_path,
    min_cells = 0,
    min_features = 0
  )

  expect_s4_class(obj, "Seurat")
  expect_true(file.exists(output_path))
  expect_equal(ncol(obj), 3L)
  expect_equal(Seurat::Assays(obj), "RNA")
  expect_equal(as.character(obj$SampleLabel), c("Sample1", "Sample1", "Sample2"))
  expect_equal(rownames(obj), c("GeneA", "GeneB"))
})

test_that("convert_matrix_to_rds supports sample labels from a cell file column", {
  matrix_path <- tempfile(fileext = ".mtx")
  cell_path <- tempfile(fileext = ".tsv")
  gene_path <- tempfile(fileext = ".tsv")
  output_path <- tempfile(fileext = ".rds")

  counts <- Matrix::Matrix(
    c(
      6, 0, 1,
      0, 5, 0
    ),
    nrow = 2,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(counts) <- c("ENSG1", "ENSG2")
  colnames(counts) <- c("cellA", "cellB", "cellC")
  Matrix::writeMM(counts, matrix_path)

  write.table(
    data.frame(
      barcode = c("cellA", "cellB", "cellC"),
      sample_id = c("Alpha", "Beta", "Alpha"),
      stringsAsFactors = FALSE
    ),
    file = cell_path,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  write.table(
    data.frame(gene_id = c("ENSG1", "ENSG2")),
    file = gene_path,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t"
  )

  obj <- convert_matrix_to_rds(
    matrix_path = matrix_path,
    cell_path = cell_path,
    gene_path = gene_path,
    output_path = output_path,
    sample_label_mode = "from_cell_file_column",
    sample_label_column = 2,
    min_cells = 0,
    min_features = 0
  )

  expect_s4_class(obj, "Seurat")
  expect_equal(as.character(obj$SampleLabel), c("Alpha", "Beta", "Alpha"))
})

test_that("convert_matrix_to_rds rejects missing inputs and invalid label configuration", {
  expect_error(
    convert_matrix_to_rds("missing.mtx", "missing.tsv", "missing.tsv", tempfile(fileext = ".rds")),
    "Input file does not exist"
  )

  matrix_path <- tempfile(fileext = ".mtx")
  cell_path <- tempfile(fileext = ".tsv")
  gene_path <- tempfile(fileext = ".tsv")
  write_matrix_market_fixture(matrix_path, cell_path, gene_path)

  expect_error(
    convert_matrix_to_rds(
      matrix_path = matrix_path,
      cell_path = cell_path,
      gene_path = gene_path,
      output_path = tempfile(fileext = ".rds"),
      sample_label_mode = "from_cell_file_column",
      min_cells = 0,
      min_features = 0
    ),
    "sample_label_column is required"
  )
})

test_that("read_h5_to_seurat reads a minimal 10x H5 file and remaps labels", {
  h5_path <- tempfile(fileext = ".h5")
  counts <- Matrix::Matrix(
    c(
      5, 0, 1,
      0, 4, 0,
      2, 1, 0
    ),
    nrow = 3,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(counts) <- c("GeneA", "GeneB", "GeneC")
  colnames(counts) <- c("cell1-1", "cell2-1", "cell3-2")
  write_minimal_10x_h5(h5_path, counts, colnames(counts), rownames(counts))

  add_metadata <- data.frame(
    Treatment = c("Ctrl", "Ctrl", "Stim"),
    row.names = colnames(counts),
    stringsAsFactors = FALSE
  )

  obj <- read_h5_to_seurat(
    filepath = h5_path,
    project_name = "ReadH5Test",
    min_features = 0,
    min_cells = 0,
    add_metadata = add_metadata,
    sample_label_map = c(Sample1 = "Control", Sample2 = "Treatment")
  )

  expect_s4_class(obj, "Seurat")
  expect_equal(as.character(obj$SampleLabel), c("Control", "Control", "Treatment"))
  expect_equal(unname(obj$Treatment), c("Ctrl", "Ctrl", "Stim"))
})

test_that("read_h5_to_seurat validates inputs", {
  h5_path <- tempfile(fileext = ".h5")
  counts <- Matrix::Matrix(
    c(
      5, 0,
      0, 4
    ),
    nrow = 2,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(counts) <- c("GeneA", "GeneB")
  colnames(counts) <- c("cell1-1", "cell2-2")
  write_minimal_10x_h5(h5_path, counts, colnames(counts), rownames(counts))

  expect_error(
    read_h5_to_seurat(h5_path, sample_label_map = c("Control", "Treatment")),
    "sample_label_map must be a named character vector"
  )
})
