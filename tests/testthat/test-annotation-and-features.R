test_that("annotate_clusters maps annotations from a vector and a CSV file", {
  obj <- make_umap_ready_seurat()

  vector_result <- annotate_clusters(
    seurat_obj = obj,
    annotations = c("T cell", "B cell"),
    output_dir = tempfile(pattern = "annot-"),
    save_rds = FALSE,
    save_plots = FALSE
  )

  expect_s4_class(vector_result, "Seurat")
  expect_equal(as.character(vector_result$cell_annotation), c("T cell", "T cell", "T cell", "T cell", "B cell", "B cell", "B cell", "B cell"))
  expect_equal(levels(Seurat::Idents(vector_result)), c("T cell", "B cell"))

  annotation_csv <- tempfile(fileext = ".csv")
  write.csv(
    data.frame(
      seurat_clusters = c("0", "1", "2"),
      cell_annotation = c("Myeloid", "Lymphoid", "Unused"),
      stringsAsFactors = FALSE
    ),
    file = annotation_csv,
    row.names = FALSE,
    quote = FALSE
  )

  expect_warning(
    csv_result <- annotate_clusters(
      seurat_obj = obj,
      annotations = annotation_csv,
      output_dir = tempfile(pattern = "annot-csv-"),
      save_rds = FALSE,
      save_plots = FALSE
    ),
    "unused cluster"
  )
  expect_s4_class(csv_result, "Seurat")
  expect_equal(as.character(csv_result$cell_annotation)[1], "Myeloid")

  expect_error(
    annotate_clusters(
      seurat_obj = obj,
      annotations = c("OnlyOne"),
      output_dir = tempfile(pattern = "annot-bad-"),
      save_rds = FALSE,
      save_plots = FALSE
    ),
    "annotations must be either a valid CSV path or a character vector matching the number of clusters"
  )
})

test_that("plot_genes_feature returns ggplot objects and handles missing genes", {
  obj <- make_umap_ready_seurat()

  plots <- plot_genes_feature(
    seurat_obj = obj,
    genes = "Gene1, Gene2, MissingGene, Gene1",
    cell_type = "T cell",
    split_by = "SampleLabel",
    output_dir = tempfile(pattern = "feature-"),
    save_plots = FALSE
  )

  expect_type(plots, "list")
  expect_equal(names(plots), c("Gene1", "Gene2"))
  expect_true(all(vapply(plots, inherits, logical(1), what = "ggplot")))

  expect_error(
    plot_genes_feature(
      seurat_obj = obj,
      genes = character(0),
      output_dir = tempfile(pattern = "feature-empty-"),
      save_plots = FALSE
    ),
    "No valid genes were provided"
  )

  expect_error(
    plot_genes_feature(
      seurat_obj = obj,
      genes = "Gene1",
      split_by = "MissingColumn",
      output_dir = tempfile(pattern = "feature-bad-"),
      save_plots = FALSE
    ),
    "split_by not found"
  )
})
