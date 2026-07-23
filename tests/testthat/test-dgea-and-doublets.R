test_that("run_clusterwise_dgea returns markers and writes CSV outputs", {
  obj <- make_marker_ready_seurat()

  out_dir <- tempfile(pattern = "dgea-")
  markers <- run_clusterwise_dgea(
    seurat_obj = obj,
    output_dir = out_dir,
    logfc_threshold = 0,
    min_pct = 0,
    plot_options = NULL
  )

  expect_true(is.data.frame(markers))
  expect_true(nrow(markers) > 0)
  expect_true(file.exists(file.path(out_dir, "all_marker.csv")))
  expect_true(file.exists(file.path(out_dir, "top15_marker.csv")))
  expect_true(all(c("avg_log2FC", "cluster", "gene") %in% colnames(markers)))

  expect_error(
    run_clusterwise_dgea(make_toy_seurat(), output_dir = tempfile(pattern = "dgea-bad-"), plot_options = NULL),
    "seurat_clusters not found"
  )
})

test_that("remove_doublets validates its inputs before processing", {
  obj <- make_toy_seurat(
    counts = make_toy_counts(n_genes = 8, n_cells = 6, mt_genes = 2),
    sample_labels = rep(c("Sample1", "Sample2"), each = 3)
  )

  expect_error(
    remove_doublets(NULL),
    "seurat_obj must be a Seurat object"
  )
  expect_error(
    remove_doublets(obj, doublet_rate = 0),
    "doublet_rate must be a single numeric value between 0 and 1"
  )
  expect_error(
    remove_doublets(obj, pcs = 1),
    "pcs must contain at least two positive integers"
  )
})
