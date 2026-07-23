test_that("run_downstream_analysis performs the core workflow on a small object", {
  obj <- make_toy_seurat(
    counts = make_toy_counts(n_genes = 10, n_cells = 8, mt_genes = 2),
    sample_labels = rep(c("Sample1", "Sample2"), each = 4)
  )
  obj <- run_qc_metrics(obj, save_plots = FALSE)

  out <- run_downstream_analysis(
    seurat_obj = obj,
    dims = 1:2,
    resolution = 0.2,
    output_dir = tempfile(pattern = "downstream-"),
    save_rds = FALSE,
    save_plots = FALSE
  )

  expect_s4_class(out, "Seurat")
  expect_true("seurat_clusters" %in% colnames(out@meta.data))
  expect_true("umap" %in% names(out@reductions))
  expect_error(
    run_downstream_analysis(obj, dims = integer(0), save_plots = FALSE, save_rds = FALSE),
    "dims must contain at least one positive integer"
  )
})

test_that("integrate_seurat_samples integrates multiple samples and preserves metadata", {
  counts <- make_toy_counts(n_genes = 30, n_cells = 40, mt_genes = 2, structured = TRUE)
  obj <- make_toy_seurat(
    counts = counts,
    sample_labels = rep(c("Sample1", "Sample2"), each = 20)
  )

  out <- integrate_seurat_samples(
    seurat_obj = obj,
    nfeatures = 20,
    save_rds = FALSE,
    output_dir = tempfile(pattern = "integrated-"),
    plot_options = NULL
  )

  expect_s4_class(out, "Seurat")
  expect_true("integrated" %in% Seurat::Assays(out))
  expect_true("SampleLabel" %in% colnames(out@meta.data))
})
