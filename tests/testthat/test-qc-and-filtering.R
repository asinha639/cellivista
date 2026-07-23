test_that("run_qc_metrics adds percent.mt and respects validation checks", {
  obj <- make_toy_seurat(
    counts = make_toy_counts(n_genes = 6, n_cells = 4, mt_genes = 2),
    sample_labels = rep(c("Sample1", "Sample2"), each = 2)
  )

  out <- run_qc_metrics(
    seurat_obj = obj,
    mt_pattern = "^MT-",
    output_dir = tempfile(pattern = "qc-"),
    sample_split_var = "SampleLabel",
    ylimit = c(0, 100),
    save_plots = FALSE
  )

  expect_s4_class(out, "Seurat")
  expect_true("percent.mt" %in% colnames(out@meta.data))
  expect_equal(length(out$percent.mt), ncol(obj))
  expect_true(all(out$percent.mt >= 0))

  expect_error(
    run_qc_metrics(obj, mt_pattern = "", save_plots = FALSE),
    "mt_pattern must be a non-empty character string"
  )
  expect_error(
    run_qc_metrics(obj, sample_split_var = "MissingColumn", save_plots = FALSE),
    "sample_split_var not found"
  )
})

test_that("post_qc_filtering_and_plots filters cells and saves summary plots", {
  obj <- make_toy_seurat(
    counts = make_toy_counts(n_genes = 6, n_cells = 4, mt_genes = 2),
    sample_labels = rep(c("Sample1", "Sample2"), each = 2)
  )
  obj <- run_qc_metrics(obj, save_plots = FALSE)

  out_dir <- tempfile(pattern = "post-qc-")
  rds_path <- tempfile(fileext = ".rds")

  filtered <- post_qc_filtering_and_plots(
    seurat_obj = obj,
    min_features = 0,
    max_features = 10,
    max_mt_percent = 90,
    output_dir = out_dir,
    sample_split_var = "SampleLabel",
    save_rds = TRUE,
    rds_path = rds_path
  )

  expect_s4_class(filtered, "Seurat")
  expect_lte(ncol(filtered), ncol(obj))
  expect_true(file.exists(file.path(out_dir, "post_qc_vplot.jpg")))
  expect_true(file.exists(file.path(out_dir, "post_qc_densityplot.jpg")))
  expect_true(file.exists(file.path(out_dir, "cellcount_piechart_post.jpg")))
  expect_true(file.exists(rds_path))

  expect_error(
    post_qc_filtering_and_plots(obj, min_features = 10, max_features = 5, save_rds = FALSE),
    "min_features must be less than max_features"
  )
})
