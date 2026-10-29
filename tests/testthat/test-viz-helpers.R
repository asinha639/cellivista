test_that("detect_cluster_column prefers seurat_clusters, then *_snn_res.* columns", {
  obj <- make_toy_seurat(clusters = rep(c("0", "1"), length.out = 6))
  expect_equal(detect_cluster_column(obj), "seurat_clusters")

  obj$seurat_clusters <- NULL
  expect_null(detect_cluster_column(obj))

  obj$RNA_snn_res.0.5 <- factor(rep(c("0", "1"), length.out = 6))
  expect_equal(detect_cluster_column(obj), "RNA_snn_res.0.5")
})

test_that("detect_sample_column returns SampleLabel only when it has values", {
  obj <- make_toy_seurat()
  expect_equal(detect_sample_column(obj), "SampleLabel")

  obj$SampleLabel <- ""
  expect_null(detect_sample_column(obj))

  obj$SampleLabel <- NULL
  expect_null(detect_sample_column(obj))
})

test_that("group columns are categorical and the default prefers cell_annotation", {
  obj <- make_toy_seurat(clusters = rep(c("0", "1"), length.out = 6))
  obj$score <- seq_len(6)
  groups <- list_group_columns(obj)
  expect_true(all(c("seurat_clusters", "SampleLabel") %in% groups))
  expect_false("score" %in% groups)
  expect_equal(default_group_column(obj), "seurat_clusters")

  obj$cell_annotation <- factor(rep(c("T cells", "B cells"), length.out = 6))
  expect_equal(default_group_column(obj), "cell_annotation")

  obj$cell_annotation <- NULL
  obj$seurat_clusters <- NULL
  expect_equal(default_group_column(obj), "ident")
})

test_that("save_viz_plot writes png, pdf and svg with the requested extension", {
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  out_dir <- tempfile(pattern = "viz-save-")

  for (fmt in c("png", "pdf")) {
    path <- save_viz_plot(p, file.path(out_dir, "plot.jpeg"), format = fmt, width = 4, height = 3)
    expect_equal(tools::file_ext(path), fmt)
    expect_true(file.exists(path))
  }
  skip_if_not_installed("svglite")
  path <- save_viz_plot(p, file.path(out_dir, "plot"), format = "svg", width = 4, height = 3)
  expect_true(file.exists(path))

  expect_error(save_viz_plot(p, file.path(out_dir, "bad"), format = "png", width = -1),
               "width and height must be positive")
})
