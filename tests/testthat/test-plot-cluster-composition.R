test_that("plot_cluster_composition draws 100% bars per sample or group", {
  obj <- make_marker_ready_seurat()
  obj$SampleLabel <- rep(c("A", "B"), length.out = ncol(obj))

  res <- plot_cluster_composition(obj)
  expect_s3_class(res$plot, "ggplot")
  expect_null(res$notice)
  expect_setequal(as.character(unique(res$data$bar)), c("A", "B"))
  expect_equal(as.numeric(tapply(res$data$percent, res$data$bar, sum)), c(100, 100))

  by_group <- plot_cluster_composition(obj, bars = "groups")
  expect_setequal(as.character(unique(by_group$data$fill)), c("A", "B"))

  obj$cell_annotation <- factor(ifelse(obj$seurat_clusters == "0", "T cells", "B cells"))
  annotated <- plot_cluster_composition(obj)
  expect_setequal(as.character(unique(annotated$data$fill)), c("T cells", "B cells"))
  expect_equal(annotated$plot$labels$fill, "cell_annotation")

  obj$SampleLabel <- NULL
  single <- plot_cluster_composition(obj, bars = "groups")
  expect_equal(as.character(unique(single$data$bar)), "All cells")
  expect_equal(sum(single$data$percent), 100)
  expect_match(single$notice, "No sample column")
})
