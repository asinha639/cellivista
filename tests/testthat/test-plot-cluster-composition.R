test_that("plot_cluster_composition draws 100% bars per sample or cluster", {
  obj <- make_marker_ready_seurat()
  obj$SampleLabel <- rep(c("A", "B"), length.out = ncol(obj))

  res <- plot_cluster_composition(obj)
  expect_s3_class(res$plot, "ggplot")
  expect_null(res$notice)
  expect_setequal(as.character(unique(res$data$bar)), c("A", "B"))
  expect_equal(as.numeric(tapply(res$data$percent, res$data$bar, sum)), c(100, 100))

  by_cluster <- plot_cluster_composition(obj, bars = "clusters")
  expect_setequal(as.character(unique(by_cluster$data$fill)), c("A", "B"))

  obj$SampleLabel <- NULL
  single <- plot_cluster_composition(obj, bars = "clusters")
  expect_equal(as.character(unique(single$data$bar)), "All cells")
  expect_equal(sum(single$data$percent), 100)
  expect_match(single$notice, "No sample column")
})
