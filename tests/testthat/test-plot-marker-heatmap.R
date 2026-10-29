test_that("plot_marker_heatmap downsamples cells and rejects mismatched markers", {
  obj <- make_marker_ready_seurat()
  markers <- data.frame(
    cluster = rep(c("0", "1", "9"), each = 6),
    gene = paste0("Gene", c(1:6, 7:12, 1:6)),
    avg_log2FC = rep(6:1, 3)
  )

  res <- quiet_toy(plot_marker_heatmap(obj, markers, n_genes = 3, cells_per_cluster = 5))
  expect_s3_class(res$plot, "ggplot")
  expect_equal(res$n_cells, 10)
  expect_equal(res$genes, c("Gene1", "Gene2", "Gene3", "Gene7", "Gene8", "Gene9"))
  expect_match(res$notice, "Dropped marker cluster")

  wrong_genes <- transform(markers, gene = paste0("Other", seq_len(nrow(markers))))
  expect_error(plot_marker_heatmap(obj, wrong_genes), "don't match this object")

  wrong_clusters <- transform(markers, cluster = "42")
  expect_error(plot_marker_heatmap(obj, wrong_clusters), "don't match this object")

  expect_error(plot_marker_heatmap(obj, NULL), "Run the Cluster-wise DGEA tab first")
})
