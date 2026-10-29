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
  expect_equal(res$ranking, "avg_log2FC")
  expect_match(res$notice, "Dropped marker cluster")

  wrong_genes <- transform(markers, gene = paste0("Other", seq_len(nrow(markers))))
  expect_error(plot_marker_heatmap(obj, wrong_genes), "don't match this object")

  wrong_clusters <- transform(markers, cluster = "42")
  expect_error(plot_marker_heatmap(obj, wrong_clusters), "don't match this object")

  expect_error(plot_marker_heatmap(obj, NULL), "Run the Cluster-wise DGEA tab first")
})

test_that("top15-style tables keep file order and report fewer genes than requested", {
  obj <- make_marker_ready_seurat()
  top15 <- data.frame(cluster = rep(c("0", "1"), each = 2), gene = c("Gene3", "Gene1", "Gene9", "Gene7"))

  res <- quiet_toy(plot_marker_heatmap(obj, top15, n_genes = 5))
  expect_equal(res$ranking, "file order")
  expect_equal(res$genes, c("Gene3", "Gene1", "Gene9", "Gene7"))
  expect_match(res$notice, "Requested top 5 genes per cluster, but the marker table has only 2")
})

test_that("read_marker_table validates columns and choose_marker_source prefers uploads", {
  dir <- tempfile(pattern = "markers-")
  dir.create(dir)
  good <- file.path(dir, "top15_marker.csv")
  bad <- file.path(dir, "bad.csv")
  session <- file.path(dir, "all_marker.csv")
  utils::write.csv(data.frame(cluster = 0, gene = "Gene1"), good, row.names = FALSE)
  utils::write.csv(data.frame(cluster = 0, symbol = "Gene1"), bad, row.names = FALSE)
  utils::write.csv(data.frame(cluster = 0, gene = "Gene2", avg_log2FC = 1), session, row.names = FALSE)

  expect_identical(read_marker_table(good)$cluster, "0")
  expect_error(read_marker_table(bad), "missing column\\(s\\): gene")

  expect_equal(choose_marker_source(good, session), list(path = good, source = "uploaded"))
  expect_equal(choose_marker_source(NULL, session), list(path = session, source = "session"))
  expect_error(choose_marker_source(NULL, file.path(dir, "absent.csv")), "or upload the marker table")
})

test_that("clusters map to cell types many-to-one, and one-to-many falls back to clusters", {
  mapping <- map_clusters_to_groups(c("0", "1", "2", "2"), c("T cells", "T cells", "B cells", "B cells"))
  expect_equal(mapping$map, c("0" = "T cells", "1" = "T cells", "2" = "B cells"))
  expect_length(mapping$ambiguous, 0)

  obj <- make_marker_ready_seurat()
  markers <- data.frame(cluster = rep(c("0", "1"), each = 3), gene = paste0("Gene", c(1:3, 7:9)), avg_log2FC = 3:1)
  obj$cell_annotation <- factor(ifelse(obj$seurat_clusters == "0", "T cells", "B cells"))
  mapped <- quiet_toy(plot_marker_heatmap(obj, markers, n_genes = 3))
  expect_match(mapped$grouping, "2 cluster\\(s\\) mapped to 2 cell_annotation group\\(s\\)")
  expect_null(mapped$notice)

  obj$cell_annotation <- factor(rep(c("T cells", "NK cells"), length.out = ncol(obj)))
  fallback <- quiet_toy(plot_marker_heatmap(obj, markers, n_genes = 3))
  expect_match(fallback$notice, "more than one cell_annotation value, so the heatmap is grouped by cluster")
  expect_equal(fallback$grouping, "Grouped by Cluster.")
})
