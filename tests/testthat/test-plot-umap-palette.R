test_that("plot_umap_palette applies palettes and falls back from Okabe-Ito above 8 clusters", {
  obj <- make_umap_ready_seurat(counts = make_toy_counts(n_genes = 10, n_cells = 10))

  default_res <- plot_umap_palette(obj, palette = "default")
  expect_s3_class(default_res$plot, "ggplot")
  expect_equal(default_res$palette_used, "default")
  expect_null(default_res$notice)

  okabe_res <- plot_umap_palette(obj, palette = "okabe_ito")
  expect_equal(okabe_res$palette_used, "okabe_ito")
  expect_null(okabe_res$notice)

  obj$seurat_clusters <- factor(paste0("c", seq_len(ncol(obj))))
  fallback <- plot_umap_palette(obj, palette = "okabe_ito")
  expect_equal(fallback$palette_used, "viridis")
  expect_match(fallback$notice, "only 8 colors")

  no_umap <- make_toy_seurat()
  expect_error(plot_umap_palette(no_umap), "No 'umap' reduction found")
})
