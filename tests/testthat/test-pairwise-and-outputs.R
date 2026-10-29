test_that("run_pairwise_dgea compares two groups and writes a CSV", {
  obj <- make_marker_ready_seurat()
  obj <- Seurat::NormalizeData(obj, assay = "RNA", verbose = FALSE)

  out_path <- tempfile(fileext = ".csv")
  res <- run_pairwise_dgea(
    seurat_obj = obj,
    group_by = "seurat_clusters",
    ident1 = "0",
    ident2 = "1",
    output_path = out_path,
    logfc_threshold = 0,
    min_pct = 0
  )

  expect_true(is.data.frame(res))
  expect_true(all(c("avg_log2FC", "p_val", "gene") %in% colnames(res)))
  expect_true(file.exists(out_path))
  expect_equal(nrow(readr::read_csv(out_path, show_col_types = FALSE)), nrow(res))

  # Genes 1-6 are high in cluster 0, so they should be up in ident1.
  up <- res$gene[res$avg_log2FC > 0]
  expect_true(all(paste0("Gene", 1:6) %in% up))
})

test_that("run_pairwise_dgea errors on an unknown grouping column", {
  obj <- make_marker_ready_seurat()
  expect_error(
    run_pairwise_dgea(
      obj,
      group_by = "NotAColumn",
      ident1 = "a",
      ident2 = "b",
      output_path = tempfile(fileext = ".csv")
    )
  )
})

test_that("plot_genes_feature saves one file per valid gene", {
  obj <- make_umap_ready_seurat()
  out_dir <- tempfile(pattern = "feature-save-")

  plots <- plot_genes_feature(
    seurat_obj = obj,
    genes = c("Gene1", "Gene2"),
    output_dir = out_dir,
    save_plots = TRUE
  )

  expect_equal(names(plots), c("Gene1", "Gene2"))
  saved <- list.files(out_dir, recursive = TRUE)
  expect_length(saved, 2)
})

test_that("plot_genes_feature validates cell types and genes", {
  obj <- make_umap_ready_seurat()

  expect_error(
    plot_genes_feature(obj, genes = "Gene1", cell_type = "Neuron",
                       output_dir = tempfile(), save_plots = FALSE),
    "Cell type not found"
  )
  expect_error(
    plot_genes_feature(obj, genes = "NotAGene",
                       output_dir = tempfile(), save_plots = FALSE),
    "None of the requested genes were found"
  )
  expect_error(
    plot_genes_feature(NULL, genes = "Gene1",
                       output_dir = tempfile(), save_plots = FALSE),
    "seurat_obj is NULL"
  )

  no_annotation <- obj
  no_annotation$cell_annotation <- NULL
  expect_error(
    plot_genes_feature(no_annotation, genes = "Gene1", cell_type = "T cell",
                       output_dir = tempfile(), save_plots = FALSE),
    "cell_annotation"
  )
})

test_that("run_qc_metrics saves plots to the output directory", {
  obj <- make_toy_seurat(
    counts = make_toy_counts(n_genes = 6, n_cells = 4, mt_genes = 2),
    sample_labels = rep(c("Sample1", "Sample2"), each = 2)
  )
  out_dir <- tempfile(pattern = "qc-plots-")

  out <- suppressWarnings(run_qc_metrics(obj, output_dir = out_dir, save_plots = TRUE))

  expect_true("percent.mt" %in% colnames(out@meta.data))
  expect_gt(length(list.files(out_dir)), 0)
})

test_that("run_downstream_analysis saves the RDS and plots when requested", {
  obj <- make_toy_seurat(
    counts = make_toy_counts(n_genes = 10, n_cells = 8, mt_genes = 2),
    sample_labels = rep(c("Sample1", "Sample2"), each = 4)
  )
  obj <- suppressWarnings(run_qc_metrics(obj, save_plots = FALSE))
  out_dir <- tempfile(pattern = "downstream-save-")
  rds_path <- file.path(tempfile(pattern = "rds-"), "downstream.rds")

  out <- suppressWarnings(run_downstream_analysis(
    seurat_obj = obj,
    dims = 1:2,
    resolution = 0.2,
    output_dir = out_dir,
    save_rds = TRUE,
    rds_path = rds_path,
    save_plots = TRUE
  ))

  expect_true(file.exists(rds_path))
  expect_s4_class(readRDS(rds_path), "Seurat")
  expect_gt(length(list.files(out_dir)), 0)
  expect_equal(ncol(readRDS(rds_path)), ncol(out))
})

test_that("annotate_clusters rejects annotations that do not match the clusters", {
  obj <- make_umap_ready_seurat()

  expect_error(
    annotate_clusters(
      obj,
      annotations = c("OnlyOne", "Two", "Three"),
      output_dir = tempfile(),
      save_rds = FALSE,
      save_plots = FALSE
    ),
    "annotations must be either a valid CSV path"
  )
})

test_that("annotate_clusters saves the annotated object", {
  obj <- make_umap_ready_seurat()
  rds_path <- file.path(tempfile(pattern = "annot-"), "annotated.rds")

  out <- annotate_clusters(
    obj,
    annotations = c("Type A", "Type B"),
    output_dir = tempfile(),
    rds_path = rds_path,
    save_rds = TRUE,
    save_plots = FALSE
  )

  expect_true(file.exists(rds_path))
  expect_setequal(as.character(unique(out$cell_annotation)), c("Type A", "Type B"))
  expect_equal(as.character(Seurat::Idents(out)), as.character(out$cell_annotation))
})

test_that("run_clusterwise_dgea writes plots when plot_options is supplied", {
  obj <- make_marker_ready_seurat()
  out_dir <- tempfile(pattern = "dgea-plots-")

  markers <- suppressWarnings(run_clusterwise_dgea(
    seurat_obj = obj,
    output_dir = out_dir,
    logfc_threshold = 0,
    min_pct = 0,
    plot_options = list(
      violin = list(save = TRUE, format = "png"),
      dot = list(save = TRUE, format = "png"),
      feature = list(save = FALSE)
    )
  ))

  expect_true(is.data.frame(markers))
  expect_true(dir.exists(file.path(out_dir, "Plots")))
  expect_gt(length(list.files(file.path(out_dir, "Plots"), recursive = TRUE)), 0)
})

test_that("run_clusterwise_dgea requires an SCT assay and non-NULL input", {
  obj <- make_toy_seurat(clusters = rep(c("0", "1"), length.out = 6))

  expect_error(
    run_clusterwise_dgea(obj, output_dir = tempfile()),
    "SCT assay not found"
  )
  expect_error(
    run_clusterwise_dgea(NULL, output_dir = tempfile()),
    "seurat_obj is NULL"
  )
})

test_that("remove_doublets requires a Seurat object with SampleLabel", {
  expect_error(remove_doublets(list()), "seurat_obj must be a Seurat object")

  obj <- make_toy_seurat()
  obj$SampleLabel <- NULL
  expect_error(remove_doublets(obj), "SampleLabel")
})

test_that("integrate_seurat_samples validates its input", {
  expect_error(
    integrate_seurat_samples(NULL),
    "requires a valid Seurat object"
  )

  obj <- make_toy_seurat()
  obj$SampleLabel <- NULL
  expect_error(integrate_seurat_samples(obj), "SampleLabel")
})
