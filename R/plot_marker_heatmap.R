#' Plot a Heatmap of Top Cluster Markers
#'
#' Draws a `Seurat::DoHeatmap()` of the top marker genes per cluster, taken
#' from existing cluster-wise DGEA results (the `all_marker.csv` written by
#' [run_clusterwise_dgea()]). Markers are never recomputed here. Cells are
#' downsampled per cluster before scaling so large datasets stay within memory.
#'
#' Marker clusters or genes that are absent from `seurat_obj` are dropped with
#' a notice. If no clusters match, or fewer than half of the selected marker
#' genes are present in the object, the function stops instead of drawing a
#' misleading heatmap.
#'
#' @param seurat_obj A clustered Seurat object.
#' @param markers A data frame of marker results with columns `cluster`,
#'   `gene`, and `avg_log2FC`.
#' @param n_genes Number of top positive markers per cluster. Default is 5.
#' @param cells_per_cluster Maximum number of cells drawn per cluster.
#'   Default is 200.
#' @param cluster_col Metadata column holding cluster labels. Default is NULL,
#'   which uses [detect_cluster_column()] and falls back to the active
#'   identities.
#' @param assay Assay to plot. Default is NULL, which uses "SCT" when present
#'   (matching the DGEA tab) and otherwise the default assay.
#' @param seed Random seed for downsampling. The global RNG state is restored
#'   afterwards. Default is 1.
#'
#' @return A list with `plot` (a ggplot object), `genes` (the plotted genes),
#'   `n_cells` (cells drawn), and `notice` (a character message, or NULL).
#' @export
#'
#' @examples
#' \dontrun{
#' markers <- readr::read_csv("DGEA/all_marker.csv")
#' res <- plot_marker_heatmap(seurat_obj, markers, n_genes = 10)
#' res$plot
#' }
plot_marker_heatmap <- function(seurat_obj,
                                markers,
                                n_genes = 5,
                                cells_per_cluster = 200,
                                cluster_col = NULL,
                                assay = NULL,
                                seed = 1) {
  if (is.null(seurat_obj) || !inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  if (is.null(markers) || !is.data.frame(markers) || nrow(markers) == 0) {
    stop("No marker results available. Run the Cluster-wise DGEA tab first.")
  }
  required_cols <- c("cluster", "gene", "avg_log2FC")
  if (!all(required_cols %in% colnames(markers))) {
    stop("markers must contain columns: ", paste(required_cols, collapse = ", "))
  }
  n_genes <- as.integer(n_genes)
  cells_per_cluster <- as.integer(cells_per_cluster)
  if (length(n_genes) != 1 || is.na(n_genes) || n_genes < 1) {
    stop("n_genes must be a positive integer.")
  }
  if (length(cells_per_cluster) != 1 || is.na(cells_per_cluster) || cells_per_cluster < 1) {
    stop("cells_per_cluster must be a positive integer.")
  }

  if (is.null(cluster_col)) {
    cluster_col <- detect_cluster_column(seurat_obj)
  } else if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop("Cluster column not found in metadata: ", cluster_col)
  }
  groups <- if (is.null(cluster_col)) {
    Seurat::Idents(seurat_obj)
  } else {
    seurat_obj@meta.data[[cluster_col]]
  }
  group_levels <- if (is.factor(groups)) {
    levels(droplevels(groups))
  } else {
    sort(unique(as.character(groups)))
  }
  groups <- factor(as.character(groups), levels = group_levels)
  names(groups) <- colnames(seurat_obj)

  if (is.null(assay)) {
    assay <- if ("SCT" %in% SeuratObject::Assays(seurat_obj)) "SCT" else SeuratObject::DefaultAssay(seurat_obj)
  } else if (!assay %in% SeuratObject::Assays(seurat_obj)) {
    stop("Assay not found in seurat_obj: ", assay)
  }

  # Top positive markers per cluster, ranked as in top15_marker.csv.
  markers <- as.data.frame(markers)
  markers$cluster <- as.character(markers$cluster)
  markers <- markers[!is.na(markers$avg_log2FC) & markers$avg_log2FC > 0, , drop = FALSE]
  markers <- markers[order(markers$cluster, -markers$avg_log2FC), , drop = FALSE]
  top <- do.call(rbind, lapply(split(markers, markers$cluster), utils::head, n = n_genes))
  if (is.null(top) || nrow(top) == 0) {
    stop("No positive marker genes found in the marker results.")
  }

  notices <- character(0)
  marker_clusters <- unique(top$cluster)
  matched_clusters <- intersect(group_levels, marker_clusters)
  mismatch_msg <- "These markers don't match this object, re-run the DGEA tab."
  if (length(matched_clusters) == 0) {
    stop(mismatch_msg)
  }
  dropped_clusters <- setdiff(marker_clusters, matched_clusters)
  if (length(dropped_clusters) > 0) {
    notices <- c(notices, paste0("Dropped marker cluster(s) not in this object: ", paste(dropped_clusters, collapse = ", "), "."))
  }

  requested_genes <- unique(top$gene)
  present_genes <- requested_genes[requested_genes %in% rownames(seurat_obj[[assay]])]
  if (length(present_genes) == 0 || length(present_genes) < length(requested_genes) / 2) {
    stop(mismatch_msg)
  }
  missing_genes <- setdiff(requested_genes, present_genes)
  if (length(missing_genes) > 0) {
    notices <- c(notices, paste0("Dropped ", length(missing_genes), " marker gene(s) not in this object."))
  }

  top <- top[top$cluster %in% matched_clusters & top$gene %in% present_genes, , drop = FALSE]
  top <- top[order(match(top$cluster, group_levels)), , drop = FALSE]
  genes <- unique(top$gene)

  # Downsample cells per cluster without disturbing the caller's RNG state.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  keep_cells <- unlist(lapply(group_levels, function(lvl) {
    cells <- names(groups)[groups == lvl]
    if (length(cells) > cells_per_cluster) sample(cells, cells_per_cluster) else cells
  }), use.names = FALSE)

  sub_obj <- subset(seurat_obj, cells = keep_cells)
  SeuratObject::DefaultAssay(sub_obj) <- assay
  Seurat::Idents(sub_obj) <- groups[colnames(sub_obj)]
  sub_obj <- Seurat::ScaleData(sub_obj, features = genes, assay = assay, verbose = FALSE)

  plot_obj <- Seurat::DoHeatmap(sub_obj, features = genes, assay = assay)

  list(
    plot = plot_obj,
    genes = genes,
    n_cells = length(keep_cells),
    notice = if (length(notices) > 0) paste(notices, collapse = " ") else NULL
  )
}
