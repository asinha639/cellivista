#' Plot Cluster Composition as Stacked Percent Bars
#'
#' Draws stacked bars showing the percent of cells in each group. With
#' `bars = "samples"` (default) there is one bar per sample, filled by
#' cluster, so each bar shows that sample's cluster make-up and sums to 100%.
#' With `bars = "clusters"` there is one bar per cluster, filled by sample.
#'
#' If the object has no usable sample column, a single "All cells" bar filled
#' by cluster is drawn and a notice explains why.
#'
#' @param seurat_obj A clustered Seurat object.
#' @param cluster_col Metadata column holding cluster labels. Default is NULL,
#'   which uses [detect_cluster_column()] and falls back to the active
#'   identities.
#' @param sample_col Metadata column holding sample labels. Default is NULL,
#'   which uses [detect_sample_column()].
#' @param bars Either "samples" (one bar per sample) or "clusters" (one bar
#'   per cluster). Default is "samples".
#'
#' @return A list with `plot` (a ggplot object), `data` (a data frame with
#'   columns `bar`, `fill`, `n_cells`, and `percent`), and `notice` (a
#'   character message, or NULL).
#' @export
#'
#' @examples
#' \dontrun{
#' res <- plot_cluster_composition(seurat_obj, bars = "samples")
#' res$plot
#' }
plot_cluster_composition <- function(seurat_obj,
                                     cluster_col = NULL,
                                     sample_col = NULL,
                                     bars = c("samples", "clusters")) {
  bars <- match.arg(bars)
  if (is.null(seurat_obj) || !inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  meta <- seurat_obj@meta.data

  if (is.null(cluster_col)) {
    cluster_col <- detect_cluster_column(seurat_obj)
  } else if (!cluster_col %in% colnames(meta)) {
    stop("Cluster column not found in metadata: ", cluster_col)
  }
  clusters <- if (is.null(cluster_col)) Seurat::Idents(seurat_obj) else meta[[cluster_col]]
  cluster_levels <- if (is.factor(clusters)) {
    levels(droplevels(clusters))
  } else {
    sort(unique(as.character(clusters)))
  }
  clusters <- factor(as.character(clusters), levels = cluster_levels)

  if (is.null(sample_col)) {
    sample_col <- detect_sample_column(seurat_obj)
  } else if (!sample_col %in% colnames(meta)) {
    stop("Sample column not found in metadata: ", sample_col)
  }

  notice <- NULL
  if (is.null(sample_col)) {
    samples <- factor(rep("All cells", ncol(seurat_obj)))
    bars <- "samples"
    notice <- "No sample column found, so all cells are shown as a single bar."
  } else {
    samples <- as.character(meta[[sample_col]])
    samples[is.na(samples) | !nzchar(trimws(samples))] <- "Unlabeled"
    samples <- factor(samples)
  }

  if (identical(bars, "samples")) {
    bar <- samples
    fill <- clusters
    x_lab <- "Sample"
    fill_lab <- "Cluster"
  } else {
    bar <- clusters
    fill <- samples
    x_lab <- "Cluster"
    fill_lab <- "Sample"
  }

  counts <- as.data.frame(table(bar = bar, fill = fill), responseName = "n_cells")
  totals <- tapply(counts$n_cells, counts$bar, sum)
  counts$percent <- 100 * counts$n_cells / as.numeric(totals[as.character(counts$bar)])
  counts <- counts[counts$n_cells > 0, , drop = FALSE]

  plot_obj <- ggplot2::ggplot(
    counts,
    ggplot2::aes(x = .data$bar, y = .data$percent, fill = .data$fill)
  ) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
    ggplot2::labs(x = x_lab, y = "Percent of cells", fill = fill_lab) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  list(plot = plot_obj, data = counts, notice = notice)
}
