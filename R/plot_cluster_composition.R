#' Plot Cluster Composition as Stacked Percent Bars
#'
#' Draws stacked bars showing the percent of cells in each group. With
#' `bars = "samples"` (default) there is one bar per sample, filled by group,
#' so each bar shows that sample's make-up and sums to 100%. With
#' `bars = "groups"` there is one bar per group, filled by sample.
#'
#' If the object has no usable sample column, a single "All cells" bar filled
#' by group is drawn and a notice explains why.
#'
#' @param seurat_obj A clustered Seurat object.
#' @param group_by Metadata column used to group cells, or `"ident"` for the
#'   active identity. Default is NULL, which uses [default_group_column()]
#'   (the cell-type annotation when present, otherwise the clusters).
#' @param sample_col Metadata column holding sample labels. Default is NULL,
#'   which uses [detect_sample_column()].
#' @param bars Either "samples" (one bar per sample) or "groups" (one bar per
#'   group). Default is "samples".
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
                                     group_by = NULL,
                                     sample_col = NULL,
                                     bars = c("samples", "groups")) {
  bars <- match.arg(bars)
  if (is.null(seurat_obj) || !inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  meta <- seurat_obj@meta.data
  groups <- resolve_viz_groups(seurat_obj, group_by)

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
    fill <- groups$values
    x_lab <- "Sample"
    fill_lab <- groups$label
  } else {
    bar <- groups$values
    fill <- samples
    x_lab <- groups$label
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
