#' Plot a Cluster UMAP with a Selectable Color Palette
#'
#' Draws the UMAP colored by group using Seurat's default colors, the
#' colorblind-safe Okabe-Ito palette, or viridis. Okabe-Ito has only 8
#' colors, so with more than 8 groups the plot falls back to viridis and
#' returns a notice explaining why.
#'
#' @param seurat_obj A Seurat object with a UMAP reduction.
#' @param palette One of "default", "okabe_ito", or "viridis". Default is
#'   "default".
#' @param group_by Metadata column used to color cells, or `"ident"` for the
#'   active identity. Default is NULL, which uses [default_group_column()]
#'   (the cell-type annotation when present, otherwise the clusters).
#' @param label Logical; draw group names on the plot, repelled so they do not
#'   overlap. Default is FALSE.
#' @param reduction Name of the reduction to plot. Default is "umap".
#'
#' @return A list with `plot` (a ggplot object), `palette_used` (the palette
#'   actually applied), and `notice` (a character message, or NULL).
#' @export
#'
#' @examples
#' \dontrun{
#' res <- plot_umap_palette(seurat_obj, palette = "okabe_ito", label = TRUE)
#' res$plot
#' }
plot_umap_palette <- function(seurat_obj,
                              palette = c("default", "okabe_ito", "viridis"),
                              group_by = NULL,
                              label = FALSE,
                              reduction = "umap") {
  palette <- match.arg(palette)
  if (is.null(seurat_obj) || !inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop("No '", reduction, "' reduction found. Run the Downstream Analysis tab first.")
  }

  groups <- resolve_viz_groups(seurat_obj, group_by)
  group_levels <- groups$levels
  n_groups <- length(group_levels)

  palette_used <- palette
  notice <- NULL
  okabe_ito <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000"
  )
  if (identical(palette, "okabe_ito") && n_groups > length(okabe_ito)) {
    palette_used <- "viridis"
    notice <- paste0(
      "Okabe-Ito has only 8 colors but there are ", n_groups,
      " groups, so viridis was used instead."
    )
  }

  cols <- switch(
    palette_used,
    default = NULL,
    okabe_ito = stats::setNames(okabe_ito[seq_len(n_groups)], group_levels),
    viridis = stats::setNames(viridisLite::viridis(n_groups), group_levels)
  )

  plot_obj <- Seurat::DimPlot(
    seurat_obj,
    reduction = reduction,
    group.by = groups$column,
    cols = cols,
    label = isTRUE(label),
    repel = isTRUE(label)
  )

  list(plot = plot_obj, palette_used = palette_used, notice = notice)
}
