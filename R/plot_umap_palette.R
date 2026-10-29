#' Plot a Cluster UMAP with a Selectable Color Palette
#'
#' Draws the UMAP colored by cluster using Seurat's default colors, the
#' colorblind-safe Okabe-Ito palette, or viridis. Okabe-Ito has only 8
#' colors, so with more than 8 clusters the plot falls back to viridis and
#' returns a notice explaining why.
#'
#' @param seurat_obj A Seurat object with a UMAP reduction.
#' @param palette One of "default", "okabe_ito", or "viridis". Default is
#'   "default".
#' @param cluster_col Metadata column holding cluster labels. Default is NULL,
#'   which uses [detect_cluster_column()] and falls back to the active
#'   identities.
#' @param reduction Name of the reduction to plot. Default is "umap".
#'
#' @return A list with `plot` (a ggplot object), `palette_used` (the palette
#'   actually applied), and `notice` (a character message, or NULL).
#' @export
#'
#' @examples
#' \dontrun{
#' res <- plot_umap_palette(seurat_obj, palette = "okabe_ito")
#' res$plot
#' }
plot_umap_palette <- function(seurat_obj,
                              palette = c("default", "okabe_ito", "viridis"),
                              cluster_col = NULL,
                              reduction = "umap") {
  palette <- match.arg(palette)
  if (is.null(seurat_obj) || !inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop("No '", reduction, "' reduction found. Run the Downstream Analysis tab first.")
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
      " clusters, so viridis was used instead."
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
    group.by = cluster_col,
    cols = cols
  )

  list(plot = plot_obj, palette_used = palette_used, notice = notice)
}
