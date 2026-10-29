#' Read a Cluster-wise DGEA Marker Table
#'
#' Reads a marker CSV written by the Cluster-wise DGEA tab, either
#' `all_marker.csv` (with `avg_log2FC`) or `top15_marker.csv` (only
#' `cluster` and `gene`, already ranked by `avg_log2FC`).
#'
#' @param path Path to the CSV file.
#'
#' @return A data frame with character `cluster` and `gene` columns, plus any
#'   other columns in the file.
#' @export
#'
#' @examples
#' \dontrun{
#' markers <- read_marker_table("DGEA/top15_marker.csv")
#' }
read_marker_table <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop("Marker table file not found.")
  }
  markers <- tryCatch(
    as.data.frame(readr::read_csv(path, show_col_types = FALSE)),
    error = function(e) stop("Could not read the marker table as CSV: ", conditionMessage(e))
  )
  missing_cols <- setdiff(c("cluster", "gene"), colnames(markers))
  if (length(missing_cols) > 0) {
    stop(
      "Marker table is missing column(s): ", paste(missing_cols, collapse = ", "),
      ". Use all_marker.csv or top15_marker.csv from the Cluster-wise DGEA tab."
    )
  }
  if (nrow(markers) == 0) {
    stop("Marker table has no rows.")
  }
  markers$cluster <- as.character(markers$cluster)
  markers$gene <- as.character(markers$gene)
  markers
}

#' Choose the Marker Table Source for the Heatmap
#'
#' An uploaded marker table takes priority over the session's Cluster-wise
#' DGEA results.
#'
#' @param uploaded Path to an uploaded marker CSV, or NULL.
#' @param session Path to the session's `all_marker.csv`, or NULL.
#'
#' @return A list with `path` and `source` (`"uploaded"` or `"session"`).
#' @export
#'
#' @examples
#' \dontrun{
#' choose_marker_source(uploaded = NULL, session = "DGEA/all_marker.csv")
#' }
choose_marker_source <- function(uploaded = NULL, session = NULL) {
  usable <- function(path) !is.null(path) && length(path) == 1 && nzchar(path) && file.exists(path)
  if (usable(uploaded)) {
    return(list(path = uploaded, source = "uploaded"))
  }
  if (usable(session)) {
    return(list(path = session, source = "session"))
  }
  stop("No marker results found. Run the Cluster-wise DGEA tab first or upload the marker table.")
}

#' Map Clusters to Another Cell Grouping
#'
#' For each cluster, collects the values of another grouping (for example the
#' cell-type annotation) among that cluster's cells. Several clusters may
#' share one value; a cluster whose cells carry more than one value is
#' reported as ambiguous. Missing values count as a separate value.
#'
#' @param clusters Per-cell cluster labels.
#' @param groups Per-cell values of the other grouping, in the same order.
#'
#' @return A list with `map` (a named character vector, cluster to group, for
#'   unambiguous clusters) and `ambiguous` (a named list, cluster to the
#'   values found, for clusters with more than one value).
#' @export
#'
#' @examples
#' map_clusters_to_groups(c("0", "0", "1"), c("T cells", "T cells", "B cells"))
map_clusters_to_groups <- function(clusters, groups) {
  if (length(clusters) != length(groups)) {
    stop("clusters and groups must have the same length.")
  }
  groups <- as.character(groups)
  groups[is.na(groups)] <- "(missing)"
  values_by_cluster <- lapply(split(groups, as.character(clusters)), unique)
  n_values <- lengths(values_by_cluster)
  list(
    map = vapply(values_by_cluster[n_values == 1], `[`, character(1), 1),
    ambiguous = values_by_cluster[n_values > 1]
  )
}

#' Plot a Heatmap of Top Cluster Markers
#'
#' Draws a `Seurat::DoHeatmap()` of the top marker genes per cluster, taken
#' from existing cluster-wise DGEA results (see [read_marker_table()]).
#' Markers are never recomputed here. Cells are downsampled per cluster
#' before scaling so large datasets stay within memory.
#'
#' Selection rule: the top `n_genes` positive markers per cluster, ranked by
#' `avg_log2FC`. When the table has no `avg_log2FC` column
#' (`top15_marker.csv`), its row order is kept, since the DGEA tab already
#' ranked it by `avg_log2FC`.
#'
#' Markers are keyed by cluster. When `group_by` is another grouping such as
#' the cell-type annotation, each cluster is mapped to its value in that
#' grouping (several clusters may share one). If any cluster maps to more than
#' one value, the heatmap is grouped by cluster instead and a notice says so.
#'
#' Marker clusters or genes that are absent from `seurat_obj` are dropped with
#' a notice. If no clusters match, or fewer than half of the selected marker
#' genes are present in the object, the function stops instead of drawing a
#' misleading heatmap.
#'
#' @param seurat_obj A clustered Seurat object.
#' @param markers A data frame of marker results with columns `cluster` and
#'   `gene`, and optionally `avg_log2FC`.
#' @param n_genes Number of top positive markers per cluster. Default is 5.
#' @param cells_per_cluster Maximum number of cells drawn per cluster.
#'   Default is 200.
#' @param group_by Grouping shown on the heatmap: a metadata column, or
#'   `"ident"` for the active identity. Default is NULL, which uses
#'   [default_group_column()].
#' @param cluster_col Metadata column the marker clusters refer to. Default is
#'   NULL, which uses [detect_cluster_column()] and falls back to the active
#'   identities.
#' @param assay Assay to plot. Default is NULL, which uses "SCT" when present
#'   (matching the DGEA tab) and otherwise the default assay.
#' @param seed Random seed for downsampling. The global RNG state is restored
#'   afterwards. Default is 1.
#'
#' @return A list with `plot` (a ggplot object), `genes` (the plotted genes),
#'   `n_cells` (cells drawn), `ranking` (`"avg_log2FC"` or `"file order"`),
#'   `grouping` (a description of how cells are grouped), and `notice` (a
#'   character message, or NULL).
#' @export
#'
#' @examples
#' \dontrun{
#' markers <- read_marker_table("DGEA/all_marker.csv")
#' res <- plot_marker_heatmap(seurat_obj, markers, n_genes = 10)
#' res$plot
#' }
plot_marker_heatmap <- function(seurat_obj,
                                markers,
                                n_genes = 5,
                                cells_per_cluster = 200,
                                group_by = NULL,
                                cluster_col = NULL,
                                assay = NULL,
                                seed = 1) {
  if (is.null(seurat_obj) || !inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  if (is.null(markers) || !is.data.frame(markers) || nrow(markers) == 0) {
    stop("No marker results available. Run the Cluster-wise DGEA tab first or upload the marker table.")
  }
  missing_cols <- setdiff(c("cluster", "gene"), colnames(markers))
  if (length(missing_cols) > 0) {
    stop("Marker table is missing column(s): ", paste(missing_cols, collapse = ", "))
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
  key_name <- if (is.null(cluster_col)) "ident" else cluster_col
  cluster_groups <- resolve_viz_groups(seurat_obj, key_name)
  groups <- cluster_groups$values
  group_levels <- cluster_groups$levels

  if (is.null(assay)) {
    assay <- if ("SCT" %in% SeuratObject::Assays(seurat_obj)) "SCT" else SeuratObject::DefaultAssay(seurat_obj)
  } else if (!assay %in% SeuratObject::Assays(seurat_obj)) {
    stop("Assay not found in seurat_obj: ", assay)
  }

  # Top positive markers per cluster, ranked as in top15_marker.csv.
  markers <- as.data.frame(markers)
  markers$cluster <- as.character(markers$cluster)
  markers$gene <- as.character(markers$gene)
  if ("avg_log2FC" %in% colnames(markers)) {
    ranking <- "avg_log2FC"
    markers <- markers[!is.na(markers$avg_log2FC) & markers$avg_log2FC > 0, , drop = FALSE]
    markers <- markers[order(markers$cluster, -markers$avg_log2FC), , drop = FALSE]
  } else {
    ranking <- "file order"
    markers <- markers[order(markers$cluster, seq_len(nrow(markers))), , drop = FALSE]
  }
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

  per_cluster <- table(factor(top$cluster, levels = matched_clusters))
  short <- per_cluster[per_cluster < n_genes]
  if (length(short) > 0) {
    used <- if (min(short) == max(short)) min(short) else paste0(min(short), "-", max(short))
    notices <- c(notices, paste0(
      "Requested top ", n_genes, " genes per cluster, but the marker table has only ", used,
      " for ", length(short), " cluster(s); those were used."
    ))
  }

  requested_genes <- unique(top$gene[top$cluster %in% matched_clusters])
  present_genes <- requested_genes[requested_genes %in% rownames(seurat_obj[[assay]])]
  if (length(present_genes) == 0 || length(present_genes) < length(requested_genes) / 2) {
    stop(mismatch_msg)
  }
  missing_genes <- setdiff(requested_genes, present_genes)
  if (length(missing_genes) > 0) {
    notices <- c(notices, paste0("Dropped ", length(missing_genes), " marker gene(s) not in this object."))
  }

  # Map clusters to the requested grouping; fall back to clusters if ambiguous.
  display_col <- if (is.null(group_by)) default_group_column(seurat_obj) else group_by
  display <- groups
  grouping <- paste0("Grouped by ", cluster_groups$label, ".")
  if (!identical(display_col, key_name)) {
    target <- resolve_viz_groups(seurat_obj, display_col)
    mapping <- map_clusters_to_groups(groups, target$values)
    if (length(mapping$ambiguous) > 0) {
      details <- vapply(names(mapping$ambiguous), function(cl) {
        paste0(cl, " (", paste(mapping$ambiguous[[cl]], collapse = ", "), ")")
      }, character(1))
      notices <- c(notices, paste0(
        "Cluster(s) ", paste(details, collapse = "; "), " contain more than one ", target$label,
        " value, so the heatmap is grouped by cluster."
      ))
    } else {
      mapped <- unname(mapping$map[as.character(groups)])
      display <- factor(mapped, levels = intersect(c(target$levels, "(missing)"), unique(mapped)))
      names(display) <- names(groups)
      used_map <- mapping$map[matched_clusters]
      grouping <- paste0(
        "Grouped by ", target$label, ": ", length(used_map), " cluster(s) mapped to ",
        length(unique(used_map)), " ", target$label, " group(s)."
      )
    }
  }

  top <- top[top$cluster %in% matched_clusters & top$gene %in% present_genes, , drop = FALSE]
  cluster_display <- vapply(top$cluster, function(cl) as.character(display[groups == cl][1]), character(1))
  top <- top[order(match(cluster_display, levels(display)), match(top$cluster, group_levels)), , drop = FALSE]
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
    cells <- names(groups)[!is.na(groups) & groups == lvl]
    if (length(cells) > cells_per_cluster) sample(cells, cells_per_cluster) else cells
  }), use.names = FALSE)

  sub_obj <- subset(seurat_obj, cells = keep_cells)
  SeuratObject::DefaultAssay(sub_obj) <- assay
  Seurat::Idents(sub_obj) <- droplevels(display[colnames(sub_obj)])
  sub_obj <- Seurat::ScaleData(sub_obj, features = genes, assay = assay, verbose = FALSE)

  plot_obj <- Seurat::DoHeatmap(sub_obj, features = genes, assay = assay)

  list(
    plot = plot_obj,
    genes = genes,
    n_cells = length(keep_cells),
    ranking = ranking,
    grouping = grouping,
    notice = if (length(notices) > 0) paste(notices, collapse = " ") else NULL
  )
}
