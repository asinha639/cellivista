#' Doublet Removal Using DoubletFinder per SampleLabel
#'
#' Performs doublet detection and removal for each sample in a Seurat object,
#' using DoubletFinder.
#'
#' @param seurat_obj A Seurat object with a `SampleLabel` column in metadata.
#' @param output_dir Directory to save plots. Default is "./Plots/DoubletFinder".
#' @param save_rds Logical; whether to save the final filtered Seurat object.
#'   Default is TRUE.
#' @param rds_path Path to save filtered Seurat object.
#'   Default is "seurat.obj_doubletfiltered.rds".
#' @param doublet_rate Expected doublet rate (e.g., 0.075 for 7.5%).
#'   Default is 0.075.
#' @param save_plots Logical; whether to save DoubletFinder plots.
#'   Default is TRUE.
#' @param plot_options Optional list with `umap` and `summary` export settings.
#'
#' @return A Seurat object with doublets removed.
#' @export
remove_doublets <- function(seurat_obj,
                            output_dir = "./Plots/DoubletFinder",
                            save_rds = TRUE,
                            rds_path = "RDS_Files/seurat.obj_doubletfiltered.rds",
                            doublet_rate = 0.075,
                            save_plots = TRUE,
                            plot_options = NULL,
                            pcs = 1:15,
                            sct = FALSE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required.")
  }
  if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
    stop("DoubletFinder is required.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr is required.")
  }
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("purrr is required.")
  }

  if (is.null(seurat_obj) || !inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  meta <- seurat_obj[[]]
  if (!"SampleLabel" %in% names(meta)) {
    stop("seurat_obj must contain a 'SampleLabel' column in metadata.")
  }
  if (!is.numeric(doublet_rate) || length(doublet_rate) != 1 || is.na(doublet_rate) ||
      doublet_rate <= 0 || doublet_rate >= 1) {
    stop("doublet_rate must be a single numeric value between 0 and 1.")
  }

  pcs <- as.integer(pcs)
  pcs <- pcs[!is.na(pcs) & pcs > 0]
  pcs <- sort(unique(pcs))
  if (length(pcs) < 2) {
    stop("pcs must contain at least two positive integers.")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  sanitize_file_component <- function(x) {
    gsub("[^A-Za-z0-9_\\-]", "_", x)
  }
  umap_opts <- plot_options$umap %||% list()
  summary_opts <- plot_options$summary %||% list()
  save_umap_plot <- isTRUE(save_plots) && (
    is.null(plot_options) || is.null(plot_options$umap) || isTRUE(umap_opts$save)
  )
  save_summary_plot <- isTRUE(save_plots) && (
    is.null(plot_options) || is.null(plot_options$summary) || isTRUE(summary_opts$save)
  )

  save_plot <- function(plot_obj,
                        opts = NULL,
                        default_filename,
                        default_title,
                        sample_id = NULL,
                        default_width = 10,
                        default_height = 8) {
    opts <- opts %||% list()
    format <- tolower(opts$format %||% "jpeg")
    device <- switch(format,
                     jpg = "jpeg",
                     jpeg = "jpeg",
                     png = "png",
                     svg = "svg",
                     format)
    filename_prefix <- tools::file_path_sans_ext(opts$filename_prefix %||% default_filename)
    title_prefix <- opts$title_prefix %||% ""
    plot_title <- if (!is.null(sample_id) && nzchar(title_prefix)) {
      paste(title_prefix, sample_id)
    } else if (!is.null(sample_id)) {
      paste(sample_id, default_title)
    } else if (nzchar(title_prefix)) {
      paste(title_prefix, default_title)
    } else {
      default_title
    }
    plot_obj <- plot_obj + ggplot2::ggtitle(plot_title)
    file_name <- if (!is.null(sample_id)) {
      paste0(filename_prefix, sanitize_file_component(sample_id), ".", format)
    } else {
      paste0(filename_prefix, ".", format)
    }
    ggplot2::ggsave(
      filename = file.path(output_dir, file_name),
      plot = plot_obj,
      device = device,
      width = as.numeric(opts$width %||% default_width),
      height = as.numeric(opts$height %||% default_height)
    )
    invisible(file.path(output_dir, file_name))
  }

  sample_col <- meta[["SampleLabel"]]
  sample_labels <- unique(as.character(sample_col))
  sample_labels <- sample_labels[!is.na(sample_labels) & nzchar(sample_labels)]
  if (length(sample_labels) == 0) {
    stop("No valid SampleLabel values found.")
  }

  extract_best_pK <- function(bcmvn, sample_id) {
    if (!is.data.frame(bcmvn) || nrow(bcmvn) == 0) {
      stop("find.pK() returned an empty or invalid object for sample '", sample_id, "'.")
    }
    required_cols <- c("pK", "BCmetric")
    if (!all(required_cols %in% names(bcmvn))) {
      stop(
        "find.pK() output is missing required columns for sample '",
        sample_id,
        "'. Expected: pK, BCmetric."
      )
    }

    metric <- suppressWarnings(as.numeric(unlist(bcmvn[["BCmetric"]], use.names = FALSE)))
    pk_raw <- unlist(bcmvn[["pK"]], use.names = FALSE)
    pk_values <- suppressWarnings(as.numeric(as.character(pk_raw)))
    if (all(is.na(pk_values))) {
      pk_values <- suppressWarnings(as.numeric(pk_raw))
    }

    finite_idx <- which(is.finite(metric))
    if (length(finite_idx) == 0) {
      stop("find.pK() produced no finite BCmetric values for sample '", sample_id, "'.")
    }

    best_idx <- finite_idx[which.max(metric[finite_idx])]
    best_pK <- pk_values[best_idx]
    if (length(best_pK) != 1 || !is.finite(best_pK)) {
      stop("Could not derive a valid numeric pK for sample '", sample_id, "'.")
    }

    best_pK
  }

  select_classification_col <- function(meta_df, pN, pK, nExp) {
    expected_col <- paste("DF.classifications", pN, pK, nExp, sep = "_")
    if (expected_col %in% names(meta_df)) {
      return(expected_col)
    }

    candidates <- grep("^DF\\.classifications", names(meta_df), value = TRUE)
    if (length(candidates) == 0) {
      return(NULL)
    }

    candidates[length(candidates)]
  }

  get_homotypic_prop <- function(sample_seurat, sample_id) {
    cluster_ids <- as.character(Seurat::Idents(sample_seurat))
    cluster_ids <- cluster_ids[!is.na(cluster_ids) & nzchar(cluster_ids)]
    if (length(unique(cluster_ids)) < 2) {
      warning(
        "Sample '",
        sample_id,
        "' has fewer than two clusters; skipping homotypic adjustment."
      )
      return(0)
    }

    DoubletFinder::modelHomotypic(cluster_ids)
  }

  singlet_list <- list()
  count_table <- list()

  for (sample_id in sample_labels) {
    message("Processing sample: ", sample_id)

    cells_in_sample <- rownames(meta)[!is.na(sample_col) & as.character(sample_col) == sample_id]
    if (length(cells_in_sample) == 0) {
      warning("Skipping sample '", sample_id, "' because no cells were found.")
      next
    }

    sample_seurat <- subset(seurat_obj, cells = cells_in_sample)

    if (ncol(sample_seurat) < 50) {
      warning("Skipping sample '", sample_id, "' because it has fewer than 50 cells.")
      next
    }

    sample_seurat <- Seurat::NormalizeData(sample_seurat, verbose = FALSE)
    sample_seurat <- Seurat::FindVariableFeatures(
      sample_seurat,
      selection.method = "vst",
      nfeatures = 3000,
      verbose = FALSE
    )
    sample_seurat <- Seurat::ScaleData(sample_seurat, verbose = FALSE)
    sample_seurat <- Seurat::RunPCA(
      sample_seurat,
      npcs = max(pcs),
      verbose = FALSE
    )

    pca_dims <- ncol(Seurat::Embeddings(sample_seurat, "pca"))
    pcs_use <- pcs[pcs <= pca_dims]
    if (length(pcs_use) < 2) {
      warning(
        "Skipping sample '",
        sample_id,
        "' because there are not enough PCA dimensions for DoubletFinder."
      )
      next
    }

    sample_seurat <- Seurat::FindNeighbors(sample_seurat, dims = pcs_use, verbose = FALSE)
    sample_seurat <- Seurat::FindClusters(sample_seurat, resolution = 0.1, verbose = FALSE)
    sample_seurat <- Seurat::RunUMAP(sample_seurat, dims = pcs_use, verbose = FALSE)

    homotypic.prop <- get_homotypic_prop(sample_seurat, sample_id)

    sweep.res.list <- tryCatch(
      DoubletFinder::paramSweep(sample_seurat, PCs = pcs_use, sct = sct),
      error = function(e) {
        warning("paramSweep() failed for sample '", sample_id, "': ", e$message)
        NULL
      }
    )
    if (is.null(sweep.res.list)) {
      next
    }

    sweep.stats <- tryCatch(
      DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE),
      error = function(e) {
        warning("summarizeSweep() failed for sample '", sample_id, "': ", e$message)
        NULL
      }
    )
    if (is.null(sweep.stats)) {
      next
    }

    bcmvn <- tryCatch(
      DoubletFinder::find.pK(sweep.stats),
      error = function(e) {
        warning("find.pK() failed for sample '", sample_id, "': ", e$message)
        NULL
      }
    )
    if (is.null(bcmvn)) {
      next
    }

    pK <- tryCatch(
      extract_best_pK(bcmvn, sample_id),
      error = function(e) {
        warning(e$message)
        NA_real_
      }
    )
    if (!is.finite(pK)) {
      next
    }

    nExp_poi <- max(1L, round(doublet_rate * ncol(sample_seurat)))
    nExp_poi.adj <- max(1L, round(nExp_poi * (1 - homotypic.prop)))

    sample_seurat <- tryCatch(
      DoubletFinder::doubletFinder(
        sample_seurat,
        PCs = pcs_use,
        pN = 0.25,
        pK = pK,
        nExp = nExp_poi.adj,
        reuse.pANN = NULL,
        sct = sct
      ),
      error = function(e) {
        warning("doubletFinder() failed for sample '", sample_id, "': ", e$message)
        NULL
      }
    )
    if (is.null(sample_seurat)) {
      next
    }

    meta_df <- sample_seurat[[]]
    df_col <- select_classification_col(meta_df, 0.25, pK, nExp_poi.adj)
    if (is.null(df_col)) {
      warning("DoubletFinder classification column not found for sample '", sample_id, "'.")
      next
    }

    df_values <- meta_df[[df_col]]
    if (is.data.frame(df_values)) {
      df_values <- df_values[[1]]
    }
    df_values <- as.character(df_values)

    if (!any(df_values == "Singlet", na.rm = TRUE)) {
      warning("No Singlet cells found for sample '", sample_id, "'.")
      next
    }

    df_plot <- Seurat::DimPlot(
      sample_seurat,
      reduction = "umap",
      group.by = df_col
    )

    if (save_umap_plot) {
      save_plot(
        plot_obj = df_plot,
        opts = umap_opts,
        default_filename = "doubletfinder_",
        default_title = "DoubletFinder",
        sample_id = sample_id,
        default_width = 12,
        default_height = 12
      )
    }

    singlet_cells <- rownames(meta_df)[df_values == "Singlet"]
    if (length(singlet_cells) == 0) {
      warning("No singlet cells retained for sample '", sample_id, "'.")
      next
    }

    singlet_list[[sample_id]] <- subset(sample_seurat, cells = singlet_cells)
    count_table[[sample_id]] <- table(df_values)
  }

  if (length(singlet_list) == 0) {
    stop("Doublet removal failed for all samples.")
  }

  df_count <- purrr::imap_dfr(
    count_table,
    function(counts, sample_id) {
      data.frame(
        Sample = sample_id,
        Category = names(counts),
        Count = as.integer(counts),
        stringsAsFactors = FALSE
      )
    }
  )

  df_long <- dplyr::group_by(df_count, .data$Sample)
  df_long <- dplyr::mutate(df_long, Percentage = .data$Count / sum(.data$Count))
  df_long <- dplyr::ungroup(df_long)

  bar_plot <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = Sample, y = Count, fill = Category)
  ) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0(Count, " (", round(Percentage * 100, 1), "%)")),
      position = ggplot2::position_dodge(width = 0.9),
      vjust = -0.5,
      size = 3,
      fontface = "bold",
      color = "gray35"
    ) +
    ggplot2::labs(
      title = "Singlet & Doublet Count per Sample",
      x = "Sample",
      y = "Count",
      fill = "Category"
    ) +
    ggplot2::theme_minimal()

  if (save_summary_plot) {
    save_plot(
      plot_obj = bar_plot,
      opts = summary_opts,
      default_filename = "singlet.doublet_count_barplot",
      default_title = "Singlet & Doublet Count per Sample",
      default_width = 8,
      default_height = 6
    )
  }

  all_singlet_cells <- unique(unlist(lapply(singlet_list, colnames), use.names = FALSE))
  singlet_merged <- subset(seurat_obj, cells = all_singlet_cells)

  if (save_rds) {
    dir.create(dirname(rds_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(singlet_merged, rds_path)
  }

  singlet_merged
}
