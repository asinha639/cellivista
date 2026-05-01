options(shiny.maxRequestSize = 25 * 1024^3)

required_packages <- c(
  "shiny",
  "Seurat",
  "Matrix",
  "ggplot2",
  "dplyr",
  "base64enc",
  "patchwork",
  "readr",
  "purrr",
  "future"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

library(shiny)
library(Seurat)
library(ggplot2)
library(dplyr)

`%||%` <- function(a, b) if (!is.null(a)) a else b

app_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
helper_dir <- file.path(app_dir, "R")
helper_files <- c(
  "convert_matrix_to_rds.R",
  "read_h5_to_seurat.R",
  "run_qc_metrics.R",
  "post_qc_filtering_and_plots.R",
  "remove_doublets.R",
  "integrate_seurat_samples.R",
  "run_downstream_analysis.R",
  "run_clusterwise_dgea.R",
  "annotate_clusters.R",
  "plot_genes_feature.R"
)

missing_helpers <- helper_files[!file.exists(file.path(helper_dir, helper_files))]
if (length(missing_helpers) > 0) {
  stop(
    "Missing helper file(s) in R/: ",
    paste(missing_helpers, collapse = ", "),
    call. = FALSE
  )
}

for (f in helper_files) source(file.path(helper_dir, f), local = FALSE)

parse_sample_map <- function(x) {
  x <- trimws(x)
  if (!nzchar(x)) return(NULL)
  parts <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) return(NULL)
  out <- lapply(parts, function(part) {
    kv <- trimws(unlist(strsplit(part, "=", fixed = TRUE)))
    if (length(kv) != 2 || !nzchar(kv[1]) || !nzchar(kv[2])) {
      stop("Sample label map must use format like Sample1=Control,Sample2=BPV")
    }
    kv
  })
  vals <- vapply(out, function(z) z[2], character(1))
  nms <- vapply(out, function(z) z[1], character(1))
  stats::setNames(vals, nms)
}

encode_img <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  mime <- switch(ext, png = "image/png", jpg = "image/jpeg", jpeg = "image/jpeg", svg = "image/svg+xml", "application/octet-stream")
  base64enc::dataURI(file = path, mime = mime)
}

safe_read_rds <- function(fileinfo) {
  if (is.null(fileinfo) || is.null(fileinfo$datapath) || !file.exists(fileinfo$datapath)) {
    stop("Uploaded RDS file not found.")
  }
  obj <- tryCatch(
    readRDS(fileinfo$datapath),
    error = function(e) {
      stop("Failed to read uploaded .rds file: ", e$message, call. = FALSE)
    }
  )
  if (!inherits(obj, "Seurat")) {
    stop("Uploaded .rds file does not contain a Seurat object.", call. = FALSE)
  }
  obj
}

normalize_sample_label_for_seurat <- function(seurat_obj, fallback_label = "Sample1") {
  if (is.null(seurat_obj) || !inherits(seurat_obj, "Seurat")) {
    stop("normalize_sample_label_for_seurat() requires a Seurat object.", call. = FALSE)
  }

  meta <- seurat_obj[[]]
  sample_label_note <- NULL

  if (!"SampleLabel" %in% names(meta)) {
    seurat_obj$SampleLabel <- rep(fallback_label, ncol(seurat_obj))
    sample_label_note <- paste0("SampleLabel was missing and was added as '", fallback_label, "' for all cells.")
  } else {
    sample_labels <- as.character(meta[["SampleLabel"]])
    missing_mask <- is.na(sample_labels) | !nzchar(trimws(sample_labels))

    if (all(missing_mask)) {
      seurat_obj$SampleLabel <- rep(fallback_label, ncol(seurat_obj))
      sample_label_note <- paste0("SampleLabel was present but empty, so '", fallback_label, "' was assigned to all cells.")
    } else {
      if (any(missing_mask)) {
        sample_labels[missing_mask] <- fallback_label
        sample_label_note <- paste0("SampleLabel had missing values; those cells were assigned '", fallback_label, "'.")
      } else {
        sample_label_note <- paste0(
          "SampleLabel found: ",
          paste(unique(sample_labels), collapse = ", "),
          "."
        )
      }
      seurat_obj$SampleLabel <- sample_labels
    }
  }

  list(
    object = seurat_obj,
    sample_label_note = sample_label_note
  )
}

load_seurat_rds_for_upload <- function(fileinfo, fallback_label = "Sample1") {
  seurat_obj <- safe_read_rds(fileinfo)
  normalize_sample_label_for_seurat(seurat_obj, fallback_label = fallback_label)
}

help_texts <- c(
  uploadInputTypeHelp = "Choose between raw 10x count data or a prebuilt Seurat object.",
  h5fileHelp = "Upload a raw 10x Genomics .h5 file compatible with Seurat::Read10X_h5().",
  rdsfileHelp = "Upload a .rds file that contains a valid Seurat object.",
  sampleMapHelp = "Optional relabeling of parsed sample IDs. Example: Sample1=Control,Sample2=BPV.",
  projectNameHelp = "Project name assigned to the Seurat object. Default: SeuratProject.",
  mtPatternHelp = "Regex used to identify mitochondrial genes. Example: ^MT- for human, ^mt- for mouse.",
  splitByHelp = "Metadata column used to split QC violin plots. Usually SampleLabel.",
  ylimitHelp = "Lower and upper y-axis bounds for the percent.mt violin plot.",
  minFeatHelp = "Cells with fewer detected genes than this value are removed.",
  maxFeatHelp = "Cells with more detected genes than this value are removed.",
  maxMTHelp = "Cells with mitochondrial percentage above this value are removed.",
  splitByPostHelp = "Metadata column used to split Post-QC plots. Usually SampleLabel.",
  doubletRateHelp = "Expected fraction of doublets. Typical range is about 0.05 to 0.10.",
  savePlotsDoubletHelp = "Save DoubletFinder UMAPs and the singlet/doublet summary bar plot.",
  doubletPlotTypesHelp = "Choose which doublet-removal plots should be saved.",
  nfeaturesIntegrationHelp = "Number of variable features used during integration.",
  savePlotsIntegrationHelp = "Save selected integration plots as image files.",
  integrationSelectedPlotsHelp = "Choose which integration plots should be saved.",
  dimsUsedHelp = "Principal-component range used for UMAP, neighbors, and clustering.",
  resolutionHelp = "Higher resolution produces more clusters.",
  savePlotsDownstreamHelp = "Save selected downstream plots as image files.",
  selectedPlotsHelp = "Choose which downstream plots should be saved.",
  logfcThresholdHelp = "Minimum log fold-change threshold for marker detection.",
  minPctHelp = "Minimum percent of cells expressing a gene for DGEA.",
  savePlotsDGEAHelp = "Save violin, dot, and/or feature plots for top markers.",
  dgeaPlotTypesHelp = "Select which plot types to generate for marker genes.",
  dgeaPlotFormatHelp = "Output image format for DGEA plots.",
  dgeaTitlePrefixHelp = "Prefix added to DGEA plot titles.",
  dgeaPlotWidthHelp = "Width of saved DGEA plots in inches.",
  dgeaPlotHeightHelp = "Height of saved DGEA plots in inches.",
  annotationInputTypeHelp = "Choose whether annotations come from a CSV file or manual comma-separated text.",
  annotationCSVHelp = "CSV should contain cluster-to-label annotations.",
  annotationTextHelp = "Comma-separated labels in cluster order.",
  savePlotsAnnotateHelp = "Save annotated UMAP, split UMAP, and pie chart images.",
  annotationPlotTypesHelp = "Choose which annotation plots should be saved.",
  genesToPlotHelp = "Comma-separated gene symbols to display with FeaturePlot.",
  saveFeaturePlotsHelp = "Save generated FeaturePlots as image files."
)

qm_icon <- function(id) {
  title_text <- if (!is.null(help_texts[[id]])) unname(help_texts[[id]]) else NULL
  if (is.null(title_text) || !nzchar(title_text)) {
    title_text <- id
  }
  tags$span(
    icon("question-circle"),
    id = id,
    class = "qm-help-icon",
    title = title_text,
    tabindex = "0",
    role = "img",
    `aria-label` = title_text,
    style = "margin-left:6px; cursor:help; color:#1c9ed8; display:inline-flex; align-items:center;"
  )
}

help_label <- function(text, id) tagList(text, qm_icon(id))

make_qc_paths <- function(base_dir) {
  list(
    violin = file.path(base_dir, "pre_qc_violin.jpg"),
    percent_mt = file.path(base_dir, "pre_qc_percent.mt.jpg"),
    scatter = file.path(base_dir, "pre_qc_FeatureScatterplot.jpg")
  )
}

make_postqc_paths <- function(base_dir) {
  list(
    vln = file.path(base_dir, "post_qc_vplot.jpg"),
    density = file.path(base_dir, "post_qc_densityplot.jpg"),
    pie = file.path(base_dir, "cellcount_piechart_post.jpg"),
    seurat = file.path(base_dir, "seurat_post_qc_filtered.rds")
  )
}

make_doublet_plot_paths <- function(base_dir, plot_options = NULL) {
  plot_options <- plot_options %||% list()
  umap_opts <- plot_options$umap %||% list()
  summary_opts <- plot_options$summary %||% list()

  umap_prefix <- tools::file_path_sans_ext(umap_opts$filename_prefix %||% "doubletfinder_")
  umap_format <- tolower(umap_opts$format %||% "jpeg")
  summary_prefix <- tools::file_path_sans_ext(summary_opts$filename_prefix %||% "singlet.doublet_count_barplot")
  summary_format <- tolower(summary_opts$format %||% "jpeg")

  list(
    umap_prefix = umap_prefix,
    umap_format = umap_format,
    summary = file.path(base_dir, paste0(summary_prefix, ".", summary_format))
  )
}

collect_doublet_plot_paths <- function(base_dir, plot_options = NULL) {
  spec <- make_doublet_plot_paths(base_dir, plot_options)
  files <- list.files(base_dir, full.names = TRUE, recursive = FALSE)
  umap_files <- files[
    startsWith(basename(files), spec$umap_prefix) &
      tolower(tools::file_ext(files)) == spec$umap_format
  ]
  umap_files <- umap_files[basename(umap_files) != basename(spec$summary)]
  list(
    umap = sort(umap_files),
    summary = if (file.exists(spec$summary)) spec$summary else NULL
  )
}

make_downstream_plot_paths <- function(base_dir, plot_options = NULL) {
  plot_options <- plot_options %||% list()

  make_plot_path <- function(plot_id, default_filename) {
    opts <- plot_options[[plot_id]]
    if (is.null(opts) || !isTRUE(opts$save)) {
      return(NULL)
    }
    format <- tolower(opts$format %||% "jpeg")
    filename <- tools::file_path_sans_ext(opts$filename %||% default_filename)
    file.path(
      base_dir,
      paste0(
        filename,
        ".",
        format
      )
    )
  }

  list(
    umap_grouped = make_plot_path("umap_grouped", "integrated_umap_grouped"),
    umap_split = make_plot_path("umap_split", "integrated_umap_split"),
    pie_chart = make_plot_path("pie_chart", "cellcount_piechart")
  )
}

make_integration_plot_paths <- function(base_dir, plot_options = NULL) {
  plot_options <- plot_options %||% list()

  make_plot_path <- function(plot_id, default_filename) {
    opts <- plot_options[[plot_id]]
    if (is.null(opts) || !isTRUE(opts$save)) {
      return(NULL)
    }
    format <- tolower(opts$format %||% "jpeg")
    filename <- tools::file_path_sans_ext(opts$filename %||% default_filename)
    file.path(
      base_dir,
      paste0(
        filename,
        ".",
        format
      )
    )
  }

  list(
    umap_grouped = make_plot_path("umap_grouped", "integrated_umap_grouped"),
    umap_split = make_plot_path("umap_split", "integrated_umap_split"),
    pie_chart = make_plot_path("pie_chart", "integration_cellcount_piechart")
  )
}

make_dgea_paths <- function(base_dir) {
  list(
    all_marker = file.path(base_dir, "all_marker.csv"),
    top15_marker = file.path(base_dir, "top15_marker.csv"),
    plot_dir = file.path(base_dir, "Plots")
  )
}

make_annotation_paths <- function(base_dir, plot_options = NULL) {
  plot_options <- plot_options %||% list()

  make_plot_path <- function(plot_id, default_filename) {
    opts <- plot_options[[plot_id]] %||% list()
    if (length(plot_options) > 0 && !isTRUE(opts$save)) {
      return(NULL)
    }
    format <- tolower(opts$format %||% "jpeg")
    filename <- tools::file_path_sans_ext(opts$filename %||% default_filename)
    file.path(base_dir, paste0(filename, ".", format))
  }

  list(
    umap_annotated = make_plot_path("umap_annotated", "umap_annotated"),
    umap_split_by_sample = make_plot_path("umap_split_by_sample", "umap_split_by_sample"),
    celltype_piechart = make_plot_path("celltype_piechart", "celltype_piechart")
  )
}

make_feature_plot_paths <- function(base_dir, genes, plot_options = NULL) {
  plot_options <- plot_options %||% list()
  format <- tolower(plot_options$format %||% "jpeg")
  filename_prefix <- plot_options$filename_prefix %||% "featureplot_"
  genes <- as.character(genes)
  stats::setNames(
    file.path(base_dir, paste0(filename_prefix, sanitize_file_component(genes), ".", format)),
    genes
  )
}

sanitize_file_component <- function(x) {
  gsub("[^A-Za-z0-9_\\-]", "_", x)
}

vignette_path <- file.path(app_dir, "vignette.pdf")

make_export_settings_panel <- function(id_prefix,
                                       heading,
                                       filename_label,
                                       default_filename,
                                       title_label,
                                       default_title,
                                       default_format = "jpeg",
                                       default_width = 10,
                                       default_height = 8) {
  wellPanel(
    h5(strong(heading)),
    textInput(paste0(id_prefix, "_filename"), filename_label, value = default_filename),
    textInput(paste0(id_prefix, "_title"), title_label, value = default_title),
    radioButtons(paste0(id_prefix, "_format"), "Format", choices = c("jpeg", "png", "svg"), selected = default_format, inline = TRUE),
    fluidRow(
      column(6, numericInput(paste0(id_prefix, "_width"), "Width (in)", value = default_width, min = 4, max = 30)),
      column(6, numericInput(paste0(id_prefix, "_height"), "Height (in)", value = default_height, min = 4, max = 30))
    )
  )
}

thumbnail_card <- function(title, path, download_id, button_label, img_height = "140px") {
  if (is.null(path) || !file.exists(path)) return(NULL)
  wellPanel(
    h5(title),
    tags$img(src = encode_img(path), height = img_height, class = "plot-inline", style = "max-width: 100%; object-fit: contain;"),
    downloadButton(download_id, button_label, class = "btn-primary qc-dl")
  )
}

get_plot_settings_ui <- function(plot_id) {
  default_name <- switch(plot_id, umap_grouped = "integrated_umap_grouped", umap_split = "integrated_umap_split", pie_chart = "cellcount_piechart", plot_id)
  default_title <- switch(plot_id, umap_grouped = "UMAP: Grouped by Samples", umap_split = "UMAP: Split by Samples", pie_chart = "Cell Count by Sample", "")
  label <- switch(plot_id, umap_grouped = "UMAP grouped", umap_split = "UMAP split", pie_chart = "Cell count pie chart", plot_id)
  wellPanel(
    h5(strong(label)),
    textInput(paste0(plot_id, "_filename"), "Filename", value = default_name),
    textInput(paste0(plot_id, "_title"), "Title", value = default_title),
    radioButtons(paste0(plot_id, "_format"), "Format", choices = c("jpeg", "png", "svg"), selected = "jpeg", inline = TRUE),
    fluidRow(
      column(6, numericInput(paste0(plot_id, "_width"), "Width (in)", value = 10, min = 4, max = 30)),
      column(6, numericInput(paste0(plot_id, "_height"), "Height (in)", value = 8, min = 4, max = 30))
    )
  )
}

get_integration_plot_settings_ui <- function(plot_id) {
  default_name <- switch(
    plot_id,
    umap_grouped = "integrated_umap_grouped",
    umap_split = "integrated_umap_split",
    pie_chart = "integration_cellcount_piechart",
    plot_id
  )
  default_title <- switch(
    plot_id,
    umap_grouped = "UMAP: Grouped by Samples",
    umap_split = "UMAP: Split by Samples",
    pie_chart = "Cell Count by Sample",
    ""
  )
  label <- switch(
    plot_id,
    umap_grouped = "UMAP grouped",
    umap_split = "UMAP split",
    pie_chart = "Cell count pie chart",
    plot_id
  )
  default_width <- switch(
    plot_id,
    umap_grouped = 12,
    umap_split = 18,
    pie_chart = 10,
    10
  )
  default_height <- switch(
    plot_id,
    umap_grouped = 12,
    umap_split = 12,
    pie_chart = 10,
    8
  )
  wellPanel(
    h5(strong(label)),
    textInput(paste0("integration_", plot_id, "_filename"), "Filename", value = default_name),
    textInput(paste0("integration_", plot_id, "_title"), "Title", value = default_title),
    radioButtons(
      paste0("integration_", plot_id, "_format"),
      "Format",
      choices = c("jpeg", "png", "svg"),
      selected = "jpeg",
      inline = TRUE
    ),
    fluidRow(
      column(6, numericInput(paste0("integration_", plot_id, "_width"), "Width (in)", value = default_width, min = 4, max = 30)),
      column(6, numericInput(paste0("integration_", plot_id, "_height"), "Height (in)", value = default_height, min = 4, max = 30))
    )
  )
}

ui <- fluidPage(
  tags$head(
    tags$link(rel = "shortcut icon", href = "logo.png"),
    tags$link(href = "https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600;700&display=swap", rel = "stylesheet"),
    tags$style(HTML("\n      :root {\n        --app-bg-1: #f5f7fb;\n        --app-bg-2: #eef3f8;\n        --app-bg-3: #ffffff;\n        --glass-bg: rgba(255, 255, 255, 0.62);\n        --glass-border: rgba(255, 255, 255, 0.38);\n        --glass-shadow: 0 10px 30px rgba(0, 0, 0, 0.08);\n        --soft-shadow: 0 8px 24px rgba(15, 23, 42, 0.08);\n        --apple-blue: #0a84ff;\n        --apple-blue-soft: #6bb6ff;\n        --text-primary: #1d1d1f;\n        --text-secondary: #5b6472;\n        --input-border: rgba(15, 23, 42, 0.10);\n      }\n\n      html, body {\n        background: linear-gradient(135deg, #e6f2f8, #ffffff, #d2ecf7) !important;\n        color: var(--text-primary) !important;\n        font-family: 'Inter', sans-serif !important;\n        -webkit-font-smoothing: antialiased;\n        -moz-osx-font-smoothing: grayscale;\n        letter-spacing: -0.01em;\n      }\n\n      body {\n        padding: 10px 14px 24px;\n      }\n\n      .container-fluid {\n        padding-left: 18px;\n        padding-right: 18px;\n      }\n\n      .navbar, .navbar-default {\n        background: rgba(255, 255, 255, 0.72) !important;\n        backdrop-filter: blur(20px);\n        -webkit-backdrop-filter: blur(20px);\n        border: 1px solid rgba(255, 255, 255, 0.35) !important;\n        border-bottom: 1px solid rgba(15, 23, 42, 0.06) !important;\n        box-shadow: 0 8px 28px rgba(15, 23, 42, 0.06);\n        border-radius: 18px;\n        margin-bottom: 14px;\n      }\n\n      .navbar .navbar-brand, .navbar-default .navbar-nav > li > a {\n        color: var(--text-primary) !important;\n        font-size: 16px;\n        font-weight: 600;\n        letter-spacing: 0.3px;\n        transition: all 0.2s ease-in-out;\n      }\n\n      .navbar-default .navbar-brand:hover,\n      .navbar-default .navbar-nav > li > a:hover {\n        color: var(--apple-blue) !important;\n      }\n\n      .navbar-default .navbar-nav > li > a {\n        border-radius: 999px;\n        margin: 8px 4px;\n        padding: 10px 14px;\n      }\n\n      .navbar-default .navbar-nav > .active > a,\n      .navbar-default .navbar-nav > .active > a:hover,\n      .navbar-default .navbar-nav > .active > a:focus {\n        background: linear-gradient(180deg, rgba(10, 132, 255, 0.16), rgba(10, 132, 255, 0.10)) !important;\n        color: var(--apple-blue) !important;\n        border: 1px solid rgba(10, 132, 255, 0.18);\n        box-shadow: inset 0 1px 0 rgba(255, 255, 255, 0.55);\n      }\n\n      .tab-content,\n      .tab-pane {\n        padding-top: 8px;\n      }\n\n      .sidebarPanel, .mainPanel, .well, .wellPanel, .panel, .panel-default {\n        background: var(--glass-bg) !important;\n        backdrop-filter: blur(20px);\n        -webkit-backdrop-filter: blur(20px);\n        border: 1px solid var(--glass-border) !important;\n        border-radius: 16px !important;\n        box-shadow: var(--glass-shadow) !important;\n      }\n\n      .sidebarPanel, .mainPanel {\n        padding: 18px 18px 16px !important;\n      }\n\n      .sidebarPanel > .form-group,\n      .mainPanel > .form-group {\n        margin-bottom: 16px;\n      }\n\n      .btn, .btn-default, .btn-primary, .btn-success, .btn-info, .btn-warning, .btn-danger, .btn-file {\n        border-radius: 12px !important;\n        border: 1px solid rgba(10, 132, 255, 0.18) !important;\n        background: linear-gradient(180deg, rgba(123, 196, 255, 0.98) 0%, rgba(10, 132, 255, 0.96) 100%) !important;\n        color: #fff !important;\n        box-shadow: 0 8px 18px rgba(10, 132, 255, 0.20);\n        transition: all 0.2s ease-in-out !important;\n        text-shadow: none !important;\n      }\n\n      .btn-primary {\n        font-weight: 600;\n        letter-spacing: 0.3px;\n      }\n\n      .btn:hover, .btn-default:hover, .btn-primary:hover, .btn-success:hover, .btn-info:hover, .btn-warning:hover, .btn-danger:hover, .btn-file:hover {\n        transform: translateY(-1px);\n        box-shadow: 0 12px 24px rgba(10, 132, 255, 0.26);\n        filter: brightness(1.01);\n      }\n\n      .btn:focus, .btn:active:focus, .btn-default:focus, .btn-primary:focus, .btn-success:focus, .btn-info:focus, .btn-warning:focus, .btn-danger:focus, .btn-file:focus {\n        outline: none !important;\n        box-shadow: 0 0 0 4px rgba(10, 132, 255, 0.14), 0 8px 18px rgba(10, 132, 255, 0.20) !important;\n      }\n\n      .form-control,\n      .selectize-input,\n      .selectize-control.single .selectize-input,\n      .selectize-control.multi .selectize-input,\n      input[type='text'],\n      input[type='number'],\n      input[type='password'],\n      input[type='email'],\n      input[type='url'],\n      textarea {\n        border-radius: 12px !important;\n        border: 1px solid var(--input-border) !important;\n        background: rgba(255, 255, 255, 0.78) !important;\n        box-shadow: inset 0 1px 2px rgba(15, 23, 42, 0.04);\n        transition: all 0.2s ease-in-out !important;\n        color: var(--text-primary) !important;\n      }\n\n      .form-control:focus,\n      .selectize-input.focus,\n      .selectize-control.single .selectize-input.focus,\n      .selectize-control.multi .selectize-input.focus,\n      input[type='text']:focus,\n      input[type='number']:focus,\n      input[type='password']:focus,\n      input[type='email']:focus,\n      input[type='url']:focus,\n      textarea:focus {\n        border-color: rgba(10, 132, 255, 0.45) !important;\n        box-shadow: 0 0 0 4px rgba(10, 132, 255, 0.12), inset 0 1px 2px rgba(15, 23, 42, 0.04) !important;\n        outline: none !important;\n      }\n\n      .checkbox, .radio {\n        transition: all 0.2s ease-in-out;\n      }\n\n      .nav-tabs {\n        border-bottom: 1px solid rgba(15, 23, 42, 0.08);\n        padding-bottom: 4px;\n        gap: 8px;\n      }\n\n      .nav-tabs > li {\n        margin-bottom: 0;\n        margin-right: 6px;\n      }\n\n      .nav-tabs > li > a {\n        border-radius: 999px !important;\n        border: 1px solid transparent !important;\n        color: var(--text-secondary) !important;\n        background: rgba(255, 255, 255, 0.52) !important;\n        box-shadow: none !important;\n        transition: all 0.2s ease-in-out;\n        margin-right: 0;\n      }\n\n      .tabbable .nav-tabs > li > a {\n        color: #1c9ed8 !important;\n        font-size: 15px;\n        font-weight: 600;\n        border-radius: 8px;\n        padding: 10px 16px;\n        transition: all 0.2s ease-in-out;\n      }\n\n      .tabbable .nav-tabs > li > a:hover {\n        background-color: #e8f4fb !important;\n        color: #157bb1 !important;\n      }\n\n      .tabbable .nav-tabs > .active > a,\n      .tabbable .nav-tabs > .active > a:hover,\n      .tabbable .nav-tabs > .active > a:focus {\n        background-color: #1c9ed8 !important;\n        color: white !important;\n        font-weight: 700;\n        box-shadow: 0 2px 6px rgba(0,0,0,0.1);\n      }\n\n      .nav-tabs > li > a:hover {\n        background: rgba(10, 132, 255, 0.08) !important;\n        color: var(--apple-blue) !important;\n      }\n\n      .nav-tabs > li.active > a,\n      .nav-tabs > li.active > a:hover,\n      .nav-tabs > li.active > a:focus {\n        background: linear-gradient(180deg, rgba(10, 132, 255, 0.18), rgba(10, 132, 255, 0.10)) !important;\n        color: var(--apple-blue) !important;\n        border-color: rgba(10, 132, 255, 0.18) !important;\n        box-shadow: inset 0 1px 0 rgba(255, 255, 255, 0.55);\n      }\n\n      .qc-dl {\n        width: 100% !important;\n        margin-bottom: 8px;\n      }\n\n      .plot-inline {\n        margin: 10px 10px 10px 0;\n        max-width: 100%;\n        border-radius: 14px;\n        box-shadow: 0 10px 24px rgba(15, 23, 42, 0.10);\n      }\n\n      .qm-help-icon {\n        opacity: 0.7;\n        transition: all 0.2s ease-in-out;\n      }\n\n      .qm-help-icon:hover,\n      .qm-help-icon:focus {\n        opacity: 1;\n        transform: translateY(-1px) scale(1.03);\n      }\n\n      .qm-help-icon .fa {\n        pointer-events: none;\n      }\n\n      .navbar img,\n      img[alt='logo'] {\n        filter: drop-shadow(0 8px 18px rgba(15, 23, 42, 0.12));\n        transition: all 0.2s ease-in-out;\n      }\n\n      .navbar img:hover,\n      img[alt='logo']:hover {\n        transform: scale(1.02);\n      }\n\n      h1, h2, h3, h4, h5, h6 {\n        color: var(--text-primary);\n        letter-spacing: -0.02em;\n      }\n\n      .form-group label,\n      .control-label,\n      .shiny-input-container label {\n        color: var(--text-secondary);\n        font-weight: 600;\n      }\n\n      * {\n        transition: all 0.2s ease-in-out;\n      }\n    "))
  ),
  tags$style(HTML("\n      :root {\n        --cv-bg-top: #eef6fc;\n        --cv-bg-mid: #f7fbff;\n        --cv-bg-bottom: #eaf4fb;\n        --cv-surface: rgba(255, 255, 255, 0.76);\n        --cv-surface-strong: rgba(255, 255, 255, 0.9);\n        --cv-border: rgba(148, 177, 205, 0.34);\n        --cv-border-soft: rgba(148, 177, 205, 0.22);\n        --cv-shadow: 0 14px 34px rgba(31, 76, 123, 0.1);\n        --cv-shadow-soft: 0 8px 22px rgba(31, 76, 123, 0.08);\n        --cv-blue: #0a84ff;\n        --cv-blue-strong: #0b72df;\n        --cv-blue-soft: #7cc0ff;\n        --cv-text: #1d2733;\n        --cv-text-muted: #5f7083;\n        --cv-focus: rgba(10, 132, 255, 0.18);\n      }\n\n      html,\n      body {\n        min-height: 100%;\n        background: linear-gradient(152deg, var(--cv-bg-top) 0%, var(--cv-bg-mid) 48%, var(--cv-bg-bottom) 100%) !important;\n        color: var(--cv-text) !important;\n      }\n\n      body::before {\n        content: \"\";\n        position: fixed;\n        inset: -14% -8% auto;\n        height: 42vh;\n        background: radial-gradient(ellipse at 20% 0%, rgba(132, 198, 255, 0.2) 0%, rgba(132, 198, 255, 0) 66%), radial-gradient(ellipse at 86% 6%, rgba(155, 210, 255, 0.2) 0%, rgba(155, 210, 255, 0) 64%);\n        pointer-events: none;\n        z-index: 0;\n      }\n\n      .container-fluid {\n        position: relative;\n        z-index: 1;\n      }\n\n      .navbar,\n      .navbar-default {\n        border-radius: 18px;\n        border: 1px solid var(--cv-border) !important;\n        background: var(--cv-surface) !important;\n        box-shadow: var(--cv-shadow-soft);\n      }\n\n      @supports ((backdrop-filter: blur(14px)) or (-webkit-backdrop-filter: blur(14px))) {\n        .navbar,\n        .navbar-default {\n          backdrop-filter: blur(14px) saturate(140%);\n          -webkit-backdrop-filter: blur(14px) saturate(140%);\n        }\n      }\n\n      .navbar-default .navbar-nav > li > a {\n        margin: 8px 5px;\n        padding: 10px 16px;\n        border-radius: 999px;\n        font-size: 14px;\n        font-weight: 600;\n        color: var(--cv-text-muted) !important;\n        border: 1px solid transparent !important;\n        transition: all 180ms ease;\n      }\n\n      .navbar-default .navbar-nav > li > a:hover {\n        background: rgba(10, 132, 255, 0.08) !important;\n        color: var(--cv-blue) !important;\n        transform: translateY(-1px);\n      }\n\n      .navbar-default .navbar-nav > .active > a,\n      .navbar-default .navbar-nav > .active > a:hover,\n      .navbar-default .navbar-nav > .active > a:focus {\n        background: linear-gradient(180deg, rgba(255, 255, 255, 0.88) 0%, rgba(225, 241, 255, 0.92) 100%) !important;\n        color: var(--cv-blue-strong) !important;\n        border: 1px solid rgba(10, 132, 255, 0.22) !important;\n        box-shadow: inset 0 1px 0 rgba(255, 255, 255, 0.85), 0 6px 16px rgba(10, 132, 255, 0.12);\n      }\n\n      .nav-tabs {\n        border-bottom: 1px solid var(--cv-border-soft);\n        padding-bottom: 6px;\n        margin-bottom: 14px;\n      }\n\n      .nav-tabs > li {\n        margin-right: 7px;\n      }\n\n      .nav-tabs > li > a,\n      .tabbable .nav-tabs > li > a {\n        border-radius: 999px !important;\n        border: 1px solid transparent !important;\n        background: rgba(255, 255, 255, 0.68) !important;\n        color: var(--cv-text-muted) !important;\n        font-weight: 600;\n        padding: 9px 15px;\n        transition: all 180ms ease;\n      }\n\n      .nav-tabs > li > a:hover,\n      .tabbable .nav-tabs > li > a:hover {\n        background: rgba(10, 132, 255, 0.08) !important;\n        border-color: rgba(10, 132, 255, 0.12) !important;\n        color: var(--cv-blue-strong) !important;\n        transform: translateY(-1px);\n      }\n\n      .nav-tabs > li.active > a,\n      .nav-tabs > li.active > a:hover,\n      .nav-tabs > li.active > a:focus,\n      .tabbable .nav-tabs > .active > a,\n      .tabbable .nav-tabs > .active > a:hover,\n      .tabbable .nav-tabs > .active > a:focus {\n        color: var(--cv-blue-strong) !important;\n        background: linear-gradient(180deg, rgba(255, 255, 255, 0.92) 0%, rgba(227, 243, 255, 0.95) 100%) !important;\n        border-color: rgba(10, 132, 255, 0.22) !important;\n        box-shadow: inset 0 1px 0 rgba(255, 255, 255, 0.85), 0 8px 16px rgba(10, 132, 255, 0.1);\n      }\n\n      .sidebarPanel,\n      .mainPanel,\n      .wellPanel,\n      .well,\n      .panel,\n      .panel-default {\n        background: var(--cv-surface) !important;\n        border: 1px solid var(--cv-border) !important;\n        border-radius: 18px !important;\n        box-shadow: var(--cv-shadow) !important;\n        padding: 18px 18px 16px !important;\n      }\n\n      @supports ((backdrop-filter: blur(14px)) or (-webkit-backdrop-filter: blur(14px))) {\n        .sidebarPanel,\n        .mainPanel,\n        .wellPanel,\n        .well,\n        .panel,\n        .panel-default {\n          backdrop-filter: blur(14px) saturate(135%);\n          -webkit-backdrop-filter: blur(14px) saturate(135%);\n        }\n      }\n\n      .wellPanel {\n        margin-bottom: 16px;\n      }\n\n      .sidebarPanel > .form-group,\n      .mainPanel > .form-group,\n      .wellPanel > .form-group,\n      .shiny-input-container {\n        margin-bottom: 14px;\n      }\n\n      .form-control,\n      .selectize-input,\n      .selectize-control.single .selectize-input,\n      .selectize-control.multi .selectize-input,\n      input[type='text'],\n      input[type='number'],\n      input[type='password'],\n      input[type='email'],\n      input[type='url'],\n      textarea {\n        min-height: 40px;\n        border-radius: 12px !important;\n        border: 1px solid rgba(130, 158, 185, 0.45) !important;\n        background: rgba(255, 255, 255, 0.84) !important;\n        color: var(--cv-text) !important;\n        box-shadow: inset 0 1px 2px rgba(31, 76, 123, 0.06);\n        transition: all 180ms ease;\n      }\n\n      .form-control:focus,\n      .selectize-input.focus,\n      input[type='text']:focus,\n      input[type='number']:focus,\n      input[type='password']:focus,\n      input[type='email']:focus,\n      input[type='url']:focus,\n      textarea:focus {\n        border-color: rgba(10, 132, 255, 0.42) !important;\n        box-shadow: 0 0 0 3px var(--cv-focus), inset 0 1px 2px rgba(31, 76, 123, 0.06) !important;\n        background: rgba(255, 255, 255, 0.94) !important;\n      }\n\n      .radio label,\n      .checkbox label,\n      .shiny-input-container label,\n      .control-label,\n      .form-group label {\n        color: var(--cv-text-muted) !important;\n        font-weight: 600;\n        letter-spacing: 0.01em;\n      }\n\n      .radio input[type='radio'],\n      .checkbox input[type='checkbox'] {\n        accent-color: var(--cv-blue);\n      }\n\n      .btn,\n      .btn-default,\n      .btn-primary,\n      .btn-success,\n      .btn-info,\n      .btn-warning,\n      .btn-danger,\n      .btn-file {\n        border-radius: 11px !important;\n        border: 1px solid rgba(10, 132, 255, 0.24) !important;\n        background: linear-gradient(180deg, #5fb2ff 0%, #0a84ff 100%) !important;\n        color: #ffffff !important;\n        font-weight: 600;\n        letter-spacing: 0.01em;\n        box-shadow: 0 8px 18px rgba(10, 132, 255, 0.2);\n        transition: all 180ms ease !important;\n      }\n\n      .btn:hover,\n      .btn:focus,\n      .btn-default:hover,\n      .btn-primary:hover,\n      .btn-success:hover,\n      .btn-info:hover,\n      .btn-warning:hover,\n      .btn-danger:hover,\n      .btn-file:hover {\n        transform: translateY(-1px);\n        box-shadow: 0 12px 24px rgba(10, 132, 255, 0.24);\n        filter: saturate(1.04);\n      }\n\n      .btn:focus,\n      .btn-default:focus,\n      .btn-primary:focus,\n      .btn-success:focus,\n      .btn-info:focus,\n      .btn-warning:focus,\n      .btn-danger:focus,\n      .btn-file:focus {\n        box-shadow: 0 0 0 3px var(--cv-focus), 0 12px 24px rgba(10, 132, 255, 0.24) !important;\n      }\n\n      .btn-file {\n        overflow: hidden;\n      }\n\n      .qc-dl {\n        width: auto !important;\n        min-width: 150px;\n        margin-top: 10px;\n        margin-bottom: 2px;\n        padding: 8px 14px !important;\n      }\n\n      .well .qc-dl,\n      .wellPanel .qc-dl,\n      .panel .qc-dl {\n        display: inline-flex;\n        align-items: center;\n        justify-content: center;\n      }\n\n      .plot-inline {\n        margin: 10px 0 8px;\n        border-radius: 14px;\n        border: 1px solid var(--cv-border-soft);\n        background: var(--cv-surface-strong);\n        box-shadow: 0 10px 24px rgba(31, 76, 123, 0.1);\n      }\n\n      .table {\n        background: rgba(255, 255, 255, 0.9);\n        border: 1px solid var(--cv-border-soft);\n        border-radius: 14px;\n        overflow: hidden;\n        box-shadow: 0 8px 18px rgba(31, 76, 123, 0.07);\n      }\n\n      .table > thead > tr > th {\n        background: rgba(10, 132, 255, 0.08);\n        color: var(--cv-text);\n        border-bottom: 1px solid var(--cv-border-soft) !important;\n        font-weight: 700;\n      }\n\n      .table > tbody > tr > td {\n        border-top: 1px solid rgba(130, 158, 185, 0.18) !important;\n      }\n\n      .table-hover > tbody > tr:hover {\n        background: rgba(10, 132, 255, 0.05);\n      }\n\n      .help-block,\n      .text-muted,\n      .shiny-text-output,\n      .shiny-html-output {\n        color: var(--cv-text-muted);\n      }\n\n      .alert,\n      .well pre,\n      pre {\n        border-radius: 12px;\n        border: 1px solid var(--cv-border-soft);\n        background: rgba(255, 255, 255, 0.84);\n      }\n\n      .qm-help-icon {\n        width: 20px;\n        height: 20px;\n        border-radius: 999px;\n        background: rgba(10, 132, 255, 0.1);\n        border: 1px solid rgba(10, 132, 255, 0.2);\n        color: var(--cv-blue) !important;\n        opacity: 0.9;\n        justify-content: center;\n        transition: all 180ms ease;\n      }\n\n      .qm-help-icon:hover,\n      .qm-help-icon:focus {\n        opacity: 1;\n        background: rgba(10, 132, 255, 0.16);\n        transform: translateY(-1px);\n      }\n\n      hr {\n        border-top: 1px solid var(--cv-border-soft);\n      }\n\n      @media (max-width: 992px) {\n        .sidebarPanel,\n        .mainPanel,\n        .wellPanel,\n        .well,\n        .panel,\n        .panel-default {\n          border-radius: 16px !important;\n          padding: 16px 16px 14px !important;\n        }\n\n        .qc-dl {\n          width: 100% !important;\n        }\n      }\n    ")),
  tags$div(style = "text-align:center; margin: 10px 0;", tags$img(src = "logo.png", height = "100px")),
  tags$style(HTML("\n      .floating-guide-widget {
        position: fixed;
        right: 18px;
        bottom: 18px;
        z-index: 9999;
        width: 140px;
        overflow: hidden;
        border-radius: 22px;
        background: rgba(255, 255, 255, 0.86);
        border: 1px solid rgba(15, 23, 42, 0.10);
        box-shadow: 0 12px 30px rgba(15, 23, 42, 0.12);
        backdrop-filter: blur(18px);
        -webkit-backdrop-filter: blur(18px);
        transition: width 0.34s cubic-bezier(0.22, 1, 0.36, 1), transform 0.34s cubic-bezier(0.22, 1, 0.36, 1), box-shadow 0.34s cubic-bezier(0.22, 1, 0.36, 1), background 0.34s ease;
      }

      .floating-guide-widget:hover,
      .floating-guide-widget:focus-within {
        width: 360px;
        transform: translateY(-2px);
        background: rgba(255, 255, 255, 0.94);
        box-shadow: 0 18px 42px rgba(15, 23, 42, 0.16);
      }

      .floating-guide-widget .guide-shell {
        display: flex;
        flex-direction: column;
        gap: 10px;
        padding: 12px 14px;
        min-height: 64px;
      }

      .floating-guide-widget .guide-header {
        display: flex;
        align-items: center;
        gap: 12px;
        min-width: 0;
      }

      .floating-guide-widget .guide-icon {
        flex: 0 0 auto;
        width: 40px;
        height: 40px;
        border-radius: 14px;
        display: inline-flex;
        align-items: center;
        justify-content: center;
        background: linear-gradient(180deg, rgba(34, 149, 255, 0.98) 0%, rgba(10, 132, 255, 0.98) 100%);
        color: #fff;
        box-shadow: 0 10px 18px rgba(10, 132, 255, 0.22);
      }

      .floating-guide-widget .guide-copy {
        min-width: 0;
        display: flex;
        flex-direction: column;
        justify-content: center;
      }

      .floating-guide-widget .guide-kicker {
        display: block;
        font-size: 11px;
        line-height: 1;
        font-weight: 700;
        letter-spacing: 0.08em;
        text-transform: uppercase;
        color: #0a84ff;
        margin-bottom: 4px;
        opacity: 0.88;
      }

      .floating-guide-widget .guide-title {
        display: block;
        font-size: 16px;
        font-weight: 700;
        line-height: 1.12;
        letter-spacing: -0.02em;
        color: #1b1f24;
      }

      .floating-guide-widget .guide-body {
        margin-left: 52px;
        max-height: 0;
        opacity: 0;
        overflow: hidden;
        transform: translateY(-4px);
        transition: max-height 0.30s cubic-bezier(0.22, 1, 0.36, 1), opacity 0.18s ease, transform 0.18s ease, margin-top 0.18s ease;
      }

      .floating-guide-widget:hover .guide-body,
      .floating-guide-widget:focus-within .guide-body {
        max-height: 160px;
        opacity: 1;
        transform: translateY(0);
      }

      .floating-guide-widget .guide-text {
        display: block;
        margin-top: 0;
        font-size: 12px;
        line-height: 1.45;
        color: #5f6875;
      }

      .floating-guide-widget .guide-download {
        display: inline-flex;
        align-items: center;
        justify-content: center;
        margin-top: 10px;
        padding: 8px 14px;
        border-radius: 999px;
        font-size: 12px;
        font-weight: 600;
        text-decoration: none !important;
        color: #fff !important;
        background: linear-gradient(180deg, rgba(34, 149, 255, 0.98) 0%, rgba(10, 132, 255, 0.98) 100%);
        border: 1px solid rgba(10, 132, 255, 0.18);
        box-shadow: 0 8px 18px rgba(10, 132, 255, 0.18);
        opacity: 0;
        transform: translateY(-4px);
        pointer-events: none;
        transition: opacity 0.18s ease, transform 0.18s ease, box-shadow 0.18s ease;
      }

      .floating-guide-widget:hover .guide-download,
      .floating-guide-widget:focus-within .guide-download {
        opacity: 1;
        transform: translateY(0);
        pointer-events: auto;
      }

      .floating-guide-widget .guide-download:hover,
      .floating-guide-widget .guide-download:focus {
        box-shadow: 0 12px 24px rgba(10, 132, 255, 0.26);
        transform: translateY(-1px);
        outline: none;
      }

      @media (max-width: 768px) {
        .floating-guide-widget {
          right: 12px;
          bottom: 12px;
          width: 130px;
        }

        .floating-guide-widget:hover,
        .floating-guide-widget:focus-within {
          width: min(330px, calc(100vw - 24px));
        }

        .floating-guide-widget .guide-body {
          margin-left: 0;
        }
      }
    ")),
  tags$div(
    class = "floating-guide-widget",
    tags$div(
      class = "guide-shell",
      tags$div(
        class = "guide-header",
        tags$div(class = "guide-icon", icon("book")),
        tags$div(
          class = "guide-copy",
          tags$span(class = "guide-kicker", "Guide"),
          tags$span(class = "guide-title", "Need help?")
        )
      ),
      tags$div(
        class = "guide-body",
        tags$span(class = "guide-text", "Try the vignette/user guide and download the example dataset for a quick Cellivista walkthrough."),
        downloadLink(
          "downloadVignette",
          "Download PDF",
          class = "guide-download"
        ),
        downloadLink(
          "downloadExampleDataset",
          "Download example dataset",
          class = "guide-download"
        )
      )
    )
  ),
  navbarPage(
    title = NULL,
    id = "mainNav",
    tabPanel(
      tags$span("Upload Data", class = "no-rds"),
      tabsetPanel(
        id = "uploadSubTabs",
        selected = "Load Data",
        tabPanel(
          "Load Data",
          sidebarLayout(
            sidebarPanel(
              radioButtons(
                "uploadInputType",
                label = help_label("Input type", "uploadInputTypeHelp"),
                choices = c("10x .h5 file" = "h5", "Seurat .rds object" = "rds"),
                selected = "h5",
                inline = TRUE
              ),
              tags$small(
                class = "text-muted",
                "Use .h5 for raw 10x data. Use .rds only if it already contains a valid Seurat object."
              ),
              conditionalPanel(
                condition = "input.uploadInputType == 'h5'",
                fileInput("h5file", label = help_label("Upload 10x .h5 file", "h5fileHelp")),
                textInput("sampleMap", label = help_label("Sample label map (e.g. Sample1=Control, Sample2=BPV)", "sampleMapHelp"), value = "Sample1=Sample1, Sample2=Sample2", placeholder = "Sample1=Sample1, Sample2=Sample2"),
                textInput("projectName", label = help_label("Project name", "projectNameHelp"), value = "SeuratProject")
              ),
              conditionalPanel(
                condition = "input.uploadInputType == 'rds'",
                fileInput("rdsfile", label = help_label("Upload Seurat .rds file", "rdsfileHelp")),
                tags$small(class = "text-muted", "The .rds must contain a valid Seurat object.")
              ),
              actionButton("readH5", "Load Data", class = "btn-primary")
            ),
            mainPanel(verbatimTextOutput("uploadStatus"))
          )
        ),
        tabPanel(
          "Convert Matrix to RDS",
          sidebarLayout(
            sidebarPanel(
              fileInput("matrixFile", "Upload count matrix (.mtx or .mtx.gz)"),
              fileInput("cellFile", "Upload cell/barcode file (.tsv, .csv, .txt, or .gz)"),
              fileInput("geneFile", "Upload gene/features file (.tsv, .csv, .txt, or .gz)"),
              textInput("convertProjectName", "Project name", value = "SeuratProject"),
              numericInput("convertMinCells", "Minimum cells per gene", value = 3, min = 0),
              numericInput("convertMinFeatures", "Minimum features per cell", value = 200, min = 0),
              radioButtons(
                "sampleLabelMode",
                "SampleLabel mode",
                choices = c(
                  "Extract from barcode suffix (-1, -2, etc.)" = "barcode_suffix",
                  "Assign one label to all cells" = "single_sample",
                  "Use column from cell/barcode file" = "from_cell_file_column"
                ),
                selected = "barcode_suffix"
              ),
              conditionalPanel(
                condition = "input.sampleLabelMode == 'barcode_suffix' || input.sampleLabelMode == 'single_sample'",
                textInput("singleSampleLabel", "Single/fallback sample label", value = "Sample1")
              ),
              conditionalPanel(
                condition = "input.sampleLabelMode == 'from_cell_file_column'",
                textInput("sampleLabelColumn", "SampleLabel column name or 1-based index", value = "")
              ),
              actionButton("convertMatrixToRDS", "Convert to RDS", class = "btn-primary"),
              tags$hr(),
              conditionalPanel(
                condition = "output.showConvertedRDSDownload",
                downloadButton("downloadConvertedRDS", "Download converted Seurat (.rds)", class = "btn-primary qc-dl")
              )
            ),
            mainPanel(verbatimTextOutput("convertStatus"))
          )
        )
      )
    ),

    tabPanel(
      "Run QC",
      sidebarLayout(
        sidebarPanel(
          textInput("mtPattern", help_label("Mitochondrial gene pattern", "mtPatternHelp"), value = "^MT-"),
          textInput("splitBy", help_label("Metadata column to split violin plots", "splitByHelp"), value = "SampleLabel"),
          sliderInput("ylimit", help_label("Y-limit for percent.mt plot", "ylimitHelp"), min = 0, max = 100, value = c(0, 80)),
          actionButton("runQC", "Run QC", class = "btn-primary"),
          tags$hr(),
          conditionalPanel(
            condition = "output.showQCDownloads",
            tagList(
              h4("Download QC plots"),
              downloadButton("downloadViolin", "Download violin plot", class = "btn-primary qc-dl"),
              downloadButton("downloadPercentMT", "Download percent.mt plot", class = "btn-primary qc-dl"),
              downloadButton("downloadScatter", "Download feature scatter", class = "btn-primary qc-dl")
            )
          )
        ),
        mainPanel(uiOutput("qcPlotsDisplay"))
      )
    ),

    tabPanel(
      "Post-QC Filtering",
      sidebarLayout(
        sidebarPanel(
          numericInput("minFeat", help_label("Minimum features per cell", "minFeatHelp"), value = 200),
          numericInput("maxFeat", help_label("Maximum features per cell", "maxFeatHelp"), value = 7500),
          numericInput("maxMT", help_label("Maximum percent mitochondrial", "maxMTHelp"), value = 25),
          textInput("splitByPost", help_label("Split plots by", "splitByPostHelp"), value = "SampleLabel"),
          actionButton("postQC", "Run Post-QC Filtering", class = "btn-primary"),
          tags$hr(),
          conditionalPanel(
            condition = "output.showPostQCDownloads",
            tagList(
              h4("Download Post-QC outputs"),
              downloadButton("downloadPostQCRDS", "Download filtered Seurat (.rds)", class = "btn-primary qc-dl"),
              downloadButton("downloadPostQCVln", "Download violin plot", class = "btn-primary qc-dl"),
              downloadButton("downloadPostQCDensity", "Download density plot", class = "btn-primary qc-dl"),
              downloadButton("downloadPostQCPie", "Download pie chart", class = "btn-primary qc-dl")
            )
          )
        ),
        mainPanel(verbatimTextOutput("postQCStatus"), uiOutput("postQCPlotsDisplay"))
      )
    ),

    tabPanel(
      "Remove Doublets",
      sidebarLayout(
        sidebarPanel(
          checkboxInput("loadRDSInput", "Load Seurat object from .rds file", value = FALSE),
          conditionalPanel(condition = "input.loadRDSInput == true", fileInput("inputRDSfile", "Select .rds file")),
          numericInput("doubletRate", help_label("Expected doublet rate", "doubletRateHelp"), value = 0.075, min = 0.001, max = 0.5, step = 0.005),
          checkboxInput("savePlotsDoublet", help_label("Save plots", "savePlotsDoubletHelp"), value = TRUE),
          conditionalPanel(
            condition = "input.savePlotsDoublet == true",
            checkboxGroupInput(
              "doubletPlotTypes",
              help_label("Choose plot types to save", "doubletPlotTypesHelp"),
              choices = list(
                "DoubletFinder UMAPs" = "umap",
                "Singlet & doublet count summary" = "summary"
              ),
              selected = c("umap", "summary")
            ),
            uiOutput("doubletPlotSettingsUI")
          ),
          actionButton("removeDoublets", "Run Doublet Removal", class = "btn-primary"),
          tags$hr()
        ),
        mainPanel(verbatimTextOutput("doubletStatus"), uiOutput("doubletDownloadsUI"))
      )
    ),

    tabPanel(
      "Integration",
      sidebarLayout(
        sidebarPanel(
          checkboxInput("loadRDSIntegration", "Load Seurat object from .rds file", value = FALSE),
          conditionalPanel(condition = "input.loadRDSIntegration == true", fileInput("inputRDSIntegration", "Select .rds file")),
          numericInput("nfeaturesIntegration", help_label("Number of integration features", "nfeaturesIntegrationHelp"), value = 3000, min = 500, step = 500),
          checkboxInput("savePlotsIntegration", help_label("Save plots", "savePlotsIntegrationHelp"), value = TRUE),
          tags$hr(),
          conditionalPanel(
            condition = "input.savePlotsIntegration == true",
            checkboxGroupInput(
              "selectedIntegrationPlots",
              help_label("Select integration plots to save", "integrationSelectedPlotsHelp"),
              choices = list(
                "UMAP grouped by samples" = "umap_grouped",
                "UMAP split by samples" = "umap_split",
                "Cell count pie chart" = "pie_chart"
              ),
              selected = c("umap_grouped", "umap_split", "pie_chart")
            ),
            uiOutput("integrationPlotSettingsUI")
          ),
          actionButton("runIntegration", "Run Integration", class = "btn-primary"),
          tags$hr(),
          conditionalPanel(
            condition = "output.showIntegrationDownloads",
            tagList(
              wellPanel(
                h4("Seurat object"),
                downloadButton("downloadIntegrationRDS", "Download integrated Seurat (.rds)", class = "btn-primary qc-dl")
              ),
              uiOutput("integrationSidebarPlotDownloadsUI")
            )
          )
        ),
        mainPanel(
          verbatimTextOutput("integrationStatus"),
          uiOutput("integrationImageDownloadsUI")
        )
      )
    ),

    tabPanel(
      "Downstream Analysis",
      sidebarLayout(
        sidebarPanel(
          checkboxInput("loadRDSDownstream", "Load Seurat object from .rds file", value = FALSE),
          conditionalPanel(condition = "input.loadRDSDownstream == true", fileInput("inputRDSDownstream", "Select .rds file")),
          sliderInput("dimsUsed", help_label("PCA/UMAP dimensions", "dimsUsedHelp"), min = 1, max = 50, value = c(1, 15)),
          numericInput("resolution", help_label("Clustering resolution", "resolutionHelp"), value = 0.1, min = 0.01, step = 0.05),
          checkboxInput("savePlotsDownstream", help_label("Save plots", "savePlotsDownstreamHelp"), value = TRUE),
          tags$hr(),
          conditionalPanel(
            condition = "input.savePlotsDownstream == true",
            checkboxGroupInput(
              "selectedPlots",
              help_label("Select plots to save", "selectedPlotsHelp"),
              choices = list(
                "UMAP grouped by samples" = "umap_grouped",
                "UMAP split by samples" = "umap_split",
                "Cell count pie chart" = "pie_chart"
              ),
              selected = c("umap_grouped", "umap_split", "pie_chart")
            ),
            uiOutput("plotSettingsUI")
          ),
          actionButton("runDownstream", "Run Downstream Analysis", class = "btn-primary")
        ),
        mainPanel(
          verbatimTextOutput("downstreamStatus"),
          uiOutput("downstreamOutputsUI")
        )
      )
    ),

    tabPanel(
      "Cluster-wise DGEA",
      sidebarLayout(
        sidebarPanel(
          checkboxInput("loadRDSDGEA", "Load Seurat object from .rds file", value = FALSE),
          conditionalPanel(condition = "input.loadRDSDGEA == true", fileInput("inputRDSDGEA", "Select .rds file")),
          numericInput("logfcThreshold", help_label("LogFC threshold", "logfcThresholdHelp"), value = 1),
          numericInput("minPct", help_label("Minimum percent expression", "minPctHelp"), value = 0.25, min = 0, max = 1, step = 0.05),
          checkboxInput("savePlotsDGEA", help_label("Save plots", "savePlotsDGEAHelp"), value = TRUE),
          conditionalPanel(
            condition = "input.savePlotsDGEA == true",
            checkboxGroupInput(
              "dgeaPlotTypes",
              help_label("Choose plot types to save", "dgeaPlotTypesHelp"),
              choices = list("Violin plot" = "violin", "Dot plot" = "dot", "One FeaturePlot per gene" = "feature"),
              selected = c("violin", "dot", "feature")
            ),
            uiOutput("dgeaPlotSettingsUI")
          ),
          actionButton("runDGEA", "Run Cluster-wise DGEA", class = "btn-primary")
        ),
        mainPanel(
          verbatimTextOutput("dgeaStatus"),
          wellPanel(
            h4("Results tables"),
            tableOutput("dgeaPreviewTable"),
            uiOutput("dgeaTableDownloadsUI")
          ),
          wellPanel(
            h4("Plots"),
            uiOutput("dgeaPlotDownloadsUI")
          )
        )
      )
    ),

    tabPanel(
      "Annotate Clusters",
      sidebarLayout(
        sidebarPanel(
          checkboxInput("loadRDSAnnotate", "Load Seurat object from .rds file", value = FALSE),
          conditionalPanel(condition = "input.loadRDSAnnotate == true", fileInput("inputRDSAnnotate", "Select Seurat .rds file")),
          radioButtons("annotationInputType", help_label("Annotation input type", "annotationInputTypeHelp"), choices = c("Upload CSV" = "csv", "Enter manually" = "manual"), selected = "csv"),
          conditionalPanel(condition = "input.annotationInputType == 'csv'", fileInput("annotationCSV", help_label("Upload CSV file for annotations", "annotationCSVHelp"))),
          conditionalPanel(condition = "input.annotationInputType == 'manual'", textAreaInput("annotationText", help_label("Enter annotations (comma-separated, order matters)", "annotationTextHelp"), value = "", rows = 5)),
          checkboxInput("savePlotsAnnotate", help_label("Save plots", "savePlotsAnnotateHelp"), value = TRUE),
          conditionalPanel(
            condition = "input.savePlotsAnnotate == true",
            checkboxGroupInput(
              "annotationPlotTypes",
              help_label("Choose plot types to save", "annotationPlotTypesHelp"),
              choices = list(
                "Annotated UMAP" = "umap_annotated",
                "Split-by-sample annotated UMAP" = "umap_split_by_sample",
                "Cell-type composition pie chart" = "celltype_piechart"
              ),
              selected = c("umap_annotated", "umap_split_by_sample", "celltype_piechart")
            ),
            uiOutput("annotationPlotSettingsUI")
          ),
          actionButton("runAnnotation", "Annotate Clusters", class = "btn-primary"),
          tags$hr()
        ),
        mainPanel(
          verbatimTextOutput("annotationStatus"),
          uiOutput("annotationDownloadsUI")
        )
      )
    ),

    tabPanel(
      "FeaturePlot Genes",
      sidebarLayout(
        sidebarPanel(
          checkboxInput("loadRDSFeaturePlot", "Load Seurat object from .rds file", value = FALSE),
          conditionalPanel(condition = "input.loadRDSFeaturePlot == true", fileInput("inputRDSFeaturePlot", "Select Seurat .rds file")),
          textAreaInput("genesToPlot", help_label("Enter genes to plot (comma-separated)", "genesToPlotHelp"), value = "", rows = 5),
          checkboxInput("saveFeaturePlots", help_label("Save FeaturePlots as images", "saveFeaturePlotsHelp"), value = TRUE),
          conditionalPanel(condition = "input.saveFeaturePlots == true", uiOutput("featurePlotSettingsUI")),
          actionButton("runFeaturePlot", "Plot Genes", class = "btn-primary")
        ),
        mainPanel(verbatimTextOutput("featurePlotStatus"), uiOutput("featurePlotDownloadsUI"))
      )
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    seurat = NULL,
    qc_plot_paths = NULL,
    qc_ready = FALSE,
    qc_counter = 0,
    postqc_plot_paths = NULL,
    postqc_rds_path = NULL,
    postqc_ready = FALSE,
    postqc_counter = 0,
    doublet_rds_path = NULL,
    doublet_plot_paths = NULL,
    doublet_download_specs = NULL,
    doublet_ready = FALSE,
    integration_rds_path = NULL,
    integration_plot_paths = NULL,
    integration_download_specs = NULL,
    integration_ready = FALSE,
    downstream_rds_path = NULL,
    downstream_plot_paths = NULL,
    downstream_ready = FALSE,
    dgea_dir = NULL,
    dgea_csv_all = NULL,
    dgea_csv_top = NULL,
    dgea_plot_paths = NULL,
    dgea_plot_groups = NULL,
    dgea_preview = NULL,
    dgea_download_specs = NULL,
    dgea_ready = FALSE,
    annotation_plot_paths = NULL,
    annotation_rds_path = NULL,
    annotation_ready = FALSE,
    feature_plot_paths = NULL,
    feature_download_specs = NULL,
    feature_ready = FALSE,
    converted_rds_path = NULL,
    converted_ready = FALSE
  )

  session_tmp <- file.path(tempdir(), paste0("cellivista_", Sys.getpid(), "_", as.integer(Sys.time())))
  dir.create(session_tmp, recursive = TRUE, showWarnings = FALSE)

  get_current_or_uploaded <- function(load_flag, file_input) {
    if (isTRUE(load_flag)) safe_read_rds(file_input) else { req(rv$seurat); rv$seurat }
  }

  output$convertStatus <- renderText("")
  output$showConvertedRDSDownload <- reactive(
    isTRUE(rv$converted_ready) &&
      !is.null(rv$converted_rds_path) &&
      file.exists(rv$converted_rds_path)
  )
  outputOptions(output, "showConvertedRDSDownload", suspendWhenHidden = FALSE)

  output$plotSettingsUI <- renderUI({
    if (!isTRUE(input$savePlotsDownstream)) return(NULL)
    selected <- input$selectedPlots %||% character(0)
    if (length(selected) == 0) return(NULL)
    do.call(tagList, lapply(selected, get_plot_settings_ui))
  })
  output$integrationPlotSettingsUI <- renderUI({
    if (!isTRUE(input$savePlotsIntegration)) return(NULL)
    selected <- input$selectedIntegrationPlots %||% character(0)
    if (length(selected) == 0) return(NULL)
    do.call(tagList, lapply(selected, get_integration_plot_settings_ui))
  })

  output$showQCDownloads <- reactive(isTRUE(rv$qc_ready))
  outputOptions(output, "showQCDownloads", suspendWhenHidden = FALSE)
  output$showPostQCDownloads <- reactive(isTRUE(rv$postqc_ready))
  outputOptions(output, "showPostQCDownloads", suspendWhenHidden = FALSE)
  output$showDoubletDownloads <- reactive(isTRUE(rv$doublet_ready) && !is.null(rv$doublet_rds_path) && file.exists(rv$doublet_rds_path))
  outputOptions(output, "showDoubletDownloads", suspendWhenHidden = FALSE)
  output$doubletPlotSettingsUI <- renderUI({
    if (!isTRUE(input$savePlotsDoublet)) return(NULL)
    selected <- input$doubletPlotTypes %||% character(0)
    panels <- list()
    if ("umap" %in% selected) {
      panels[[length(panels) + 1]] <- make_export_settings_panel(
        id_prefix = "doublet_umap",
        heading = "DoubletFinder UMAPs",
        filename_label = "Filename prefix",
        default_filename = "doubletfinder_",
        title_label = "Title prefix",
        default_title = "DoubletFinder:",
        default_width = 12,
        default_height = 12
      )
    }
    if ("summary" %in% selected) {
      panels[[length(panels) + 1]] <- make_export_settings_panel(
        id_prefix = "doublet_summary",
        heading = "Singlet & doublet count summary",
        filename_label = "Filename prefix",
        default_filename = "singlet.doublet_count_barplot",
        title_label = "Title prefix",
        default_title = "Doublet summary:",
        default_width = 8,
        default_height = 6
      )
    }
    if (length(panels) == 0) return(NULL)
    do.call(tagList, panels)
  })
  outputOptions(output, "doubletPlotSettingsUI", suspendWhenHidden = FALSE)
  output$doubletDownloadsUI <- renderUI({
    req(rv$doublet_ready)

    rds_section <- if (!is.null(rv$doublet_rds_path) && file.exists(rv$doublet_rds_path)) {
      wellPanel(
        h4("Seurat object"),
        downloadButton("downloadDoubletRDS", "Download filtered Seurat (.rds)", class = "btn-primary qc-dl")
      )
    } else {
      NULL
    }

    plot_paths <- rv$doublet_plot_paths
    umap_cards <- list()
    summary_card <- NULL
    if (!is.null(plot_paths)) {
      if (!is.null(plot_paths$umap) && length(plot_paths$umap) > 0) {
        umap_cards <- lapply(seq_along(plot_paths$umap), function(i) {
          path <- plot_paths$umap[[i]]
          title <- tools::file_path_sans_ext(basename(path))
          button_label <- paste("Download", basename(path))
          download_id <- paste0("downloadDoubletUMAP_", i)
          thumbnail_card(title, path, download_id, button_label)
        })
      }
      if (!is.null(plot_paths$summary) && file.exists(plot_paths$summary)) {
        summary_card <- thumbnail_card(
          "Singlet & doublet count summary",
          plot_paths$summary,
          "downloadDoubletSummary",
          paste("Download", basename(plot_paths$summary)),
          img_height = "180px"
        )
      }
    }
    umap_cards <- umap_cards[!vapply(umap_cards, is.null, logical(1))]

    plot_section <- NULL
    if (length(umap_cards) > 0 || !is.null(summary_card)) {
      plot_blocks <- list()
      if (length(umap_cards) > 0) {
        umap_rows <- split(seq_along(umap_cards), ceiling(seq_along(umap_cards) / 3))
        umap_blocks <- lapply(umap_rows, function(idx) {
          n_cols <- length(idx)
          do.call(fluidRow, lapply(idx, function(i) column(width = floor(12 / n_cols), umap_cards[[i]])))
        })
        plot_blocks[[length(plot_blocks) + 1]] <- tagList(h4("DoubletFinder UMAPs"), do.call(tagList, umap_blocks))
      }
      if (!is.null(summary_card)) {
        plot_blocks[[length(plot_blocks) + 1]] <- tagList(h4("Summary bar plot"), summary_card)
      }
      plot_section <- wellPanel(h4("Plot images"), do.call(tagList, plot_blocks))
    }

    sections <- list(rds_section, plot_section)
    sections <- sections[!vapply(sections, is.null, logical(1))]
    if (length(sections) == 0) return(NULL)
    do.call(tagList, sections)
  })
  outputOptions(output, "doubletDownloadsUI", suspendWhenHidden = FALSE)
  output$showIntegrationDownloads <- reactive(isTRUE(rv$integration_ready) && !is.null(rv$integration_rds_path) && file.exists(rv$integration_rds_path))
  outputOptions(output, "showIntegrationDownloads", suspendWhenHidden = FALSE)
  build_integration_image_cards <- function() {
    if (!isTRUE(rv$integration_ready)) return(NULL)
    specs <- rv$integration_download_specs
    if (is.null(specs) || length(specs) == 0) return(NULL)

    cards <- list()
    for (spec in specs) {
      cards[[length(cards) + 1]] <- thumbnail_card(spec$title, spec$path, spec$id, spec$label)
    }
    cards <- cards[!vapply(cards, is.null, logical(1))]
    if (length(cards) == 0) return(NULL)

    rows <- split(seq_along(cards), ceiling(seq_along(cards) / 3))
    row_blocks <- lapply(rows, function(idx) {
      n_cols <- length(idx)
      do.call(fluidRow, lapply(idx, function(i) column(width = floor(12 / n_cols), cards[[i]])))
    })
    do.call(tagList, row_blocks)
  }
  output$integrationImageDownloadsUI <- renderUI({
    cards <- build_integration_image_cards()
    if (is.null(cards)) return(NULL)
    wellPanel(
      h4("Plot images"),
      cards
    )
  })
  outputOptions(output, "integrationImageDownloadsUI", suspendWhenHidden = FALSE)
  build_integration_sidebar_buttons <- function() {
    if (!isTRUE(rv$integration_ready)) return(NULL)
    specs <- rv$integration_download_specs
    if (is.null(specs) || length(specs) == 0) return(NULL)

    buttons <- list()
    for (spec in specs) {
      sidebar_id <- switch(
        spec$id,
        downloadIntegrationUmapGrouped = "downloadIntegrationUmapGroupedSidebar",
        downloadIntegrationUmapSplit = "downloadIntegrationUmapSplitSidebar",
        downloadIntegrationPieChart = "downloadIntegrationPieChartSidebar",
        NULL
      )
      if (is.null(sidebar_id)) next
      buttons[[length(buttons) + 1]] <- downloadButton(sidebar_id, spec$label, class = "btn-primary qc-dl")
    }
    if (length(buttons) == 0) return(NULL)
    do.call(tagList, buttons)
  }
  output$integrationSidebarDownloadsUI <- renderUI({
    build_integration_sidebar_buttons()
  })
  outputOptions(output, "integrationSidebarDownloadsUI", suspendWhenHidden = FALSE)
  output$integrationSidebarPlotDownloadsUI <- renderUI({
    buttons <- build_integration_sidebar_buttons()
    if (is.null(buttons)) return(NULL)
    wellPanel(
      h4("Plot images"),
      buttons
    )
  })
  outputOptions(output, "integrationSidebarPlotDownloadsUI", suspendWhenHidden = FALSE)
  output$downstreamOutputsUI <- renderUI({
    if (!isTRUE(rv$downstream_ready)) return(NULL)
    sections <- list(
      wellPanel(
        h4("Seurat object"),
        downloadButton("downloadDownstreamRDS", "Download downstream Seurat (.rds)", class = "btn-primary qc-dl")
      )
    )
    plot_paths <- rv$downstream_plot_paths
    plot_files <- if (is.null(plot_paths)) character(0) else unlist(plot_paths, use.names = FALSE)
    plot_files <- plot_files[!is.na(plot_files) & nzchar(plot_files) & file.exists(plot_files)]
    if (length(plot_files) > 0) {
      sections[[length(sections) + 1]] <- wellPanel(
        h4("Plot images"),
        uiOutput("downstreamImageDownloadsUI")
      )
    }
    do.call(tagList, sections)
  })
  outputOptions(output, "downstreamOutputsUI", suspendWhenHidden = FALSE)
  output$downstreamImageDownloadsUI <- renderUI({
    req(rv$downstream_ready)
    plot_paths <- rv$downstream_plot_paths
    if (is.null(plot_paths)) return(NULL)

    cards <- list(
      thumbnail_card("UMAP grouped", plot_paths$umap_grouped, "downloadDownstreamUmapGrouped", "Download UMAP grouped"),
      thumbnail_card("UMAP split", plot_paths$umap_split, "downloadDownstreamUmapSplit", "Download UMAP split"),
      thumbnail_card("Cell count pie chart", plot_paths$pie_chart, "downloadDownstreamPieChart", "Download cell count pie chart")
    )
    cards <- cards[!vapply(cards, is.null, logical(1))]
    if (length(cards) == 0) return(NULL)

    n_cols <- min(3, length(cards))
    columns <- lapply(seq_along(cards), function(i) {
      column(width = floor(12 / n_cols), cards[[i]])
    })
    do.call(fluidRow, columns)
  })
  outputOptions(output, "downstreamImageDownloadsUI", suspendWhenHidden = FALSE)
  output$dgeaPlotSettingsUI <- renderUI({
    req(input$savePlotsDGEA, input$dgeaPlotTypes)
    panels <- lapply(input$dgeaPlotTypes, function(plot_id) {
      switch(
        plot_id,
        violin = make_export_settings_panel(
          id_prefix = "dgea_violin",
          heading = "Violin plots",
          filename_label = "Filename prefix",
          default_filename = "violin_cluster",
          title_label = "Title prefix",
          default_title = "Cluster",
          default_width = 10,
          default_height = 8
        ),
        dot = make_export_settings_panel(
          id_prefix = "dgea_dot",
          heading = "Dot plots",
          filename_label = "Filename prefix",
          default_filename = "dotplot_cluster",
          title_label = "Title prefix",
          default_title = "Cluster",
          default_width = 10,
          default_height = 8
        ),
        feature = make_export_settings_panel(
          id_prefix = "dgea_feature",
          heading = "Feature plots",
          filename_label = "Filename prefix",
          default_filename = "featureplot_cluster",
          title_label = "Title prefix",
          default_title = "Cluster",
          default_width = 10,
          default_height = 8
        ),
        NULL
      )
    })
    panels <- panels[!vapply(panels, is.null, logical(1))]
    if (length(panels) == 0) return(NULL)
    do.call(tagList, panels)
  })
  outputOptions(output, "dgeaPlotSettingsUI", suspendWhenHidden = FALSE)
  output$dgeaPreviewTable <- renderTable({
    req(rv$dgea_ready)
    rv$dgea_preview
  }, striped = TRUE, hover = TRUE, width = "100%")
  output$dgeaTableDownloadsUI <- renderUI({
    req(rv$dgea_ready)
    specs <- rv$dgea_download_specs
    if (is.null(specs) || length(specs) == 0) return(NULL)

    buttons <- list()
    for (spec in specs) {
      if (!identical(spec$type, "table")) next
      if (!is.null(spec$path) && file.exists(spec$path)) {
        buttons[[length(buttons) + 1]] <- downloadButton(spec$id, spec$label, class = "btn-primary qc-dl")
      }
    }
    if (length(buttons) == 0) return(NULL)
    do.call(tagList, buttons)
  })
  outputOptions(output, "dgeaTableDownloadsUI", suspendWhenHidden = FALSE)
  dgea_normalize_prefix <- function(x) {
    tools::file_path_sans_ext(as.character(x %||% ""))
  }
  dgea_allowed_exts <- function(fmt) {
    fmt <- tolower(as.character(fmt %||% ""))
    if (identical(fmt, "jpg") || identical(fmt, "jpeg")) {
      return(c("jpg", "jpeg"))
    }
    fmt
  }
  dgea_filter_group_paths <- function(paths, prefix, fmt) {
    if (is.null(paths) || length(paths) == 0) return(character(0))
    valid_paths <- paths[file.exists(paths)]
    if (length(valid_paths) == 0) return(character(0))
    prefix_norm <- dgea_normalize_prefix(prefix)
    exts <- dgea_allowed_exts(fmt)
    keep <- startsWith(basename(valid_paths), prefix_norm) &
      (tolower(tools::file_ext(valid_paths)) %in% exts)
    valid_paths[keep]
  }
  write_dgea_group_zip <- function(paths, file) {
    valid_paths <- paths[file.exists(paths)]
    if (length(valid_paths) == 0) {
      stop("No files are available for this DGEA plot group.")
    }
    tmp_zip <- tempfile(fileext = ".zip")
    on.exit(unlink(tmp_zip), add = TRUE)
    utils::zip(
      zipfile = tmp_zip,
      files = normalizePath(valid_paths, winslash = "/", mustWork = TRUE),
      flags = "-j"
    )
    file.copy(tmp_zip, file, overwrite = TRUE)
  }
  output$dgeaPlotDownloadsUI <- renderUI({
    req(rv$dgea_ready)
    specs <- rv$dgea_download_specs
    if (is.null(specs) || length(specs) == 0) return(NULL)

    plot_groups <- rv$dgea_plot_groups %||% list()
    zip_buttons <- list()
    if (length(plot_groups$violin %||% character(0)) > 0) {
      zip_buttons[[length(zip_buttons) + 1]] <- downloadButton(
        "downloadDGEAViolinZip",
        "Download all violin plots (.zip)",
        class = "btn-primary qc-dl"
      )
    }
    if (length(plot_groups$dot %||% character(0)) > 0) {
      zip_buttons[[length(zip_buttons) + 1]] <- downloadButton(
        "downloadDGEADotZip",
        "Download all dot plots (.zip)",
        class = "btn-primary qc-dl"
      )
    }
    if (length(plot_groups$feature %||% character(0)) > 0) {
      zip_buttons[[length(zip_buttons) + 1]] <- downloadButton(
        "downloadDGEAFeatureZip",
        "Download all FeaturePlots (.zip)",
        class = "btn-primary qc-dl"
      )
    }
    zip_panel <- if (length(zip_buttons) > 0) {
      wellPanel(
        h5("Download all plot groups"),
        do.call(tagList, zip_buttons)
      )
    } else {
      NULL
    }

    cards <- list()
    for (spec in specs) {
      if (!identical(spec$type, "plot")) next
      cards[[length(cards) + 1]] <- thumbnail_card(spec$title, spec$path, spec$id, spec$label)
    }
    cards <- cards[!vapply(cards, is.null, logical(1))]
    if (length(cards) == 0) return(zip_panel)

    rows <- split(seq_along(cards), ceiling(seq_along(cards) / 3))
    row_blocks <- lapply(rows, function(idx) {
      n_cols <- length(idx)
      do.call(fluidRow, lapply(idx, function(i) column(width = floor(12 / n_cols), cards[[i]])))
    })
    sections <- c(
      if (!is.null(zip_panel)) list(zip_panel) else list(),
      list(do.call(tagList, row_blocks))
    )
    do.call(tagList, sections)
  })
  outputOptions(output, "dgeaPlotDownloadsUI", suspendWhenHidden = FALSE)
  output$downloadDGEAViolinZip <- downloadHandler(
    filename = function() "dgea_violin_plots.zip",
    content = function(file) {
      paths <- rv$dgea_plot_groups$violin %||% character(0)
      write_dgea_group_zip(paths, file)
    }
  )
  output$downloadDGEADotZip <- downloadHandler(
    filename = function() "dgea_dot_plots.zip",
    content = function(file) {
      paths <- rv$dgea_plot_groups$dot %||% character(0)
      write_dgea_group_zip(paths, file)
    }
  )
  output$downloadDGEAFeatureZip <- downloadHandler(
    filename = function() "dgea_featureplots.zip",
    content = function(file) {
      paths <- rv$dgea_plot_groups$feature %||% character(0)
      write_dgea_group_zip(paths, file)
    }
  )
  output$annotationPlotSettingsUI <- renderUI({
    if (!isTRUE(input$savePlotsAnnotate)) return(NULL)
    selected <- input$annotationPlotTypes %||% character(0)
    panels <- list()
    if ("umap_annotated" %in% selected) {
      panels[[length(panels) + 1]] <- make_export_settings_panel(
        id_prefix = "annotation_umap_annotated",
        heading = "Annotated UMAP",
        filename_label = "Filename",
        default_filename = "umap_annotated",
        title_label = "Title",
        default_title = "UMAP: Annotated Cell Types",
        default_width = 10,
        default_height = 8
      )
    }
    if ("umap_split_by_sample" %in% selected) {
      panels[[length(panels) + 1]] <- make_export_settings_panel(
        id_prefix = "annotation_umap_split_by_sample",
        heading = "Split-by-sample annotated UMAP",
        filename_label = "Filename",
        default_filename = "umap_split_by_sample",
        title_label = "Title",
        default_title = "UMAP: Annotated Cell Types by Sample",
        default_width = 18,
        default_height = 8
      )
    }
    if ("celltype_piechart" %in% selected) {
      panels[[length(panels) + 1]] <- make_export_settings_panel(
        id_prefix = "annotation_celltype_piechart",
        heading = "Cell-type composition pie chart",
        filename_label = "Filename",
        default_filename = "celltype_piechart",
        title_label = "Title",
        default_title = "Cell Type Composition",
        default_width = 10,
        default_height = 10
      )
    }
    if (length(panels) == 0) return(NULL)
    do.call(tagList, panels)
  })
  outputOptions(output, "annotationPlotSettingsUI", suspendWhenHidden = FALSE)
  output$annotationDownloadsUI <- renderUI({
    req(rv$annotation_ready)
    plot_paths <- rv$annotation_plot_paths

    rds_section <- if (!is.null(rv$annotation_rds_path) && file.exists(rv$annotation_rds_path)) {
      wellPanel(
        h4("Seurat object"),
        downloadButton("downloadAnnotationRDS", "Download annotated Seurat (.rds)", class = "btn-primary qc-dl")
      )
    } else {
      NULL
    }

    cards <- list(
      thumbnail_card("Annotated UMAP", plot_paths$umap_annotated, "downloadAnnotationUMAP", "Download annotated UMAP"),
      thumbnail_card("Split by sample", plot_paths$umap_split_by_sample, "downloadAnnotationUMAPSplit", "Download split annotated UMAP"),
      thumbnail_card("Cell type composition", plot_paths$celltype_piechart, "downloadAnnotationPieChart", "Download cell type pie chart")
    )
    cards <- cards[!vapply(cards, is.null, logical(1))]

    plot_section <- if (length(cards) > 0) {
      rows <- split(seq_along(cards), ceiling(seq_along(cards) / 3))
      row_blocks <- lapply(rows, function(idx) {
        n_cols <- length(idx)
        do.call(fluidRow, lapply(idx, function(i) column(width = floor(12 / n_cols), cards[[i]])))
      })
      wellPanel(
        h4("Annotation plots"),
        do.call(tagList, row_blocks)
      )
    } else {
      NULL
    }

    sections <- list(rds_section, plot_section)
    sections <- sections[!vapply(sections, is.null, logical(1))]
    if (length(sections) == 0) return(NULL)
    do.call(tagList, sections)
  })
  outputOptions(output, "annotationDownloadsUI", suspendWhenHidden = FALSE)
  output$featurePlotSettingsUI <- renderUI({
    tagList(
      make_export_settings_panel(
        id_prefix = "featureplot",
        heading = "FeaturePlot export",
        filename_label = "Filename prefix",
        default_filename = "featureplot_",
        title_label = "Title prefix",
        default_title = "FeaturePlot:",
        default_width = 7,
        default_height = 6
      )
    )
  })
  outputOptions(output, "featurePlotSettingsUI", suspendWhenHidden = FALSE)
  output$featurePlotDownloadsUI <- renderUI({
    req(rv$feature_ready)
    specs <- rv$feature_download_specs
    if (is.null(specs) || length(specs) == 0) return(NULL)
    paths <- rv$feature_plot_paths %||% character(0)
    existing_paths <- paths[file.exists(paths)]

    cards <- list()
    for (spec in specs) {
      cards[[length(cards) + 1]] <- thumbnail_card(spec$title, spec$path, spec$id, spec$label)
    }
    cards <- cards[!vapply(cards, is.null, logical(1))]
    zip_button <- if (length(existing_paths) > 0) {
      downloadButton("downloadFeaturePlotsZip", "Download all FeaturePlots (.zip)", class = "btn-primary qc-dl")
    } else {
      NULL
    }
    if (length(cards) == 0 && is.null(zip_button)) return(NULL)

    row_blocks <- list()
    if (length(cards) > 0) {
      rows <- split(seq_along(cards), ceiling(seq_along(cards) / 3))
      row_blocks <- lapply(rows, function(idx) {
        n_cols <- length(idx)
        do.call(fluidRow, lapply(idx, function(i) column(width = floor(12 / n_cols), cards[[i]])))
      })
    }

    panel_content <- c(
      if (!is.null(zip_button)) list(zip_button) else list(),
      row_blocks
    )
    wellPanel(
      h4("FeaturePlot downloads"),
      do.call(tagList, panel_content)
    )
  })
  outputOptions(output, "featurePlotDownloadsUI", suspendWhenHidden = FALSE)
  output$downloadFeaturePlotsZip <- downloadHandler(
    filename = function() "featureplots.zip",
    content = function(file) {
      paths <- rv$feature_plot_paths %||% character(0)
      existing_paths <- paths[file.exists(paths)]
      if (length(existing_paths) == 0) {
        stop("No FeaturePlot image files are available to zip.")
      }
      tmp_zip <- tempfile(fileext = ".zip")
      on.exit(unlink(tmp_zip), add = TRUE)
      utils::zip(
        zipfile = tmp_zip,
        files = normalizePath(existing_paths, winslash = "/", mustWork = TRUE),
        flags = "-j"
      )
      file.copy(tmp_zip, file, overwrite = TRUE)
    }
  )
  output$showAnnotationDownloads <- reactive(isTRUE(rv$annotation_ready) && !is.null(rv$annotation_rds_path) && file.exists(rv$annotation_rds_path))
  outputOptions(output, "showAnnotationDownloads", suspendWhenHidden = FALSE)

  output$qcPlotsDisplay <- renderUI({
    req(rv$qc_ready)
    rv$qc_counter
    tagList(
      h4("QC Plots"),
      tags$img(src = encode_img(rv$qc_plot_paths$violin), height = "400px", class = "plot-inline"),
      tags$img(src = encode_img(rv$qc_plot_paths$percent_mt), height = "400px", class = "plot-inline"),
      tags$img(src = encode_img(rv$qc_plot_paths$scatter), height = "400px", class = "plot-inline")
    )
  })

  output$postQCPlotsDisplay <- renderUI({
    req(rv$postqc_ready)
    rv$postqc_counter
    tagList(
      h4("Post-QC Plots"),
      tags$img(src = encode_img(rv$postqc_plot_paths$vln), height = "400px", class = "plot-inline"),
      tags$img(src = encode_img(rv$postqc_plot_paths$density), height = "400px", class = "plot-inline"),
      tags$img(src = encode_img(rv$postqc_plot_paths$pie), height = "400px", class = "plot-inline")
    )
  })

  observeEvent(input$convertMatrixToRDS, {
    output$convertStatus <- renderText("Converting matrix files to Seurat RDS...")
    tryCatch({
      rv$converted_ready <- FALSE
      rv$converted_rds_path <- NULL

      if (is.null(input$matrixFile) || is.null(input$matrixFile$datapath) || !file.exists(input$matrixFile$datapath)) {
        stop("Please upload a valid matrix file (.mtx or .mtx.gz).")
      }
      if (is.null(input$cellFile) || is.null(input$cellFile$datapath) || !file.exists(input$cellFile$datapath)) {
        stop("Please upload a valid cell/barcode file.")
      }
      if (is.null(input$geneFile) || is.null(input$geneFile$datapath) || !file.exists(input$geneFile$datapath)) {
        stop("Please upload a valid gene/features file.")
      }

      conversion_dir <- file.path(session_tmp, "converted_rds")
      dir.create(conversion_dir, recursive = TRUE, showWarnings = FALSE)
      output_path <- file.path(conversion_dir, "converted_seurat.rds")

      sample_column_input <- trimws(input$sampleLabelColumn %||% "")
      sample_label_column <- NULL
      if (identical(input$sampleLabelMode, "from_cell_file_column")) {
        if (!nzchar(sample_column_input)) {
          stop("Please provide SampleLabel column name or 1-based index.")
        }
        if (grepl("^[0-9]+$", sample_column_input)) {
          sample_label_column <- as.integer(sample_column_input)
        } else {
          sample_label_column <- sample_column_input
        }
      }

      seurat_obj <- convert_matrix_to_rds(
        matrix_path = input$matrixFile$datapath,
        cell_path = input$cellFile$datapath,
        gene_path = input$geneFile$datapath,
        output_path = output_path,
        project_name = input$convertProjectName,
        sample_label_mode = input$sampleLabelMode,
        single_sample_label = input$singleSampleLabel,
        sample_label_column = sample_label_column,
        min_cells = input$convertMinCells,
        min_features = input$convertMinFeatures
      )

      rv$seurat <- seurat_obj
      rv$converted_rds_path <- output_path
      rv$converted_ready <- TRUE

      sample_labels <- sort(unique(as.character(seurat_obj$SampleLabel)))
      status_lines <- c(
        "Conversion completed successfully.",
        paste0("Cells: ", ncol(seurat_obj)),
        paste0("Genes/features: ", nrow(seurat_obj)),
        paste0("SampleLabel values: ", paste(sample_labels, collapse = ", ")),
        paste0("Output file: ", basename(output_path))
      )
      output$convertStatus <- renderText(paste(status_lines, collapse = "\n"))
    }, error = function(e) {
      output$convertStatus <- renderText(paste("Conversion failed:", conditionMessage(e)))
    })
  })
  output$downloadConvertedRDS <- downloadHandler(
    filename = function() {
      path <- rv$converted_rds_path
      if (is.null(path)) "converted_seurat.rds" else basename(path)
    },
    content = function(file) {
      path <- rv$converted_rds_path
      if (is.null(path) || !file.exists(path)) {
        stop("Converted Seurat RDS is not available for download.")
      }
      file.copy(path, file, overwrite = TRUE)
    }
  )

  observeEvent(input$readH5, {
    output$uploadStatus <- renderText("")
    upload_type <- input$uploadInputType %||% "h5"

    if (identical(upload_type, "rds")) {
      if (is.null(input$rdsfile) || is.null(input$rdsfile$datapath) || !nzchar(input$rdsfile$datapath)) {
        output$uploadStatus <- renderText("Read failed: no .rds file was uploaded. Please upload a Seurat .rds object first.")
        return()
      }
      output$uploadStatus <- renderText("Reading Seurat .rds object...")
      tryCatch({
        loaded <- load_seurat_rds_for_upload(input$rdsfile, fallback_label = "Sample1")
        rv$seurat <- loaded$object
        status_lines <- c(
          "Loaded Seurat object successfully.",
          paste0("Cells: ", ncol(rv$seurat)),
          paste0("Genes/features: ", nrow(rv$seurat)),
          loaded$sample_label_note
        )
        output$uploadStatus <- renderText(paste(status_lines, collapse = "\n"))
      }, error = function(e) {
        output$uploadStatus <- renderText(paste("Read failed:", e$message))
      })
      return()
    }

    if (is.null(input$h5file) || is.null(input$h5file$datapath) || !nzchar(input$h5file$datapath)) {
      output$uploadStatus <- renderText("Read failed: no .h5 file was uploaded. Please upload a 10x .h5 file first.")
      return()
    }
    output$uploadStatus <- renderText("Reading .h5 file...")
    map_text <- trimws(input$sampleMap)
    sample_map <- NULL
    if (nzchar(map_text)) {
      tryCatch({ sample_map <- parse_sample_map(map_text) }, error = function(e) { output$uploadStatus <- renderText(paste("Read failed: invalid sample map format.", e$message)); sample_map <<- NULL })
      if (is.null(sample_map) && nzchar(map_text)) return()
    }
    tryCatch({
      rv$seurat <- read_h5_to_seurat(filepath = input$h5file$datapath, project_name = input$projectName, sample_label_map = sample_map)
      output$uploadStatus <- renderText(paste0("Read successful.\nCells: ", ncol(rv$seurat), "\nGenes: ", nrow(rv$seurat), "\nSamples detected: ", paste(unique(as.character(rv$seurat$SampleLabel)), collapse = ", ")))
    }, error = function(e) output$uploadStatus <- renderText(paste("Read failed:", e$message)))
  })

  observeEvent(input$runQC, {
    tryCatch({
      req(rv$seurat)
      qc_dir <- file.path(session_tmp, "qc")
      dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
      rv$seurat <- run_qc_metrics(seurat_obj = rv$seurat, mt_pattern = input$mtPattern, output_dir = qc_dir, sample_split_var = input$splitBy, ylimit = input$ylimit, save_plots = TRUE)
      rv$qc_plot_paths <- make_qc_paths(qc_dir)
      if (!all(file.exists(unlist(rv$qc_plot_paths)))) stop("QC plot generation did not produce all expected files.")
      rv$qc_ready <- TRUE
      rv$qc_counter <- rv$qc_counter + 1
      showNotification("QC completed successfully.", type = "message")
    }, error = function(e) showNotification(paste("QC failed:", e$message), type = "error"))
  })

  observeEvent(input$postQC, {
    output$postQCStatus <- renderText("Running post-QC filtering...")
    tryCatch({
      req(rv$seurat)
      postqc_dir <- file.path(session_tmp, "post_qc")
      dir.create(postqc_dir, recursive = TRUE, showWarnings = FALSE)
      paths <- make_postqc_paths(postqc_dir)
      rv$seurat <- post_qc_filtering_and_plots(seurat_obj = rv$seurat, min_features = input$minFeat, max_features = input$maxFeat, max_mt_percent = input$maxMT, output_dir = postqc_dir, sample_split_var = input$splitByPost, save_rds = TRUE, rds_path = paths$seurat)
      rv$postqc_plot_paths <- paths[c("vln", "density", "pie")]
      rv$postqc_rds_path <- paths$seurat
      if (!all(file.exists(unlist(paths)))) stop("Post-QC output generation did not produce all expected files.")
      rv$postqc_ready <- TRUE
      rv$postqc_counter <- rv$postqc_counter + 1
      output$postQCStatus <- renderText("Post-QC filtering completed successfully.")
    }, error = function(e) output$postQCStatus <- renderText(paste("Post-QC failed:", e$message)))
  })

  output$downloadViolin <- downloadHandler(filename = function() "pre_qc_violin.jpg", content = function(file) file.copy(rv$qc_plot_paths$violin, file, overwrite = TRUE))
  output$downloadPercentMT <- downloadHandler(filename = function() "pre_qc_percent.mt.jpg", content = function(file) file.copy(rv$qc_plot_paths$percent_mt, file, overwrite = TRUE))
  output$downloadScatter <- downloadHandler(filename = function() "pre_qc_FeatureScatterplot.jpg", content = function(file) file.copy(rv$qc_plot_paths$scatter, file, overwrite = TRUE))
  output$downloadPostQCRDS <- downloadHandler(filename = function() "seurat_post_qc_filtered.rds", content = function(file) file.copy(rv$postqc_rds_path, file, overwrite = TRUE))
  output$downloadPostQCVln <- downloadHandler(filename = function() "post_qc_vplot.jpg", content = function(file) file.copy(rv$postqc_plot_paths$vln, file, overwrite = TRUE))
  output$downloadPostQCDensity <- downloadHandler(filename = function() "post_qc_densityplot.jpg", content = function(file) file.copy(rv$postqc_plot_paths$density, file, overwrite = TRUE))
  output$downloadPostQCPie <- downloadHandler(filename = function() "cellcount_piechart_post.jpg", content = function(file) file.copy(rv$postqc_plot_paths$pie, file, overwrite = TRUE))
  output$downloadVignette <- downloadHandler(
    filename = function() "vignette.pdf",
    content = function(file) {
      if (!file.exists(vignette_path)) {
        stop("vignette.pdf was not found next to app.R.")
      }
      file.copy(vignette_path, file, overwrite = TRUE)
    }
  )
  output$downloadExampleDataset <- downloadHandler(
    filename = function() "GSE132044_example_dataset.zip",
    content = function(file) {
      example_dir <- file.path(app_dir, "GSE132044")
      required_relative_files <- c(
        file.path("Raw Data", "GSE132044_pbmc_hg38_count_matrix.mtx.gz"),
        file.path("Raw Data", "GSE132044_pbmc_hg38_cell.tsv.gz"),
        file.path("Raw Data", "GSE132044_pbmc_hg38_gene.tsv.gz"),
        "rds_convert.R",
        "GSE132044_pbmc_hg38.rds"
      )
      required_full_paths <- file.path(example_dir, required_relative_files)
      missing_idx <- which(!file.exists(required_full_paths))
      if (length(missing_idx) > 0) {
        missing_files <- required_relative_files[missing_idx]
        stop(
          paste0(
            "Missing required example dataset file(s):\n- ",
            paste(missing_files, collapse = "\n- ")
          )
        )
      }

      staging_parent <- file.path(tempdir(), paste0("GSE132044_example_dataset_", as.integer(Sys.time())))
      staging_root <- file.path(staging_parent, "GSE132044")
      dir.create(staging_root, recursive = TRUE, showWarnings = FALSE)

      for (i in seq_along(required_relative_files)) {
        rel_path <- required_relative_files[[i]]
        src <- required_full_paths[[i]]
        dest <- file.path(staging_root, rel_path)
        dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
        file.copy(src, dest, overwrite = TRUE)
      }

      old_wd <- getwd()
      on.exit(setwd(old_wd), add = TRUE)
      setwd(staging_parent)
      files_to_zip <- list.files("GSE132044", recursive = TRUE, all.files = FALSE)
      files_to_zip <- file.path("GSE132044", files_to_zip)
      utils::zip(zipfile = file, files = files_to_zip)
    }
  )

  observeEvent(input$removeDoublets, {
    output$doubletStatus <- renderText("Running doublet removal...")
    tryCatch({
      seurat_in <- get_current_or_uploaded(input$loadRDSInput, input$inputRDSfile)
      doublet_dir <- file.path(session_tmp, "doublet")
      dir.create(doublet_dir, recursive = TRUE, showWarnings = FALSE)
      old_doublet_files <- list.files(doublet_dir, full.names = TRUE, recursive = FALSE)
      if (length(old_doublet_files) > 0) unlink(old_doublet_files, recursive = TRUE, force = TRUE)
      rv$doublet_rds_path <- file.path(doublet_dir, "seurat_doublet_filtered.rds")
      rv$doublet_ready <- FALSE
      rv$doublet_plot_paths <- NULL
      rv$doublet_download_specs <- NULL
      save_doublet_plots <- isTRUE(input$savePlotsDoublet)
      selected_doublet_plots <- if (save_doublet_plots) (input$doubletPlotTypes %||% character(0)) else character(0)
      doublet_plot_options <- if (save_doublet_plots) list(
        umap = list(
          save = "umap" %in% selected_doublet_plots,
          filename_prefix = input$doublet_umap_filename,
          title_prefix = input$doublet_umap_title,
          format = input$doublet_umap_format,
          width = input$doublet_umap_width,
          height = input$doublet_umap_height
        ),
        summary = list(
          save = "summary" %in% selected_doublet_plots,
          filename_prefix = input$doublet_summary_filename,
          title_prefix = input$doublet_summary_title,
          format = input$doublet_summary_format,
          width = input$doublet_summary_width,
          height = input$doublet_summary_height
        )
      ) else NULL
      rv$seurat <- remove_doublets(
        seurat_obj = seurat_in,
        output_dir = doublet_dir,
        save_rds = TRUE,
        rds_path = rv$doublet_rds_path,
        doublet_rate = input$doubletRate,
        save_plots = save_doublet_plots,
        plot_options = doublet_plot_options
      )
      rv$doublet_plot_paths <- if (save_doublet_plots) collect_doublet_plot_paths(doublet_dir, doublet_plot_options) else NULL
      rv$doublet_download_specs <- list()
      if (!is.null(rv$doublet_rds_path) && file.exists(rv$doublet_rds_path)) {
        rv$doublet_download_specs[[length(rv$doublet_download_specs) + 1]] <- list(
          id = "downloadDoubletRDS",
          path = rv$doublet_rds_path
        )
      }
      if (!is.null(rv$doublet_plot_paths)) {
        if (!is.null(rv$doublet_plot_paths$summary) && file.exists(rv$doublet_plot_paths$summary)) {
          rv$doublet_download_specs[[length(rv$doublet_download_specs) + 1]] <- list(
            id = "downloadDoubletSummary",
            path = rv$doublet_plot_paths$summary
          )
        }
        if (!is.null(rv$doublet_plot_paths$umap) && length(rv$doublet_plot_paths$umap) > 0) {
          for (i in seq_along(rv$doublet_plot_paths$umap)) {
            rv$doublet_download_specs[[length(rv$doublet_download_specs) + 1]] <- list(
              id = paste0("downloadDoubletUMAP_", i),
              path = rv$doublet_plot_paths$umap[[i]]
            )
          }
        }
      }
      for (spec in rv$doublet_download_specs) {
        local({
          spec_local <- spec
          output[[spec_local$id]] <- downloadHandler(
            filename = function() basename(spec_local$path),
            content = function(file) {
              if (is.null(spec_local$path) || !file.exists(spec_local$path)) {
                stop("Requested doublet output is not available for download.")
              }
              file.copy(spec_local$path, file, overwrite = TRUE)
            }
          )
        })
      }
      rv$doublet_ready <- TRUE
      output$doubletStatus <- renderText("Doublets removed successfully.")
    }, error = function(e) output$doubletStatus <- renderText(paste("Doublet removal failed:", e$message)))
  })
  output$downloadDoubletRDS <- downloadHandler(filename = function() "seurat_doublet_filtered.rds", content = function(file) file.copy(rv$doublet_rds_path, file, overwrite = TRUE))

  observeEvent(input$runIntegration, {
    output$integrationStatus <- renderText("Running integration...")
    tryCatch({
      rv$integration_ready <- FALSE
      rv$integration_plot_paths <- NULL
      rv$integration_download_specs <- NULL
      seurat_input <- get_current_or_uploaded(input$loadRDSIntegration, input$inputRDSIntegration)
      if (is.null(seurat_input) || !inherits(seurat_input, "Seurat")) {
        stop("Integration requires a valid Seurat object. Please load a Seurat .rds file first.")
      }
      save_integration_plots <- isTRUE(input$savePlotsIntegration)
      selected_plots <- if (save_integration_plots) (input$selectedIntegrationPlots %||% character(0)) else character(0)
      plot_options <- NULL
      if (save_integration_plots) {
        plot_options <- list()
        for (plot_id in selected_plots) {
          plot_options[[plot_id]] <- list(
            save = TRUE,
            filename = input[[paste0("integration_", plot_id, "_filename")]],
            format = input[[paste0("integration_", plot_id, "_format")]],
            width = input[[paste0("integration_", plot_id, "_width")]],
            height = input[[paste0("integration_", plot_id, "_height")]],
            title = input[[paste0("integration_", plot_id, "_title")]]
          )
        }
      }
      integration_dir <- file.path(session_tmp, "integration")
      dir.create(integration_dir, recursive = TRUE, showWarnings = FALSE)
      rv$integration_rds_path <- file.path(integration_dir, "integrated_seurat.rds")
      integration_future_max_size_gb <- 25
      integration_future_max_size <- integration_future_max_size_gb * 1024^3
      withProgress(message = "Integrating samples...", value = 0, {
        options(future.globals.maxSize = integration_future_max_size)
        incProgress(0.2, detail = "Running integration")
        integrated <- integrate_seurat_samples(
          seurat_obj = seurat_input,
          output_path = rv$integration_rds_path,
          nfeatures = input$nfeaturesIntegration,
          future_max_size = integration_future_max_size,
          save_rds = TRUE,
          output_dir = integration_dir,
          plot_options = plot_options
        )
        rv$seurat <- integrated
        incProgress(1)
      })
      rv$integration_plot_paths <- if (save_integration_plots) make_integration_plot_paths(integration_dir, plot_options) else NULL
      integration_files <- unlist(rv$integration_plot_paths, use.names = FALSE)
      integration_files <- integration_files[!is.na(integration_files) & nzchar(integration_files)]
      if (save_integration_plots && length(selected_plots) > 0 && length(integration_files) > 0 && !all(file.exists(integration_files))) {
        stop("Integration completed but selected plot files were not generated.")
      }
      rv$integration_download_specs <- list()
      if (!is.null(rv$integration_plot_paths$umap_grouped) && file.exists(rv$integration_plot_paths$umap_grouped)) {
        rv$integration_download_specs[[length(rv$integration_download_specs) + 1]] <- list(
          id = "downloadIntegrationUmapGrouped",
          label = "Download UMAP grouped",
          path = rv$integration_plot_paths$umap_grouped,
          title = "UMAP grouped"
        )
      }
      if (!is.null(rv$integration_plot_paths$umap_split) && file.exists(rv$integration_plot_paths$umap_split)) {
        rv$integration_download_specs[[length(rv$integration_download_specs) + 1]] <- list(
          id = "downloadIntegrationUmapSplit",
          label = "Download UMAP split",
          path = rv$integration_plot_paths$umap_split,
          title = "UMAP split"
        )
      }
      if (!is.null(rv$integration_plot_paths$pie_chart) && file.exists(rv$integration_plot_paths$pie_chart)) {
        rv$integration_download_specs[[length(rv$integration_download_specs) + 1]] <- list(
          id = "downloadIntegrationPieChart",
          label = "Download cell count pie chart",
          path = rv$integration_plot_paths$pie_chart,
          title = "Cell count pie chart"
        )
      }
      rv$integration_ready <- TRUE
      output$integrationStatus <- renderText("Integration completed successfully.")
    }, error = function(e) output$integrationStatus <- renderText(paste("Integration failed:", e$message)))
  })
  output$downloadIntegrationRDS <- downloadHandler(filename = function() "integrated_seurat.rds", content = function(file) file.copy(rv$integration_rds_path, file, overwrite = TRUE))
  output$downloadIntegrationUmapGrouped <- downloadHandler(
    filename = function() {
      path <- rv$integration_plot_paths$umap_grouped
      if (is.null(path)) "integrated_umap_grouped" else basename(path)
    },
    content = function(file) {
      path <- rv$integration_plot_paths$umap_grouped
      if (is.null(path) || !file.exists(path)) stop("Grouped UMAP plot is not available for download.")
      file.copy(path, file, overwrite = TRUE)
    }
  )
  output$downloadIntegrationUmapSplit <- downloadHandler(
    filename = function() {
      path <- rv$integration_plot_paths$umap_split
      if (is.null(path)) "integrated_umap_split" else basename(path)
    },
    content = function(file) {
      path <- rv$integration_plot_paths$umap_split
      if (is.null(path) || !file.exists(path)) stop("Split UMAP plot is not available for download.")
      file.copy(path, file, overwrite = TRUE)
    }
  )
  output$downloadIntegrationPieChart <- downloadHandler(
    filename = function() {
      path <- rv$integration_plot_paths$pie_chart
      if (is.null(path)) "integration_cellcount_piechart" else basename(path)
    },
    content = function(file) {
      path <- rv$integration_plot_paths$pie_chart
      if (is.null(path) || !file.exists(path)) stop("Cell count pie chart is not available for download.")
      file.copy(path, file, overwrite = TRUE)
    }
  )
  output$downloadIntegrationUmapGroupedSidebar <- downloadHandler(
    filename = function() {
      path <- rv$integration_plot_paths$umap_grouped
      if (is.null(path)) "integrated_umap_grouped" else basename(path)
    },
    content = function(file) {
      path <- rv$integration_plot_paths$umap_grouped
      if (is.null(path) || !file.exists(path)) stop("Grouped UMAP plot is not available for download.")
      file.copy(path, file, overwrite = TRUE)
    }
  )
  output$downloadIntegrationUmapSplitSidebar <- downloadHandler(
    filename = function() {
      path <- rv$integration_plot_paths$umap_split
      if (is.null(path)) "integrated_umap_split" else basename(path)
    },
    content = function(file) {
      path <- rv$integration_plot_paths$umap_split
      if (is.null(path) || !file.exists(path)) stop("Split UMAP plot is not available for download.")
      file.copy(path, file, overwrite = TRUE)
    }
  )
  output$downloadIntegrationPieChartSidebar <- downloadHandler(
    filename = function() {
      path <- rv$integration_plot_paths$pie_chart
      if (is.null(path)) "integration_cellcount_piechart" else basename(path)
    },
    content = function(file) {
      path <- rv$integration_plot_paths$pie_chart
      if (is.null(path) || !file.exists(path)) stop("Cell count pie chart is not available for download.")
      file.copy(path, file, overwrite = TRUE)
    }
  )

  observeEvent(input$runDownstream, {
    output$downstreamStatus <- renderText("Running downstream analysis...")
    tryCatch({
      rv$downstream_ready <- FALSE
      rv$downstream_plot_paths <- NULL
      seurat_input <- get_current_or_uploaded(input$loadRDSDownstream, input$inputRDSDownstream)
      dims <- seq(input$dimsUsed[1], input$dimsUsed[2])
      save_downstream_plots <- isTRUE(input$savePlotsDownstream)
      selected_downstream_plots <- if (save_downstream_plots) (input$selectedPlots %||% character(0)) else character(0)
      plot_options <- NULL
      if (save_downstream_plots) {
        plot_options <- list()
        for (plot_id in selected_downstream_plots) {
          plot_options[[plot_id]] <- list(
            save = TRUE,
            filename = input[[paste0(plot_id, "_filename")]],
            format = input[[paste0(plot_id, "_format")]],
            width = input[[paste0(plot_id, "_width")]],
            height = input[[paste0(plot_id, "_height")]],
            title = input[[paste0(plot_id, "_title")]]
          )
        }
      }
      downstream_dir <- file.path(session_tmp, "downstream")
      dir.create(downstream_dir, recursive = TRUE, showWarnings = FALSE)
      rv$downstream_rds_path <- file.path(downstream_dir, "downstream_analyzed.rds")
      withProgress(message = "Running downstream analysis...", value = 0, {
        incProgress(0.3, detail = "Clustering and plotting")
        seurat_input <- run_downstream_analysis(
          seurat_obj = seurat_input,
          dims = dims,
          resolution = input$resolution,
          output_dir = downstream_dir,
          save_rds = TRUE,
          rds_path = rv$downstream_rds_path,
          save_plots = save_downstream_plots,
          plot_options = plot_options
        )
        rv$seurat <- seurat_input
        incProgress(1)
      })
      rv$downstream_plot_paths <- if (save_downstream_plots) make_downstream_plot_paths(downstream_dir, plot_options) else NULL
      downstream_files <- unlist(rv$downstream_plot_paths, use.names = FALSE)
      downstream_files <- downstream_files[!is.na(downstream_files) & nzchar(downstream_files)]
      if (save_downstream_plots && length(downstream_files) > 0 && !all(file.exists(downstream_files))) {
        stop("Downstream analysis did not produce all expected plot files.")
      }
      rv$downstream_ready <- TRUE
      output$downstreamStatus <- renderText("Downstream analysis completed successfully.")
    }, error = function(e) output$downstreamStatus <- renderText(paste("Downstream analysis failed:", e$message)))
  })
  output$downloadDownstreamRDS <- downloadHandler(filename = function() "downstream_analyzed.rds", content = function(file) file.copy(rv$downstream_rds_path, file, overwrite = TRUE))
  output$downloadDownstreamUmapGrouped <- downloadHandler(
    filename = function() {
      path <- rv$downstream_plot_paths$umap_grouped
      if (is.null(path)) "integrated_umap_grouped" else basename(path)
    },
    content = function(file) {
      path <- rv$downstream_plot_paths$umap_grouped
      if (is.null(path) || !file.exists(path)) stop("Grouped UMAP plot is not available for download.")
      file.copy(path, file, overwrite = TRUE)
    }
  )
  output$downloadDownstreamUmapSplit <- downloadHandler(
    filename = function() {
      path <- rv$downstream_plot_paths$umap_split
      if (is.null(path)) "integrated_umap_split" else basename(path)
    },
    content = function(file) {
      path <- rv$downstream_plot_paths$umap_split
      if (is.null(path) || !file.exists(path)) stop("Split UMAP plot is not available for download.")
      file.copy(path, file, overwrite = TRUE)
    }
  )
  output$downloadDownstreamPieChart <- downloadHandler(
    filename = function() {
      path <- rv$downstream_plot_paths$pie_chart
      if (is.null(path)) "cellcount_piechart" else basename(path)
    },
    content = function(file) {
      path <- rv$downstream_plot_paths$pie_chart
      if (is.null(path) || !file.exists(path)) stop("Cell count pie chart is not available for download.")
      file.copy(path, file, overwrite = TRUE)
    }
  )

  observeEvent(input$runDGEA, {
    output$dgeaStatus <- renderText("Running DGEA...")
    tryCatch({
      rv$dgea_ready <- FALSE
      rv$dgea_dir <- NULL
      rv$dgea_csv_all <- NULL
      rv$dgea_csv_top <- NULL
      rv$dgea_plot_paths <- NULL
      rv$dgea_plot_groups <- NULL
      rv$dgea_preview <- NULL
      rv$dgea_download_specs <- NULL
      seurat_input <- get_current_or_uploaded(input$loadRDSDGEA, input$inputRDSDGEA)
      plot_options <- if (!isTRUE(input$savePlotsDGEA)) NULL else list(
        violin = list(
          save = "violin" %in% input$dgeaPlotTypes,
          filename_prefix = input$dgea_violin_filename,
          title_prefix = input$dgea_violin_title,
          format = input$dgea_violin_format,
          width = input$dgea_violin_width,
          height = input$dgea_violin_height
        ),
        dot = list(
          save = "dot" %in% input$dgeaPlotTypes,
          filename_prefix = input$dgea_dot_filename,
          title_prefix = input$dgea_dot_title,
          format = input$dgea_dot_format,
          width = input$dgea_dot_width,
          height = input$dgea_dot_height
        ),
        feature = list(
          save = "feature" %in% input$dgeaPlotTypes,
          filename_prefix = input$dgea_feature_filename,
          title_prefix = input$dgea_feature_title,
          format = input$dgea_feature_format,
          width = input$dgea_feature_width,
          height = input$dgea_feature_height
        )
      )
      dgea_dir <- file.path(session_tmp, "dgea")
      dir.create(dgea_dir, recursive = TRUE, showWarnings = FALSE)
      paths <- make_dgea_paths(dgea_dir)
      withProgress(message = "Running cluster-wise DGEA...", value = 0, {
        incProgress(0.3, detail = "Finding markers")
        run_clusterwise_dgea(seurat_obj = seurat_input, output_dir = dgea_dir, logfc_threshold = input$logfcThreshold, min_pct = input$minPct, plot_options = plot_options)
        incProgress(1)
      })
      rv$dgea_dir <- dgea_dir
      rv$dgea_csv_all <- paths$all_marker
      rv$dgea_csv_top <- paths$top15_marker
      rv$dgea_plot_paths <- if (isTRUE(input$savePlotsDGEA)) {
        plot_dir <- paths$plot_dir
        if (dir.exists(plot_dir)) {
          list.files(plot_dir, full.names = TRUE, recursive = FALSE)
        } else {
          character(0)
        }
      } else {
        character(0)
      }
      rv$dgea_plot_groups <- if (isTRUE(input$savePlotsDGEA)) {
        list(
          violin = dgea_filter_group_paths(rv$dgea_plot_paths, input$dgea_violin_filename, input$dgea_violin_format),
          dot = dgea_filter_group_paths(rv$dgea_plot_paths, input$dgea_dot_filename, input$dgea_dot_format),
          feature = dgea_filter_group_paths(rv$dgea_plot_paths, input$dgea_feature_filename, input$dgea_feature_format)
        )
      } else {
        list(violin = character(0), dot = character(0), feature = character(0))
      }
      rv$dgea_download_specs <- list(
        list(id = "downloadDGEAAllMarkers", type = "table", label = "Download full marker CSV", path = rv$dgea_csv_all, title = "Full marker CSV"),
        list(id = "downloadDGEATopMarkers", type = "table", label = "Download top markers CSV", path = rv$dgea_csv_top, title = "Top markers CSV")
      )
      if (length(rv$dgea_plot_paths) > 0) {
        plot_specs <- lapply(seq_along(rv$dgea_plot_paths), function(i) {
          path <- rv$dgea_plot_paths[[i]]
          list(id = paste0("downloadDGEAPlot", i), type = "plot", label = paste("Download", basename(path)), path = path, title = tools::file_path_sans_ext(basename(path)))
        })
        rv$dgea_download_specs <- c(rv$dgea_download_specs, plot_specs)
      }
      for (spec in rv$dgea_download_specs) {
        local({
          spec_local <- spec
          output[[spec_local$id]] <- downloadHandler(
            filename = function() basename(spec_local$path),
            content = function(file) {
              if (is.null(spec_local$path) || !file.exists(spec_local$path)) {
                stop("Requested DGEA file is not available for download.")
              }
              file.copy(spec_local$path, file, overwrite = TRUE)
            }
          )
        })
      }
      if (!is.null(rv$dgea_csv_all) && file.exists(rv$dgea_csv_all)) {
        rv$dgea_preview <- tryCatch(readr::read_csv(rv$dgea_csv_all, show_col_types = FALSE, n_max = 5), error = function(e) NULL)
      }
      if (is.null(rv$dgea_preview)) {
        rv$dgea_preview <- data.frame()
      }
      if (!all(file.exists(c(rv$dgea_csv_all, rv$dgea_csv_top)))) {
        stop("DGEA did not produce the expected CSV outputs.")
      }
      all_markers_tbl <- readr::read_csv(rv$dgea_csv_all, show_col_types = FALSE)
      top_markers_tbl <- readr::read_csv(rv$dgea_csv_top, show_col_types = FALSE)
      if (nrow(all_markers_tbl) == 0 || nrow(top_markers_tbl) == 0) {
        stop("No marker genes detected from SCT assay.")
      }
      if (isTRUE(input$savePlotsDGEA) && length(input$dgeaPlotTypes) > 0) {
        plot_files <- rv$dgea_plot_paths[file.exists(rv$dgea_plot_paths)]
        if (length(plot_files) == 0) {
          stop("DGEA did not generate any plot outputs.")
        }
      }
      output$dgeaStatus <- renderText("Cluster-wise DGEA completed successfully.")
      rv$dgea_ready <- TRUE
    }, error = function(e) {
      err_msg <- conditionMessage(e)
      if (startsWith(err_msg, "Cluster-wise DGEA failed:") || identical(err_msg, "No marker genes detected from SCT assay.")) {
        output$dgeaStatus <- renderText(err_msg)
      } else {
        output$dgeaStatus <- renderText(paste("DGEA failed:", err_msg))
      }
    })
  })

  observeEvent(input$runAnnotation, {
    output$annotationStatus <- renderText("Running annotation...")
    tryCatch({
      rv$annotation_ready <- FALSE
      rv$annotation_plot_paths <- NULL
      seurat_input <- get_current_or_uploaded(input$loadRDSAnnotate, input$inputRDSAnnotate)
      annotations_input <- NULL
      if (identical(input$annotationInputType, "csv")) {
        req(input$annotationCSV)
        annotations_input <- input$annotationCSV$datapath
      } else {
        annotation_vector <- trimws(unlist(strsplit(input$annotationText, ",", fixed = TRUE)))
        annotation_vector <- annotation_vector[nzchar(annotation_vector)]
        if (length(annotation_vector) == 0) stop("Please provide at least one annotation.")
        annotations_input <- annotation_vector
      }
      save_annotation_plots <- isTRUE(input$savePlotsAnnotate)
      selected_annotation_plots <- if (save_annotation_plots) (input$annotationPlotTypes %||% character(0)) else character(0)
      annotation_plot_options <- if (save_annotation_plots) list(
        umap_annotated = list(
          save = "umap_annotated" %in% selected_annotation_plots,
          filename = input$annotation_umap_annotated_filename,
          title = input$annotation_umap_annotated_title,
          format = input$annotation_umap_annotated_format,
          width = input$annotation_umap_annotated_width,
          height = input$annotation_umap_annotated_height
        ),
        umap_split_by_sample = list(
          save = "umap_split_by_sample" %in% selected_annotation_plots,
          filename = input$annotation_umap_split_by_sample_filename,
          title = input$annotation_umap_split_by_sample_title,
          format = input$annotation_umap_split_by_sample_format,
          width = input$annotation_umap_split_by_sample_width,
          height = input$annotation_umap_split_by_sample_height
        ),
        celltype_piechart = list(
          save = "celltype_piechart" %in% selected_annotation_plots,
          filename = input$annotation_celltype_piechart_filename,
          title = input$annotation_celltype_piechart_title,
          format = input$annotation_celltype_piechart_format,
          width = input$annotation_celltype_piechart_width,
          height = input$annotation_celltype_piechart_height
        )
      ) else NULL
      annotation_dir <- file.path(session_tmp, "annotation")
      dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)
      rv$annotation_rds_path <- file.path(annotation_dir, "annotated_seurat.rds")
      withProgress(message = "Annotating clusters...", value = 0, {
        incProgress(0.5, detail = "Applying annotations")
        annotated_obj <- annotate_clusters(seurat_obj = seurat_input, annotations = annotations_input, output_dir = annotation_dir, rds_path = rv$annotation_rds_path, save_rds = TRUE, save_plots = save_annotation_plots, plot_options = annotation_plot_options)
        rv$seurat <- annotated_obj
        incProgress(1)
      })
      rv$annotation_plot_paths <- if (save_annotation_plots) make_annotation_paths(annotation_dir, annotation_plot_options) else NULL
      rv$annotation_ready <- TRUE
      output$annotationStatus <- renderText("Clusters annotated successfully.")
    }, error = function(e) {
      msg <- conditionMessage(e)
      if (!nzchar(msg)) msg <- paste(class(e), collapse = ", ")
      output$annotationStatus <- renderText(paste("Annotation failed:", msg))
    })
  })
  output$downloadAnnotationRDS <- downloadHandler(filename = function() {
    path <- rv$annotation_rds_path
    if (is.null(path)) "annotated_seurat.rds" else basename(path)
  }, content = function(file) file.copy(rv$annotation_rds_path, file, overwrite = TRUE))
  output$downloadAnnotationUMAP <- downloadHandler(
    filename = function() {
      path <- rv$annotation_plot_paths$umap_annotated
      if (is.null(path)) "umap_annotated.jpeg" else basename(path)
    },
    content = function(file) {
      path <- rv$annotation_plot_paths$umap_annotated
      if (is.null(path) || !file.exists(path)) stop("Annotated UMAP is not available for download.")
      file.copy(path, file, overwrite = TRUE)
    }
  )
  output$downloadAnnotationUMAPSplit <- downloadHandler(
    filename = function() {
      path <- rv$annotation_plot_paths$umap_split_by_sample
      if (is.null(path)) "umap_split_by_sample.jpeg" else basename(path)
    },
    content = function(file) {
      path <- rv$annotation_plot_paths$umap_split_by_sample
      if (is.null(path) || !file.exists(path)) stop("Split annotated UMAP is not available for download.")
      file.copy(path, file, overwrite = TRUE)
    }
  )
  output$downloadAnnotationPieChart <- downloadHandler(
    filename = function() {
      path <- rv$annotation_plot_paths$celltype_piechart
      if (is.null(path)) "celltype_piechart.jpeg" else basename(path)
    },
    content = function(file) {
      path <- rv$annotation_plot_paths$celltype_piechart
      if (is.null(path) || !file.exists(path)) stop("Annotation pie chart is not available for download.")
      file.copy(path, file, overwrite = TRUE)
    }
  )

  observeEvent(input$runFeaturePlot, {
    output$featurePlotStatus <- renderText("Generating feature plots...")
    tryCatch({
      rv$feature_ready <- FALSE
      rv$feature_plot_paths <- NULL
      rv$feature_download_specs <- NULL
      seurat_input <- get_current_or_uploaded(input$loadRDSFeaturePlot, input$inputRDSFeaturePlot)
      genes <- trimws(unlist(strsplit(input$genesToPlot, ",", fixed = TRUE)))
      genes <- genes[nzchar(genes)]
      if (length(genes) == 0) stop("Please enter at least one gene name.")
      feature_plot_options <- if (!isTRUE(input$saveFeaturePlots)) NULL else list(
        filename_prefix = input$featureplot_filename,
        title_prefix = input$featureplot_title,
        format = input$featureplot_format,
        width = input$featureplot_width,
        height = input$featureplot_height
      )
      feature_dir <- file.path(session_tmp, "featureplots")
      dir.create(feature_dir, recursive = TRUE, showWarnings = FALSE)
      withProgress(message = "Plotting genes...", value = 0, {
        incProgress(0.4, detail = "Generating FeaturePlots")
        plot_list <- plot_genes_feature(seurat_obj = seurat_input, genes = genes, output_dir = feature_dir, save_plots = input$saveFeaturePlots, plot_options = feature_plot_options)
        incProgress(1)
      })
      rv$feature_plot_paths <- if (isTRUE(input$saveFeaturePlots)) {
        make_feature_plot_paths(feature_dir, names(plot_list), feature_plot_options)
      } else {
        NULL
      }
      rv$feature_download_specs <- if (isTRUE(input$saveFeaturePlots) && !is.null(rv$feature_plot_paths)) {
        lapply(seq_along(rv$feature_plot_paths), function(i) {
          gene <- names(rv$feature_plot_paths)[[i]]
          safe_gene <- sanitize_file_component(gene)
          list(
            id = paste0("downloadFeaturePlot_", i, "_", safe_gene),
            label = paste("Download", gene),
            path = rv$feature_plot_paths[[gene]],
            title = gene
          )
        })
      } else {
        NULL
      }
      for (spec in rv$feature_download_specs %||% list()) {
        local({
          spec_local <- spec
          output[[spec_local$id]] <- downloadHandler(
            filename = function() basename(spec_local$path),
            content = function(file) {
              if (is.null(spec_local$path) || !file.exists(spec_local$path)) stop("Requested FeaturePlot file is not available for download.")
              file.copy(spec_local$path, file, overwrite = TRUE)
            }
          )
        })
      }
      rv$feature_ready <- TRUE
      output$featurePlotStatus <- renderText("Feature plots generated successfully.")
    }, error = function(e) output$featurePlotStatus <- renderText(paste("Feature plotting failed:", e$message)))
  })
}

shinyApp(ui, server)




