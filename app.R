library(shiny)
library(Seurat)
library(ggplot2)
library(dplyr)
library(shinyBS)

## Load helper functions (server-safe)
app_dir <- normalizePath(".", winslash = "/", mustWork = FALSE)
helper_dir <- file.path(app_dir, "R")
helper_files <- c(
  "read_h5_to_seurat.R",
  "run_qc_metrics.R",
  "remove_doublets.R",
  "post_qc_filtering_and_plots.R",
  "integrate_seurat_samples.R",
  "run_downstream_analysis.R",
  "annotate_clusters.R",
  "plot_genes_feature.R",
  "run_clusterwise_dgea.R"
)
for (f in helper_files) {
  fp <- file.path(helper_dir, f)
  if (file.exists(fp)) {
    source(fp, local = FALSE)
  } else {
    warning(paste("Helper file not found:", fp))
  }
}


seurat_obj <- NULL

options(shiny.maxRequestSize = 10*1024^3)  # 10 GB


ui <- fluidPage(
  tags$head(
    tags$link(rel = "shortcut icon", href = "logo.png"),
    tags$style(HTML("
  html, body {
    background: linear-gradient(135deg, #e6f2f8, #ffffff, #d2ecf7) !important;
    margin: 0;
    padding: 0;
  }

  html::before {
    content: '';
    position: fixed;
    bottom: 20px;
    right: 20px;
    width: 500px;
    height: 500px;
    background-image: url('logo.png');
    background-repeat: no-repeat;
    background-size: contain;
    background-position: bottom right;
    opacity: 0.15;
    z-index: 0;
    pointer-events: none;
  }

  .navbar, .navbar-default {
    background-color: #1c9ed8 !important;
    border-color: #1c9ed8;
  }

  .navbar .navbar-brand,
  .navbar-default .navbar-nav > li > a {
    color: white !important;
  }

  .navbar-default .navbar-nav > li > a:hover,
  .navbar-default .navbar-nav > .active > a {
    background-color: #157bb1 !important;
    color: white !important;
  }

  .tabbable .nav-tabs > li > a {
    color: #1c9ed8 !important;
  }

  .tabbable .nav-tabs > .active > a {
    background-color: #1c9ed8 !important;
    color: white !important;
  }

  .btn-primary {
    background-color: #1c9ed8;
    border-color: #1c9ed8;
  }

  .btn-primary:hover {
    background-color: #157bb1;
    border-color: #157bb1;
  }

  .saves-rds {
    font-weight: 600 !important;
    color: #612020 !important;
    background-color: #fce4e4 !important;
    padding: 6px 10px;
    border-radius: 6px;
    box-shadow: inset 0 0 0 1px #e5bcbc;
  }

  .saves-rds:hover {
    background-color: #c3e2f2 !important;
    color: #102738 !important;
  }

  .no-rds {
    font-weight: 600 !important;
    color: #1c3f5d !important;
    background-color: #d6ecf5 !important;
    padding: 6px 10px;
    border-radius: 6px;
    box-shadow: inset 0 0 0 1px #b2d5e6;
  }

  .no-rds:hover {
    background-color: #f8d4d4 !important;
    color: #3d1515 !important;
  }

  .qc-dl {
    width: 100% !important;
    margin-bottom: 10px;
    font-weight: 600;
    border-radius: 8px;
    box-shadow: 0 2px 6px rgba(0,0,0,0.08);
  }
"))
  ),

  tags$div(
    style = "text-align:center; margin-bottom:10px; margin-top:10px;",
    tags$img(src = "logo.png", height = "120px")
  ),

  navbarPage(NULL, id = "mainNav",
             tabPanel(tags$span("1. Upload Data", class = "no-rds"),
                      sidebarLayout(
                        sidebarPanel(
                          uiOutput("uploadAlertH5"),

                          fileInput("h5file",
                                    label = tagList(
                                      "Upload .h5 File",
                                      tags$span(
                                        icon("question-circle"),
                                        id = "h5fileHelp",
                                        style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                      )
                                    )
                          ),
                          bsTooltip("h5fileHelp",
                                    "Upload a 10X Genomics .h5 formatted file. Ensure it\\'s compatible with Seurat\\'s Read10X_h5().",
                                    placement = "right", trigger = "hover"),

                          textInput("sampleMap",
                                    label = tagList(
                                      "Sample label map (e.g. Sample1=Control,Sample2=BPV)",
                                      tags$span(
                                        icon("question-circle"),
                                        id = "sampleMapHelp",
                                        style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                      )
                                    ),
                                    value = ""
                          ),
                          bsTooltip("sampleMapHelp",
                                    "Comma-separated key-value pairs mapping sample IDs to conditions. Example: \\'Sample1=Control,Sample2=BPV\\'.",
                                    placement = "right", trigger = "hover"),

                          textInput("projectName",
                                    label = tagList(
                                      "Project name",
                                      tags$span(
                                        icon("question-circle"),
                                        id = "projectNameHelp",
                                        style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                      )
                                    ),
                                    value = "SeuratProject"
                          ),
                          bsTooltip("projectNameHelp",
                                    "Project name assigned to the Seurat object. Used for internal labeling and tracking.",
                                    placement = "right", trigger = "hover"),

                          actionButton("readH5", "Read File", class = "btn-primary")
                        ),
                        mainPanel(
                          verbatimTextOutput("uploadStatus")
                        )
                      )
             ),

             tabPanel(
               tags$span("2. Run QC", class = "no-rds"),
               sidebarLayout(
                 sidebarPanel(
                   textInput("mtPattern",
                             label = tagList(
                               "Mitochondrial gene pattern",
                               tags$span(
                                 icon("question-circle"),
                                 id = "mtPatternHelpIcon",
                                 style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                               )
                             ),
                             value = "^MT-"
                   ),
                   bsTooltip("mtPatternHelpIcon",
                             "Regex pattern used to identify mitochondrial genes. Default \\'\\^MT-\\' is for human data. Use \\'\\^mt-\\' for mouse.",
                             placement = "right", trigger = "hover"),

                   textInput("splitBy",
                             label = tagList(
                               "Metadata column to split violin plots",
                               tags$span(
                                 icon("question-circle"),
                                 id = "splitByHelpIcon",
                                 style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                               )
                             ),
                             value = "SampleLabel"
                   ),
                   bsTooltip("splitByHelpIcon",
                             "Splits violin plots by a metadata column (e.g., \\'SampleLabel\\', \\'batch\\') to visualize QC per group.",
                             placement = "right", trigger = "hover"),

                   sliderInput("ylimit",
                               label = tagList(
                                 "Y-limit for percent.mt plot",
                                 tags$span(
                                   icon("question-circle"),
                                   id = "ylimitHelpIcon",
                                   style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                 )
                               ),
                               min = 0, max = 100, value = c(0, 80)
                   ),
                   bsTooltip("ylimitHelpIcon",
                             "Y-axis limit for the mitochondrial content plot (%MT). Helps focus on the desired range.",
                             placement = "right", trigger = "hover"),
                   actionButton("runQC", "Run QC", class = "btn-primary"),

                   tags$hr(),

                   # Show download buttons only after QC is run
                   conditionalPanel(
                     condition = "output.showQCDownloads",
                     tagList(
                       h4("Download QC Plots"),
                       downloadButton("downloadViolin", "Download Violin Plot", class = "btn-primary qc-dl"),
                       downloadButton("downloadPercentMT", "Download Percent MT Plot", class = "btn-primary qc-dl"),
                       downloadButton("downloadScatter", "Download Feature Scatter", class = "btn-primary qc-dl")
                     )
                   )
                 ),

                 mainPanel(
                   uiOutput("qcPlotsDisplay")
                 )
               )
             )


             ,

             tabPanel(tags$span("3. Post-QC Filtering", class = "saves-rds"),
                      sidebarLayout(
                        sidebarPanel(
                          numericInput("minFeat",
                                       label = tagList(
                                         "Minimum features per cell",
                                         tags$span(
                                           icon("question-circle"),
                                           id = "minFeatHelp",
                                           style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                         )
                                       ),
                                       value = 200
                          ),
                          bsTooltip("minFeatHelp",
                                    "Cells with fewer than this number of detected genes (features) will be filtered out. Useful to remove empty droplets or low-quality cells.",
                                    placement = "right", trigger = "hover"),

                          numericInput("maxFeat",
                                       label = tagList(
                                         "Maximum features per cell",
                                         tags$span(
                                           icon("question-circle"),
                                           id = "maxFeatHelp",
                                           style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                         )
                                       ),
                                       value = 7500
                          ),
                          bsTooltip("maxFeatHelp",
                                    "Cells with more than this number of detected genes (features) will be filtered out. Useful for removing doublets or overly complex cells.",
                                    placement = "right", trigger = "hover"),

                          numericInput("maxMT",
                                       label = tagList(
                                         "Max % mitochondrial",
                                         tags$span(
                                           icon("question-circle"),
                                           id = "maxMTHelp",
                                           style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                         )
                                       ),
                                       value = 25
                          ),
                          bsTooltip("maxMTHelp",
                                    "Cells with mitochondrial gene content above this threshold will be removed. High mitochondrial content often suggests low-quality or dying cells.",
                                    placement = "right", trigger = "hover"),

                          textInput("splitByPost",
                                    label = tagList(
                                      "Split plots by",
                                      tags$span(
                                        icon("question-circle"),
                                        id = "splitByPostHelp",
                                        style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                      )
                                    ),
                                    value = "SampleLabel"
                          ),
                          bsTooltip("splitByPostHelp",
                                    "Metadata column used to split visualizations by sample or group (e.g. SampleLabel, batch).",
                                    placement = "right", trigger = "hover"),

                          actionButton("postQC", "Run Post-QC Filtering", class = "btn-primary"),

                          tags$hr(),

                          # Show download buttons only after Post-QC is run
                          conditionalPanel(
                            condition = "output.showPostQCDownloads",
                            tagList(
                              h4("Download Post-QC Outputs"),
                              downloadButton("downloadPostQCRDS", "Download Filtered Seurat (.rds)", class = "btn-primary qc-dl"),
                              downloadButton("downloadPostQCVln", "Download Post-QC Violin Plot", class = "btn-primary qc-dl"),
                              downloadButton("downloadPostQCDensity", "Download Post-QC Density Plot", class = "btn-primary qc-dl"),
                              downloadButton("downloadPostQCPie", "Download Post-QC Cellcount Pie", class = "btn-primary qc-dl")
                            )
                          )
                        ),
                        mainPanel(
                          verbatimTextOutput("postQCStatus"),
                          uiOutput("postQCPlotsDisplay")
                        )
                      )
             ),

tabPanel(tags$span("4. Remove Doublets", class = "saves-rds"),
                      sidebarLayout(
                        sidebarPanel(
                          uiOutput("uploadAlertRDS"),

                          checkboxInput("loadRDSInput",
                                        label = tagList(
                                          "Load Seurat object from .rds file",
                                          tags$span(
                                            icon("question-circle"),
                                            id = "loadRDSInputHelp",
                                            style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                          )
                                        ),
                                        value = FALSE
                          ),
                          bsTooltip("loadRDSInputHelp",
                                    "Check this box if you want to load a previously saved Seurat object from an .rds file.",
                                    placement = "right", trigger = "hover"),

                          conditionalPanel(
                            condition = "input.loadRDSInput == true",
                            fileInput("inputRDSfile",
                                      label = tagList(
                                        "Select .rds file",
                                        tags$span(
                                          icon("question-circle"),
                                          id = "inputRDSfileHelp",
                                          style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                        )
                                      )
                            )
                          ),
                          bsTooltip("inputRDSfileHelp",
                                    "Upload an existing .rds file containing a Seurat object for doublet removal.",
                                    placement = "right", trigger = "hover"),

                          textInput("doubletOutDir",
                                    label = tagList(
                                      "Output directory",
                                      tags$span(
                                        icon("question-circle"),
                                        id = "doubletOutDirHelp",
                                        style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                      )
                                    ),
                                    value = "./Plots/DoubletFinder"
                          ),
                          bsTooltip("doubletOutDirHelp",
                                    "Directory where doublet diagnostic plots will be saved.",
                                    placement = "right", trigger = "hover"),

                          checkboxInput("saveRDSDoublet",
                                        label = tagList(
                                          "Save filtered Seurat object?",
                                          tags$span(
                                            icon("question-circle"),
                                            id = "saveRDSDoubletHelp",
                                            style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                          )
                                        ),
                                        value = TRUE
                          ),
                          bsTooltip("saveRDSDoubletHelp",
                                    "Check to save the Seurat object after doublets have been removed.",
                                    placement = "right", trigger = "hover"),

                          textInput("rdsPathDoublet",
                                    label = tagList(
                                      "RDS path",
                                      tags$span(
                                        icon("question-circle"),
                                        id = "rdsPathDoubletHelp",
                                        style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                      )
                                    ),
                                    value = "seurat.obj_doubletfiltered.rds"
                          ),
                          bsTooltip("rdsPathDoubletHelp",
                                    "File path where the filtered Seurat object (.rds) will be saved.",
                                    placement = "right", trigger = "hover"),

                          numericInput("doubletRate",
                                       label = tagList(
                                         "Expected doublet rate",
                                         tags$span(
                                           icon("question-circle"),
                                           id = "doubletRateHelp",
                                           style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                         )
                                       ),
                                       value = 0.075
                          ),
                          bsTooltip("doubletRateHelp",
                                    "Expected proportion of doublets in your dataset. This is typically 5-10% for droplet-based protocols.",
                                    placement = "right", trigger = "hover"),

                          actionButton("removeDoublets", "Run Doublet Removal", class = "btn-primary")
                        ),
                        mainPanel(
                          verbatimTextOutput("doubletStatus")
                        )
                      )
             )
             ,

             tabPanel(tags$span("5. Integration", class = "saves-rds"),
                      sidebarLayout(
                        sidebarPanel(
                          uiOutput("uploadAlertRDS"),

                          checkboxInput("loadRDSIntegration",
                                        label = tagList(
                                          "Load Seurat object from .rds file",
                                          tags$span(
                                            icon("question-circle"),
                                            id = "loadRDSIntegrationHelp",
                                            style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                          )
                                        ),
                                        value = FALSE
                          ),
                          bsTooltip("loadRDSIntegrationHelp",
                                    "Check this box if you want to load a previously saved Seurat object from an .rds file for integration.",
                                    placement = "right", trigger = "hover"),

                          conditionalPanel(
                            condition = "input.loadRDSIntegration == true",
                            fileInput("inputRDSIntegration",
                                      label = tagList(
                                        "Select .rds file",
                                        tags$span(
                                          icon("question-circle"),
                                          id = "inputRDSIntegrationHelp",
                                          style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                        )
                                      )
                            )
                          ),
                          bsTooltip("inputRDSIntegrationHelp",
                                    "Upload an .rds file containing a Seurat object to be used in integration.",
                                    placement = "right", trigger = "hover"),

                          numericInput("nfeaturesIntegration",
                                       label = tagList(
                                         "Number of integration features",
                                         tags$span(
                                           icon("question-circle"),
                                           id = "nfeaturesIntegrationHelp",
                                           style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                         )
                                       ),
                                       value = 3000, min = 500, step = 500
                          ),
                          bsTooltip("nfeaturesIntegrationHelp",
                                    "Number of top variable genes to use for integration. Typical values are between 2000 and 3000.",
                                    placement = "right", trigger = "hover"),

                          numericInput("futureMaxSize",
                                       label = tagList(
                                         "Max future memory (GB)",
                                         tags$span(
                                           icon("question-circle"),
                                           id = "futureMaxSizeHelp",
                                           style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                         )
                                       ),
                                       value = 2
                          ),
                          bsTooltip("futureMaxSizeHelp",
                                    "Maximum allowed memory (in GB) for parallel Seurat integration. Increase if you encounter memory issues.",
                                    placement = "right", trigger = "hover"),

                          checkboxInput("saveRDSIntegration",
                                        label = tagList(
                                          "Save integrated Seurat object?",
                                          tags$span(
                                            icon("question-circle"),
                                            id = "saveRDSIntegrationHelp",
                                            style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                          )
                                        ),
                                        value = TRUE
                          ),
                          bsTooltip("saveRDSIntegrationHelp",
                                    "Check to save the integrated Seurat object as an .rds file for later use.",
                                    placement = "right", trigger = "hover"),

                          textInput("rdsPathIntegration",
                                    label = tagList(
                                      "Output RDS path",
                                      tags$span(
                                        icon("question-circle"),
                                        id = "rdsPathIntegrationHelp",
                                        style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                      )
                                    ),
                                    value = "RDS_Files/integrated_seurat.obj.rds"
                          ),
                          bsTooltip("rdsPathIntegrationHelp",
                                    "File path where the integrated Seurat object will be saved (e.g., \\'RDS_Files/integrated_seurat.obj.rds\\').",
                                    placement = "right", trigger = "hover"),

                          actionButton("runIntegration", "Run Integration", class = "btn-primary")
                        ),
                        mainPanel(
                          verbatimTextOutput("integrationStatus")
                        )
                      )
             ),

             tabPanel(tags$span("6. Downstream Analysis", class = "saves-rds"),
                      sidebarLayout(
                        sidebarPanel(
                          uiOutput("uploadAlertRDS"),

                          checkboxInput("loadRDSDownstream",
                                        label = tagList(
                                          "Load Seurat object from .rds file",
                                          tags$span(
                                            icon("question-circle"),
                                            id = "loadRDSDownstreamHelp",
                                            style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                          )
                                        ),
                                        value = FALSE
                          ),
                          bsTooltip("loadRDSDownstreamHelp",
                                    "Check this box if you want to load a previously saved Seurat object from an .rds file for downstream analysis.",
                                    placement = "right", trigger = "hover"),

                          conditionalPanel(
                            condition = "input.loadRDSDownstream == true",
                            fileInput("inputRDSDownstream",
                                      label = tagList(
                                        "Select .rds file",
                                        tags$span(
                                          icon("question-circle"),
                                          id = "inputRDSDownstreamHelp",
                                          style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                        )
                                      )
                            )
                          ),
                          bsTooltip("inputRDSDownstreamHelp",
                                    "Upload a .rds file containing a Seurat object to be used for downstream steps.",
                                    placement = "right", trigger = "hover"),

                          sliderInput("dimsUsed",
                                      label = tagList(
                                        "PCA/UMAP dimensions",
                                        tags$span(
                                          icon("question-circle"),
                                          id = "dimsUsedHelp",
                                          style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                        )
                                      ),
                                      min = 1, max = 50, value = c(1, 15)
                          ),
                          bsTooltip("dimsUsedHelp",
                                    "Range of principal components (PCs) to use for UMAP, clustering, and downstream analysis.",
                                    placement = "right", trigger = "hover"),

                          numericInput("resolution",
                                       label = tagList(
                                         "Clustering resolution",
                                         tags$span(
                                           icon("question-circle"),
                                           id = "resolutionHelp",
                                           style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                         )
                                       ),
                                       value = 0.1, step = 0.1
                          ),
                          bsTooltip("resolutionHelp",
                                    "Clustering resolution parameter. Higher values lead to more clusters. Try 0.1–1.0 for exploration.",
                                    placement = "right", trigger = "hover"),

                          textInput("downstreamOutDir",
                                    label = tagList(
                                      "Output directory for plots",
                                      tags$span(
                                        icon("question-circle"),
                                        id = "downstreamOutDirHelp",
                                        style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                      )
                                    ),
                                    value = "./Plots"
                          ),
                          bsTooltip("downstreamOutDirHelp",
                                    "Directory where UMAP plots, pie charts, and other figures will be saved.",
                                    placement = "right", trigger = "hover"),

                          checkboxInput("saveRDSDownstream",
                                        label = tagList(
                                          "Save downstream Seurat object?",
                                          tags$span(
                                            icon("question-circle"),
                                            id = "saveRDSDownstreamHelp",
                                            style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                          )
                                        ),
                                        value = TRUE
                          ),
                          bsTooltip("saveRDSDownstreamHelp",
                                    "Check to save the Seurat object after downstream steps are completed.",
                                    placement = "right", trigger = "hover"),

                          textInput("rdsPathDownstream",
                                    label = tagList(
                                      "RDS path",
                                      tags$span(
                                        icon("question-circle"),
                                        id = "rdsPathDownstreamHelp",
                                        style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                      )
                                    ),
                                    value = "RDS_Files/downstream_analyzed.rds"
                          ),
                          bsTooltip("rdsPathDownstreamHelp",
                                    "File path where the downstream-analyzed Seurat object will be saved (e.g., \\'RDS_Files/downstream_analyzed.rds\\').",
                                    placement = "right", trigger = "hover"),

                          tags$hr(),
                          h4("Customize Plots to Save"),
                          checkboxGroupInput("selectedPlots", "Select plots to save:",
                                             choices = list(
                                               "UMAP grouped by samples" = "umap_grouped",
                                               "UMAP split by samples" = "umap_split",
                                               "Cell count pie chart" = "pie_chart"
                                             )
                          ),

                          uiOutput("plotSettingsUI"),

                          actionButton("runDownstream", "Run Downstream Analysis", class = "btn-primary")
                        ),
                        mainPanel(
                          verbatimTextOutput("downstreamStatus")
                        )
                      )
             ),

             tabPanel(tags$span("7. Cluster-wise DGEA", class = "no-rds"),
                      sidebarLayout(
                        sidebarPanel(
                          uiOutput("uploadAlertRDS"),

                          checkboxInput("loadRDSDGEA",
                                        label = tagList(
                                          "Load Seurat object from .rds file",
                                          tags$span(
                                            icon("question-circle"),
                                            id = "loadRDSDGEAHelp",
                                            style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                          )
                                        ),
                                        value = FALSE
                          ),
                          bsTooltip("loadRDSDGEAHelp",
                                    "Check this box to load a previously saved Seurat object from an .rds file for DGEA.",
                                    placement = "right", trigger = "hover"),

                          conditionalPanel(
                            condition = "input.loadRDSDGEA == true",
                            fileInput("inputRDSDGEA",
                                      label = tagList(
                                        "Select .rds file",
                                        tags$span(
                                          icon("question-circle"),
                                          id = "inputRDSDGEAHelp",
                                          style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                        )
                                      )
                            )
                          ),
                          bsTooltip("inputRDSDGEAHelp",
                                    "Upload a .rds file containing the Seurat object to use for differential gene expression analysis.",
                                    placement = "right", trigger = "hover"),

                          numericInput("logfcThreshold",
                                       label = tagList(
                                         "LogFC threshold",
                                         tags$span(
                                           icon("question-circle"),
                                           id = "logfcThresholdHelp",
                                           style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                         )
                                       ),
                                       value = 1
                          ),
                          bsTooltip("logfcThresholdHelp",
                                    "Minimum log fold-change for a gene to be considered differentially expressed. Default is 1.",
                                    placement = "right", trigger = "hover"),

                          numericInput("minPct",
                                       label = tagList(
                                         "Minimum percent expression",
                                         tags$span(
                                           icon("question-circle"),
                                           id = "minPctHelp",
                                           style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                         )
                                       ),
                                       value = 0.25
                          ),
                          bsTooltip("minPctHelp",
                                    "Minimum fraction of cells in either cluster expressing the gene. Used to filter low-expressed genes.",
                                    placement = "right", trigger = "hover"),

                          textInput("dgeaOutputDir",
                                    label = tagList(
                                      "Output directory for DGEA results",
                                      tags$span(
                                        icon("question-circle"),
                                        id = "dgeaOutputDirHelp",
                                        style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                      )
                                    ),
                                    value = "./DGEA"
                          ),
                          bsTooltip("dgeaOutputDirHelp",
                                    "Directory where DGEA results and plots will be saved.",
                                    placement = "right", trigger = "hover"),

                          checkboxInput("savePlotsDGEA",
                                        label = tagList(
                                          "Save plots?",
                                          tags$span(
                                            icon("question-circle"),
                                            id = "savePlotsDGEAHelp",
                                            style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                          )
                                        ),
                                        value = TRUE
                          ),
                          conditionalPanel(
                            condition = "input.savePlotsDGEA == true",
                            checkboxGroupInput("dgeaPlotTypes",
                                               label = tagList(
                                                 "Choose plot types to save",
                                                 tags$span(
                                                   icon("question-circle"),
                                                   id = "dgeaPlotTypesHelp",
                                                   style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                                 )
                                               ),
                                               choices = list(
                                                 "Violin plot" = "violin",
                                                 "Dot plot" = "dot",
                                                 "One FeaturePlot per gene" = "feature"
                                               ),
                                               selected = c("violin", "dot", "feature")
                            ),
                            bsTooltip("dgeaPlotTypesHelp",
                                      "Select which types of plots to generate for each cluster\\'s top marker genes. FeaturePlots show gene expression on UMAP; one is generated per gene.",
                                      placement = "right", trigger = "hover"),

                            radioButtons("dgeaPlotFormat",
                                         label = tagList(
                                           "Image format",
                                           tags$span(
                                             icon("question-circle"),
                                             id = "dgeaPlotFormatHelp",
                                             style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                           )
                                         ),
                                         choices = c("jpeg", "png", "svg"), inline = TRUE),
                            bsTooltip("dgeaPlotFormatHelp",
                                      "Choose the image format for exported plots. SVG is best for vector graphics; JPEG is widely supported.",
                                      placement = "right", trigger = "hover"),

                            textInput("dgeaTitlePrefix",
                                      label = tagList(
                                        "Plot title prefix",
                                        tags$span(
                                          icon("question-circle"),
                                          id = "dgeaTitlePrefixHelp",
                                          style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                        )
                                      ),
                                      value = "Cluster"),
                            bsTooltip("dgeaTitlePrefixHelp",
                                      "Optional prefix added to all plot titles. For example, if prefix = \\'MySample\\', the title becomes \\'MySample 3 - CD3D\\'.",
                                      placement = "right", trigger = "hover"),

                            fluidRow(
                              column(6,
                                     numericInput("dgeaPlotWidth",
                                                  label = tagList(
                                                    "Width (in)",
                                                    tags$span(
                                                      icon("question-circle"),
                                                      id = "dgeaPlotWidthHelp",
                                                      style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                                    )
                                                  ),
                                                  value = 10, min = 4, max = 30)
                              ),
                              column(6,
                                     numericInput("dgeaPlotHeight",
                                                  label = tagList(
                                                    "Height (in)",
                                                    tags$span(
                                                      icon("question-circle"),
                                                      id = "dgeaPlotHeightHelp",
                                                      style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                                    )
                                                  ),
                                                  value = 8, min = 4, max = 30)
                              )
                            ),
                            bsTooltip("dgeaPlotWidthHelp",
                                      "Controls the horizontal size of each plot file in inches.",
                                      placement = "right", trigger = "hover"),
                            bsTooltip("dgeaPlotHeightHelp",
                                      "Controls the vertical size of each plot file in inches.",
                                      placement = "right", trigger = "hover")
                          ),
                          bsTooltip("savePlotsDGEAHelp",
                                    "If checked, violin plots, dot plots, and feature plots for top markers will be saved.",
                                    placement = "right", trigger = "hover"),

                          actionButton("runDGEA", "Run Cluster-wise DGEA", class = "btn-primary")
                        ),
                        mainPanel(
                          verbatimTextOutput("dgeaStatus")
                        )
                      )
             ),

             tabPanel(tags$span("8. Annotate Clusters", class = "saves-rds"),
                      sidebarLayout(
                        sidebarPanel(
                          uiOutput("uploadAlertRDS"),

                          checkboxInput("loadRDSAnnotate",
                                        label = tagList(
                                          "Load Seurat object from .rds file",
                                          tags$span(
                                            icon("question-circle"),
                                            id = "loadRDSAnnotateHelp",
                                            style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                          )
                                        ),
                                        value = FALSE
                          ),
                          bsTooltip("loadRDSAnnotateHelp",
                                    "Check to load an existing annotated or clustered Seurat object from an .rds file.",
                                    placement = "right", trigger = "hover"),

                          conditionalPanel(
                            condition = "input.loadRDSAnnotate == true",
                            fileInput("inputRDSAnnotate",
                                      label = tagList(
                                        "Select Seurat .rds file",
                                        tags$span(
                                          icon("question-circle"),
                                          id = "inputRDSAnnotateHelp",
                                          style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                        )
                                      )
                            )
                          ),
                          bsTooltip("inputRDSAnnotateHelp",
                                    "Upload the .rds file containing the Seurat object to annotate clusters in.",
                                    placement = "right", trigger = "hover"),

                          radioButtons("annotationInputType",
                                       label = tagList(
                                         "Annotation input type",
                                         tags$span(
                                           icon("question-circle"),
                                           id = "annotationInputTypeHelp",
                                           style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                         )
                                       ),
                                       choices = c("Upload CSV" = "csv", "Enter manually" = "manual"),
                                       selected = "csv"
                          ),
                          bsTooltip("annotationInputTypeHelp",
                                    "Choose how to provide cluster annotations: upload a CSV file or type them manually.",
                                    placement = "right", trigger = "hover"),

                          conditionalPanel(
                            condition = "input.annotationInputType == 'csv'",
                            fileInput("annotationCSV",
                                      label = tagList(
                                        "Upload CSV file for annotations (2 columns: seurat_clusters, cell_annotation)",
                                        tags$span(
                                          icon("question-circle"),
                                          id = "annotationCSVHelp",
                                          style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                        )
                                      )
                            )
                          ),
                          bsTooltip("annotationCSVHelp",
                                    "Upload a CSV with two columns: \\\"seurat_clusters\\\" and \\\"cell_annotation\\\" to map clusters to labels.",
                                    placement = "right", trigger = "hover"),

                          conditionalPanel(
                            condition = "input.annotationInputType == 'manual'",
                            textAreaInput("annotationText",
                                          label = tagList(
                                            "Enter annotations (comma-separated, order matters)",
                                            tags$span(
                                              icon("question-circle"),
                                              id = "annotationTextHelp",
                                              style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                            )
                                          ),
                                          value = "", rows = 5
                            )
                          ),
                          bsTooltip("annotationTextHelp",
                                    "Enter annotations manually as comma-separated values. The order must match the cluster IDs.",
                                    placement = "right", trigger = "hover"),

                          textInput("annotOutputDir",
                                    label = tagList(
                                      "Output directory for annotation plots",
                                      tags$span(
                                        icon("question-circle"),
                                        id = "annotOutputDirHelp",
                                        style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                      )
                                    ),
                                    value = "./Plots/Annotation"
                          ),
                          bsTooltip("annotOutputDirHelp",
                                    "Directory where annotation-based plots (e.g., UMAP with new labels) will be saved.",
                                    placement = "right", trigger = "hover"),

                          checkboxInput("saveRDSAnnotate",
                                        label = tagList(
                                          "Save annotated Seurat object?",
                                          tags$span(
                                            icon("question-circle"),
                                            id = "saveRDSAnnotateHelp",
                                            style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                          )
                                        ),
                                        value = TRUE
                          ),
                          bsTooltip("saveRDSAnnotateHelp",
                                    "Check to save the annotated Seurat object as an .rds file.",
                                    placement = "right", trigger = "hover"),

                          textInput("rdsPathAnnotate",
                                    label = tagList(
                                      "Path to save annotated .rds",
                                      tags$span(
                                        icon("question-circle"),
                                        id = "rdsPathAnnotateHelp",
                                        style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                      )
                                    ),
                                    value = "RDS_Files/annotated_seurat.obj.rds"
                          ),
                          bsTooltip("rdsPathAnnotateHelp",
                                    "File path where the annotated Seurat object will be saved (e.g., \\'RDS_Files/annotated_seurat.obj.rds\\').",
                                    placement = "right", trigger = "hover"),

                          actionButton("runAnnotation", "Annotate Clusters")
                        ),
                        mainPanel(
                          verbatimTextOutput("annotationStatus")
                        )
                      )
             ),

             tabPanel(tags$span("9. FeaturePlot Genes", class = "no-rds"),
                      sidebarLayout(
                        sidebarPanel(
                          uiOutput("uploadAlertRDS"),

                          checkboxInput("loadRDSFeaturePlot",
                                        label = tagList(
                                          "Load Seurat object from .rds file",
                                          tags$span(
                                            icon("question-circle"),
                                            id = "loadRDSFeaturePlotHelp",
                                            style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                          )
                                        ),
                                        value = FALSE
                          ),
                          bsTooltip("loadRDSFeaturePlotHelp",
                                    "Check this box to load a previously saved Seurat object (.rds) for FeaturePlot visualization.",
                                    placement = "right", trigger = "hover"),

                          conditionalPanel(
                            condition = "input.loadRDSFeaturePlot == true",
                            fileInput("inputRDSFeaturePlot",
                                      label = tagList(
                                        "Select Seurat .rds file",
                                        tags$span(
                                          icon("question-circle"),
                                          id = "inputRDSFeaturePlotHelp",
                                          style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                        )
                                      )
                            )
                          ),
                          bsTooltip("inputRDSFeaturePlotHelp",
                                    "Upload a Seurat .rds file to use for plotting gene expression across the UMAP embedding.",
                                    placement = "right", trigger = "hover"),

                          textAreaInput("genesToPlot",
                                        label = tagList(
                                          "Enter genes to plot (comma-separated)",
                                          tags$span(
                                            icon("question-circle"),
                                            id = "genesToPlotHelp",
                                            style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                          )
                                        ),
                                        value = "", rows = 5
                          ),
                          bsTooltip("genesToPlotHelp",
                                    "Enter a comma-separated list of gene symbols (e.g., CD3D, MS4A1, LYZ) to visualize with FeaturePlot.",
                                    placement = "right", trigger = "hover"),

                          textInput("featurePlotOutputDir",
                                    label = tagList(
                                      "Output directory for FeaturePlots",
                                      tags$span(
                                        icon("question-circle"),
                                        id = "featurePlotOutputDirHelp",
                                        style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                      )
                                    ),
                                    value = "./Plots/FeaturePlots"
                          ),
                          bsTooltip("featurePlotOutputDirHelp",
                                    "Directory where FeaturePlot images will be saved.",
                                    placement = "right", trigger = "hover"),

                          checkboxInput("saveFeaturePlots",
                                        label = tagList(
                                          "Save FeaturePlots as images?",
                                          tags$span(
                                            icon("question-circle"),
                                            id = "saveFeaturePlotsHelp",
                                            style = "margin-left: 5px; cursor: pointer; color: #1c9ed8;"
                                          )
                                        ),
                                        value = TRUE
                          ),
                          bsTooltip("saveFeaturePlotsHelp",
                                    "Check to save FeaturePlot outputs as image files (e.g., JPEG) in the specified directory.",
                                    placement = "right", trigger = "hover"),

                          actionButton("runFeaturePlot", "Plot Genes")
                        ),

                        mainPanel(
                          verbatimTextOutput("featurePlotStatus"),
                          uiOutput("featurePlotOutputs")
                        )
                      )
             )



  ),

  # RAM usage bar
  # uiOutput("ramBarUI")
)

server <- function(input, output, session) {

  rv <- reactiveValues(
    seurat = NULL,
    qc_plot_paths = NULL,
    qc_ready = FALSE,
    qc_counter = 0  # 🔁 counter to force re-render
  ,
    postqc_plot_paths = NULL,
    postqc_rds_path = NULL,
    postqc_ready = FALSE,
    postqc_counter = 0
  )


  output$plotSettingsUI <- renderUI({
    req(input$selectedPlots)

    plot_ui <- lapply(input$selectedPlots, function(plot_id) {
      default_name <- switch(plot_id,
                             "umap_grouped" = "integrated_umap_grouped",
                             "umap_split"   = "integrated_umap_split",
                             "pie_chart"    = "cellcount_piechart")

      default_title <- switch(plot_id,
                              "umap_grouped" = "UMAP: Grouped by Samples",
                              "umap_split"   = "UMAP: Split by Samples",
                              "pie_chart"    = "Cell Count by Sample")
      wellPanel(
        h5(strong(switch(plot_id,
                         "umap_grouped" = "UMAP Grouped Plot",
                         "umap_split"   = "UMAP Split Plot",
                         "pie_chart"    = "Cell Count Pie Chart"
        ))),
        textInput(paste0(plot_id, "_filename"), "Filename:", value = default_name),
        textInput(paste0(plot_id, "_title"), "Plot title (optional):", value = default_title),
        radioButtons(paste0(plot_id, "_format"), "Image format:", choices = c("jpeg", "png", "svg"), inline = TRUE),
        fluidRow(
          column(6, numericInput(paste0(plot_id, "_width"), "Width (in)", value = 10, min = 4, max = 30)),
          column(6, numericInput(paste0(plot_id, "_height"), "Height (in)", value = 8, min = 4, max = 30))
        )
      )
    })

    do.call(tagList, plot_ui)
  })

  output$qcPlotsDisplay <- renderUI({
    req(rv$qc_ready)
    rv$qc_counter  # force re-render when incremented

    encode_img <- function(path) {
      if (!file.exists(path)) return(NULL)
      # Read file raw and re-encode every time
      mime <- if (grepl("\\.png$", path)) "image/png" else "image/jpeg"
      encoded <- base64enc::dataURI(file = path, mime = mime)
      encoded
    }

    tagList(
      h4("QC Plots"),
      img(src = encode_img(rv$qc_plot_paths$violin), height = "400px", style = "margin:10px;"),
      img(src = encode_img(rv$qc_plot_paths$percent_mt), height = "400px", style = "margin:10px;"),
      img(src = encode_img(rv$qc_plot_paths$scatter), height = "400px", style = "margin:10px;")
    )
  })




  output$downloadViolin <- downloadHandler(
  filename = function() { "pre_qc_violin.jpg" },
  content = function(file) {
    req(rv$qc_plot_paths$violin)
    file.copy(rv$qc_plot_paths$violin, file)
  }
)

output$downloadPercentMT <- downloadHandler(
  filename = function() { "pre_qc_percent.mt.jpg" },
  content = function(file) {
    req(rv$qc_plot_paths$percent_mt)
    file.copy(rv$qc_plot_paths$percent_mt, file)
  }
)

output$downloadScatter <- downloadHandler(
  filename = function() { "pre_qc_FeatureScatterplot.jpg" },
  content = function(file) {
    req(rv$qc_plot_paths$scatter)
    file.copy(rv$qc_plot_paths$scatter, file)
  }


  # ---------------------------
  # Post-QC: display + downloads
  # ---------------------------

  output$postQCPlotsDisplay <- renderUI({
    req(rv$postqc_ready)
    rv$postqc_counter  # force re-render when incremented

    encode_img <- function(path) {
      if (!file.exists(path)) return(NULL)
      mime <- if (grepl("\\.png$", path)) "image/png" else "image/jpeg"
      base64enc::dataURI(file = path, mime = mime)
    }

    tagList(
      h4("Post-QC Plots"),
      img(src = encode_img(rv$postqc_plot_paths$vln), height = "400px", style = "margin:10px;"),
      img(src = encode_img(rv$postqc_plot_paths$density), height = "400px", style = "margin:10px;"),
      img(src = encode_img(rv$postqc_plot_paths$pie), height = "400px", style = "margin:10px;")
    )
  })

  output$downloadPostQCRDS <- downloadHandler(
    filename = function() { "seurat_post_qc_filtered.rds" },
    content = function(file) {
      req(rv$postqc_rds_path)
      file.copy(rv$postqc_rds_path, file, overwrite = TRUE)
    }
  )

  output$downloadPostQCVln <- downloadHandler(
    filename = function() { "post_qc_vplot.jpg" },
    content = function(file) {
      req(rv$postqc_plot_paths$vln)
      file.copy(rv$postqc_plot_paths$vln, file, overwrite = TRUE)
    }
  )

  output$downloadPostQCDensity <- downloadHandler(
    filename = function() { "post_qc_densityplot.jpg" },
    content = function(file) {
      req(rv$postqc_plot_paths$density)
      file.copy(rv$postqc_plot_paths$density, file, overwrite = TRUE)
    }
  )

  output$downloadPostQCPie <- downloadHandler(
    filename = function() { "cellcount_piechart_post.jpg" },
    content = function(file) {
      req(rv$postqc_plot_paths$pie)
      file.copy(rv$postqc_plot_paths$pie, file, overwrite = TRUE)
    }
  )

)


output$showQCDownloads <- reactive({
  isTRUE(rv$qc_ready)
})
outputOptions(output, "showQCDownloads", suspendWhenHidden = FALSE)


output$showPostQCDownloads <- reactive({
  isTRUE(rv$postqc_ready)
})
outputOptions(output, "showPostQCDownloads", suspendWhenHidden = FALSE)



  # get_total_ram_mb <- function() {
  #   sysname <- Sys.info()[["sysname"]]
  #
  #   if (sysname == "Windows") {
  #     # Use WMIC to get total memory in bytes
  #     raw_output <- system("wmic computersystem get TotalPhysicalMemory", intern = TRUE)
  #     mem_line <- raw_output[2]
  #     mem_bytes <- as.numeric(gsub("[^0-9]", "", mem_line))
  #     total_ram <- mem_bytes / 1024^2  # convert to MB
  #   } else if (file.exists("/proc/meminfo")) {
  #     # Linux
  #     meminfo <- readLines("/proc/meminfo")
  #     mem_total_line <- grep("^MemTotal:", meminfo, value = TRUE)
  #     mem_kb <- as.numeric(gsub("[^0-9]", "", mem_total_line))
  #     total_ram <- mem_kb / 1024  # convert KB to MB
  #   } else if (sysname == "Darwin") {
  #     # macOS
  #     mem_bytes <- as.numeric(system("sysctl -n hw.memsize", intern = TRUE))
  #     total_ram <- mem_bytes / 1024^2  # bytes to MB
  #   } else {
  #     warning("OS not recognized; defaulting to 16 GB")
  #     total_ram <- 16000
  #   }
  #
  #   return(round(total_ram, 0))
  # }
  #
  # total_ram_mb <- get_total_ram_mb()
  #
  #
  # autoInvalidate <- reactiveTimer(3000)  # update every 3 seconds

  output$uploadAlertH5 <- renderUI({
    if (is.null(input$h5file)) return(NULL)

    tags$div(
      style = "color: #856404; background-color: #fff3cd; border: 1px solid #ffeeba; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
      tags$strong("\u26A0\ufe0f Please wait: "),
      "Wait for data upload to complete before proceeding."
    )
  })

  output$uploadAlertRDS <- renderUI({
    any_uploading <- FALSE

    if (isTruthy(input$loadRDSInput) && isTruthy(input$inputRDSfile)) any_uploading <- TRUE
    if (isTruthy(input$loadRDSIntegration) && isTruthy(input$inputRDSIntegration)) any_uploading <- TRUE
    if (isTruthy(input$loadRDSDownstream) && isTruthy(input$inputRDSDownstream)) any_uploading <- TRUE
    if (isTruthy(input$loadRDSDGEA) && isTruthy(input$inputRDSDGEA)) any_uploading <- TRUE
    if (isTruthy(input$loadRDSAnnotate) && isTruthy(input$inputRDSAnnotate)) any_uploading <- TRUE
    if (isTruthy(input$loadRDSFeaturePlot) && isTruthy(input$inputRDSFeaturePlot)) any_uploading <- TRUE

    if (!any_uploading) return(NULL)

    tags$div(
      style = "color: #856404; background-color: #fff3cd; border: 1px solid #ffeeba; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
      tags$strong("\u26A0\ufe0f Please wait: "),
      "Wait for RDS file upload to complete before proceeding."
    )
  })

  # output$ramBarUI <- renderUI({
  #   autoInvalidate()
  #   used_ram <- as.numeric(pryr::mem_used()) / 1024^2
  #   percent_used <- round(used_ram / total_ram_mb * 100)
  #
  #   # Dynamic color based on RAM usage
  #   bar_color <- if (percent_used < 50) {
  #     "green"
  #   } else if (percent_used < 80) {
  #     "orange"
  #   } else {
  #     "red"
  #   }
  #
  #   fluidRow(
  #     column(
  #       12,
  #       div(style = "position:fixed; bottom:15px; left:20px; right:20px; height:20px; z-index:9999; background-color:#e0e0e0; border-radius:10px; box-shadow: 0px 0px 8px rgba(0,0,0,0.3);",
  #           div(style = paste0(
  #             "width:", percent_used, "%;",
  #             "height:100%;",
  #             "background-color:", bar_color, ";",
  #             "text-align:center;",
  #             "color:white;",
  #             "font-weight:bold;",
  #             "font-size:12px;",
  #             "line-height:20px;"
  #           ),
  #           paste0(round(used_ram, 1), " MB (", percent_used, "% RAM)")
  #           )
  #       )
  #     )
  #   )
  # })




  observeEvent(input$readH5, {
    req(input$h5file)
    output$uploadStatus <- renderText("Reading .h5 file...")

    map_text <- input$sampleMap
    sample_map <- NULL
    if (nzchar(map_text)) {
      tryCatch({
        parts <- unlist(strsplit(map_text, ","))
        kv <- sapply(parts, function(x) unlist(strsplit(x, "=")))
        sample_map <- setNames(kv[2,], kv[1,])
      }, error = function(e) {
        output$uploadStatus <- renderText("Invalid sample map format.")
        return()
      })
    }

    tryCatch({
      rv$seurat <- read_h5_to_seurat(
        filepath = input$h5file$datapath,
        project_name = input$projectName
      )
      # Apply sample label remapping (done here to avoid server-side errors
      # from premature recoding inside the helper function)
      if (!is.null(sample_map) && length(sample_map) > 0) {
        rv$seurat@meta.data$SampleLabel <- dplyr::recode(rv$seurat@meta.data$SampleLabel, !!!sample_map)
      }
output$uploadStatus <- renderText("Seurat object successfully created.")
      output$uploadAlertH5 <- renderUI(NULL)
    }, error = function(e) {
      output$uploadStatus <- renderText(paste("Error reading file:", e$message))
    })
  })

  observeEvent(input$runQC, {
    req(rv$seurat)

    # Create QC output dir inside tempdir()
    qc_dir <- file.path(tempdir(), "qc_plots")
    dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

    # Store expected output paths in reactiveValues
    rv$qc_plot_paths <- list(
      violin = file.path(qc_dir, "pre_qc_violin.jpg"),
      percent_mt = file.path(qc_dir, "pre_qc_percent.mt.jpg"),
      scatter = file.path(qc_dir, "pre_qc_FeatureScatterplot.jpg")
    )

    # Reset status before re-running
    rv$qc_ready <- FALSE

    tryCatch({
      rv$seurat <- run_qc_metrics(
        seurat_obj = rv$seurat,
        mt_pattern = input$mtPattern,
        output_dir = qc_dir,
        sample_split_var = input$splitBy,
        ylimit = input$ylimit,
        save_plots = TRUE
      )

      # Check if all plots were saved
      all_exist <- all(file.exists(unlist(rv$qc_plot_paths)))

      if (all_exist) {
        rv$qc_ready <- TRUE
        rv$qc_counter <- rv$qc_counter + 1  # 🔁 Trigger UI refresh
      } else {
        warning("One or more QC plots are missing.")
      }

      # Optional: debug print
      print("Expected saved plots:")
      print(list.files(qc_dir, full.names = TRUE))

    }, error = function(e) {
      showNotification(paste("QC failed:", e$message), type = "error")
    })
  })




  observeEvent(input$postQC, {
    req(rv$seurat)
    output$postQCStatus <- renderText("Filtering...")

    # Create Post-QC output dir inside tempdir()
    postqc_dir <- file.path(tempdir(), "post_qc_outputs")
    dir.create(postqc_dir, showWarnings = FALSE, recursive = TRUE)

    # Save filtered Seurat to a temp RDS so the download button always works
    rv$postqc_rds_path <- file.path(postqc_dir, "seurat_post_qc_filtered.rds")

    # Expected output plot paths (written by post_qc_filtering_and_plots.R)
    rv$postqc_plot_paths <- list(
      vln = file.path(postqc_dir, "post_qc_vplot.jpg"),
      density = file.path(postqc_dir, "post_qc_densityplot.jpg"),
      pie = file.path(postqc_dir, "cellcount_piechart_post.jpg")
    )

    rv$postqc_ready <- FALSE

    tryCatch({
      rv$seurat <- post_qc_filtering_and_plots(
        seurat_obj = rv$seurat,
        min_features = input$minFeat,
        max_features = input$maxFeat,
        max_mt_percent = input$maxMT,
        output_dir = postqc_dir,
        sample_split_var = input$splitByPost,
        save_rds = TRUE,
        rds_path = rv$postqc_rds_path
      )

      all_exist <- all(file.exists(c(unlist(rv$postqc_plot_paths), rv$postqc_rds_path)))

      if (all_exist) {
        rv$postqc_ready <- TRUE
        rv$postqc_counter <- rv$postqc_counter + 1
      } else {
        warning("One or more Post-QC outputs are missing.")
      }

      output$postQCStatus <- renderText("Post-QC filtering complete.")
    }, error = function(e) {
      output$postQCStatus <- renderText(paste("Post-QC failed:", e$message))
    })
  })

observeEvent(input$removeDoublets, {
    output$doubletStatus <- renderText("Starting doublet removal...")

    tryCatch({
      # Use RDS file if selected
      if (input$loadRDSInput) {
        req(input$inputRDSfile)
        seurat_in <- readRDS(input$inputRDSfile$datapath)
      } else {
        req(rv$seurat)
        seurat_in <- rv$seurat
      }

      rv$seurat <- remove_doublets(
        seurat_obj = seurat_in,
        output_dir = input$doubletOutDir,
        save_rds = input$saveRDSDoublet,
        rds_path = input$rdsPathDoublet,
        doublet_rate = input$doubletRate
      )

      output$doubletStatus <- renderText("Doublets removed successfully.")
      output$uploadAlertRDS <- renderUI(NULL)
    }, error = function(e) {
      output$doubletStatus <- renderText(paste("Doublet removal failed:", e$message))
    })
  })

  observeEvent(input$runIntegration, {
    output$integrationStatus <- renderText("Running integration...")

    tryCatch({
      # Load Seurat object from file if checkbox is checked
      seurat_input <- NULL
      if (input$loadRDSIntegration) {
        req(input$inputRDSIntegration)
        seurat_input <- readRDS(input$inputRDSIntegration$datapath)
      } else {
        req(rv$seurat)
        seurat_input <- rv$seurat
      }

      #future_max_size_bytes <- as.numeric(input$futureMaxSize) * 1024^3

      withProgress(message = "Integrating samples...", value = 0, {

        # Apply memory limit BEFORE anything runs in parallel
        options(future.globals.maxSize = input$futureMaxSize * 1024^3)

        incProgress(0.1, detail = "Splitting Seurat object")
        seurat_list <- Seurat::SplitObject(seurat_input, split.by = "SampleLabel")

        incProgress(0.2, detail = "Running SCTransform")
        seurat_list <- lapply(seurat_list, function(x) {
          Seurat::SCTransform(x, verbose = FALSE)
        })

        incProgress(0.4, detail = "Selecting features")
        features <- Seurat::SelectIntegrationFeatures(seurat_list, nfeatures = input$nfeaturesIntegration)

        incProgress(0.5, detail = "Preparing SCT integration")
        seurat_list <- Seurat::PrepSCTIntegration(seurat_list, anchor.features = features)

        incProgress(0.7, detail = "Finding integration anchors")
        anchors <- Seurat::FindIntegrationAnchors(object.list = seurat_list,
                                                  normalization.method = "SCT",
                                                  anchor.features = features)

        incProgress(0.9, detail = "Integrating data")
        integrated <- Seurat::IntegrateData(anchorset = anchors, normalization.method = "SCT")

        if (input$saveRDSIntegration) {
          dir.create(dirname(input$rdsPathIntegration), recursive = TRUE, showWarnings = FALSE)
          saveRDS(integrated, input$rdsPathIntegration)
        }

        # Update reactive variable only if not loading from RDS
        rv$seurat <- integrated

        incProgress(1, detail = "Done")
      })

      output$integrationStatus <- renderText("Integration completed successfully.")
    }, error = function(e) {
      output$integrationStatus <- renderText(paste("Integration failed:", e$message))
    })
  })

  observeEvent(input$runDownstream, {
    output$downstreamStatus <- renderText("Running downstream analysis...")

    tryCatch({
      # Use RDS if selected
      if (input$loadRDSDownstream) {
        req(input$inputRDSDownstream)
        seurat_input <- readRDS(input$inputRDSDownstream$datapath)
      } else {
        req(rv$seurat)
        seurat_input <- rv$seurat
      }

      dims <- seq(input$dimsUsed[1], input$dimsUsed[2])

      # Build plot_options from UI inputs
      plot_options <- list()
      for (plot_id in input$selectedPlots) {
        plot_options[[plot_id]] <- list(
          save = TRUE,
          filename = input[[paste0(plot_id, "_filename")]],
          format = input[[paste0(plot_id, "_format")]],
          width = input[[paste0(plot_id, "_width")]],
          height = input[[paste0(plot_id, "_height")]],
          title = input[[paste0(plot_id, "_title")]]
        )
      }

      withProgress(message = "Running downstream analysis...", value = 0, {
        incProgress(0.2, detail = "Running clustering and plotting...")

        seurat_input <- run_downstream_analysis(
          seurat_obj = seurat_input,
          dims = dims,
          resolution = input$resolution,
          output_dir = input$downstreamOutDir,
          save_rds = input$saveRDSDownstream,
          rds_path = input$rdsPathDownstream,
          plot_options = plot_options
        )

        # Update reactive value if not loading from file
        if (!input$loadRDSDownstream) {
          rv$seurat <- seurat_input
        }

        incProgress(1)
      })

      output$downstreamStatus <- renderText("Downstream analysis complete.")
    }, error = function(e) {
      output$downstreamStatus <- renderText(paste("Downstream analysis failed:", e$message))
    })
  })


  observeEvent(input$runDGEA, {
    output$dgeaStatus <- renderText("Running DGEA...")

    tryCatch({
      # Load Seurat object if requested
      if (input$loadRDSDGEA) {
        req(input$inputRDSDGEA)
        seurat_input <- readRDS(input$inputRDSDGEA$datapath)
      } else {
        req(rv$seurat)
        seurat_input <- rv$seurat
      }

      withProgress(message = "Running cluster-wise DGEA...", value = 0, {
        incProgress(0.2, detail = "Running DGEA and extracting markers...")

        plot_options <- if (!input$savePlotsDGEA) NULL else list(
          save_violin = "violin" %in% input$dgeaPlotTypes,
          save_dot = "dot" %in% input$dgeaPlotTypes,
          save_feature = "feature" %in% input$dgeaPlotTypes,
          format = input$dgeaPlotFormat,
          width = input$dgeaPlotWidth,
          height = input$dgeaPlotHeight,
          title_prefix = input$dgeaTitlePrefix
        )

        run_clusterwise_dgea(
          seurat_obj = seurat_input,
          output_dir = input$dgeaOutputDir,
          logfc_threshold = input$logfcThreshold,
          min_pct = input$minPct,
          plot_options = plot_options
        )

        incProgress(1, detail = "Done!")
      })

      output$dgeaStatus <- renderText("Cluster-wise DGEA and plots completed successfully.")
    }, error = function(e) {
      output$dgeaStatus <- renderText(paste("DGEA failed:", e$message))
    })
  })

  observeEvent(input$runAnnotation, {
    output$annotationStatus <- renderText("Running annotation...")

    tryCatch({
      # Load Seurat object
      if (input$loadRDSAnnotate) {
        req(input$inputRDSAnnotate)
        seurat_input <- readRDS(input$inputRDSAnnotate$datapath)
      } else {
        req(rv$seurat)
        seurat_input <- rv$seurat
      }

      annotations_input <- NULL

      # Handle annotation input
      if (input$annotationInputType == "csv") {
        req(input$annotationCSV)
        annotations_input <- input$annotationCSV$datapath
      } else if (input$annotationInputType == "manual") {
        annotation_vector <- unlist(strsplit(input$annotationText, ","))
        annotation_vector <- trimws(annotation_vector)
        if (length(annotation_vector) != length(unique(seurat_input$seurat_clusters))) {
          stop(paste("Number of annotations (", length(annotation_vector),
                     ") must match number of clusters (",
                     length(unique(seurat_input$seurat_clusters)), ")."))
        }
        annotations_input <- annotation_vector
      }

      withProgress(message = "Annotating clusters...", value = 0, {
        incProgress(0.5, detail = "Applying annotations")

        annotated_obj <- annotate_clusters(
          seurat_obj = seurat_input,
          annotations = annotations_input,
          output_dir = input$annotOutputDir,
          rds_path = input$rdsPathAnnotate,
          save_rds = input$saveRDSAnnotate
        )

        # Update reactive value only if not using external RDS
        if (!input$loadRDSAnnotate) {
          rv$seurat <- annotated_obj
        }

        incProgress(1, detail = "Done")
      })

      output$annotationStatus <- renderText("Clusters annotated successfully.")
    }, error = function(e) {
      output$annotationStatus <- renderText(paste("Annotation failed:", e$message))
    })
  })

  observeEvent(input$runFeaturePlot, {
    output$featurePlotStatus <- renderText("Generating feature plots...")

    tryCatch({
      # Load Seurat object if needed
      if (input$loadRDSFeaturePlot) {
        req(input$inputRDSFeaturePlot)
        seurat_input <- readRDS(input$inputRDSFeaturePlot$datapath)
      } else {
        req(rv$seurat)
        seurat_input <- rv$seurat
      }

      # Parse genes
      genes <- unlist(strsplit(input$genesToPlot, ","))
      genes <- trimws(genes)
      genes <- genes[genes != ""]

      if (length(genes) == 0) {
        stop("Please enter at least one gene name.")
      }

      withProgress(message = "Plotting genes...", value = 0, {
        incProgress(0.3, detail = "Generating FeaturePlots")

        # Run and capture plots
        plot_list <- plot_genes_feature(
          seurat_obj = seurat_input,
          genes = genes,
          output_dir = input$featurePlotOutputDir,
          save_plots = input$saveFeaturePlots
        )

        # Render the plots in Shiny
        output$featurePlotOutputs <- renderUI({
          plot_output_list <- lapply(names(plot_list), function(gene) {
            plotname <- paste0("plot_", gene)
            plotOutput(plotname, height = "400px", width = "50%")
          })
          do.call(tagList, plot_output_list)
        })

        # Create individual renderPlot for each plot
        for (gene in names(plot_list)) {
          local({
            my_gene <- gene
            output[[paste0("plot_", my_gene)]] <- renderPlot({
              plot_list[[my_gene]]
            })
          })
        }

        incProgress(1, detail = "Done")
      })

      output$featurePlotStatus <- renderText("Feature plots generated and displayed successfully.")
    }, error = function(e) {
      output$featurePlotStatus <- renderText(paste("Feature plotting failed:", e$message))
    })
  })





}

shinyApp(ui, server)

# Run the app on all network interfaces, not just localhost
# shiny::runApp(host = "0.0.0.0", port = 5807)

# Link to share to others
# http://10.18.13.248:5807
