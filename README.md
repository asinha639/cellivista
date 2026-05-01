# Cellivista

<p align="center">
  <img src="www/logo.png" alt="Cellivista logo" width="280"/>
</p>

<p align="center">
  <b>Cellivista</b> is a Shiny application for end-to-end single-cell RNA-seq analysis using a Seurat-based workflow.
</p>

<p align="center">
  <a href="https://au-cbgm-shiny.augusta.edu/Cellivista/"><img alt="Live App" src="https://img.shields.io/badge/Live%20App-Available-1696d2"></a>
  <img alt="R" src="https://img.shields.io/badge/R-%3E%3D%204.3-276DC3">
  <img alt="Shiny" src="https://img.shields.io/badge/Shiny-Web%20Application-1f9ed8">
  <img alt="Status" src="https://img.shields.io/badge/Status-Server%20tested-success">
</p>

## Overview

Cellivista provides a guided, browser-based workflow for single-cell RNA-seq analysis, from input preparation and Seurat object creation through quality control, doublet removal, multi-sample integration, clustering, differential expression, cluster annotation, and FeaturePlot-based gene visualization. The app supports 10x Genomics HDF5 (`.h5`) files, existing Seurat `.rds` objects, and matrix/barcode/gene input files that can be converted to Seurat `.rds` format within the app.

### Main workflow modules

1. Convert matrix/barcode/gene files to a Seurat `.rds` object when needed
2. Upload 10x Genomics `.h5` data or an existing Seurat `.rds` object
3. Compute quality-control metrics and generate QC plots
4. Apply post-QC filtering and export filtered objects
5. Detect and remove doublets with DoubletFinder
6. Integrate multiple samples using SCTransform-based Seurat integration
7. Perform downstream dimensionality reduction, clustering, and UMAP visualization
8. Identify cluster-wise marker genes
9. Annotate clusters with user-provided cell type labels
10. Visualize gene expression with FeaturePlots

## Live deployment

Production instance:

**https://au-cbgm-shiny.augusta.edu/Cellivista/**

## Repository structure

## Repository structure

```text
Cellivista/
├── app.R
├── R/
│   ├── annotate_clusters.R
│   ├── convert_matrix_to_rds.R
│   ├── integrate_seurat_samples.R
│   ├── plot_genes_feature.R
│   ├── post_qc_filtering_and_plots.R
│   ├── read_h5_to_seurat.R
│   ├── remove_doublets.R
│   ├── run_clusterwise_dgea.R
│   ├── run_downstream_analysis.R
│   └── run_qc_metrics.R
├── www/
│   └── logo.png
├── inst/
│   └── extdata/
│       ├── sample_cluster_annotations.csv
│       └── GSE132044/
│           ├── celltype_annotations.csv     # example annotation file
│           ├── DOWNLOAD_INSTRUCTIONS.txt    # GEO download instructions for PBMC hg38 raw files
│           └── rds_convert.R                # RDS conversion script
├── docs/
│   └── Cellivista_vignette.pdf
```

## Requirements

### Core R packages

- shiny
- Seurat
- ggplot2
- dplyr
- shinyBS
- patchwork
- purrr
- readr
- future
- DoubletFinder
- Matrix
- shinyWidgets
- tibble
- tidyr

## Installation

### Option 1: Run locally from the repository

```r
install.packages(c(
  "shiny", "Seurat", "ggplot2", "dplyr", "shinyBS",
  "patchwork", "purrr", "readr", "future", "shinyWidgets",
  "tibble", "tidyr"
))
```

Install `DoubletFinder` according to its current repository instructions if it is not already available in your R environment.

Then launch the app:

```r
shiny::runApp(".")
```

## Input data

- Primary inputs:
  - 10x Genomics HDF5 (`.h5`) gene expression matrices
  - Previously saved Seurat `.rds` objects
- Additional supported inputs:
  - Matrix Market count matrices (`.mtx` or `.mtx.gz`) with corresponding barcode/cell and gene/features files, which can be converted to a Seurat `.rds` object within the app
- Optional inputs:
  - CSV-based cluster annotation maps for cell type labeling

A sample annotation template is included at:

```text
inst/extdata/sample_cluster_annotations.csv
```

## Outputs

Depending on user selections, Cellivista can generate:

- QC violin plots and feature-scatter plots
- Post-filtering diagnostic plots
- DoubletFinder summary plots
- Integrated Seurat objects (`.rds`)
- UMAP visualizations
- Cluster-wise DGEA result tables
- Annotation plots
- FeaturePlots for selected genes

## Citation

If you use Cellivista in your work, please cite the software using the repository information below:

> cellivista: An Interactive Seurat Based Shiny App for Single Cell RNA seq Analysis

A machine-readable citation file (`CITATION.cff`) is included in this repository.

**Note:** A formal citation (with DOI and versioned release) will be provided upon publication.

## Reproducibility and publication notes

For manuscript use, the following are recommended before public release:

- Tag a release version matching the manuscript submission
- Add final author list and affiliations to `CITATION.cff`
- Add a project license appropriate for your intended distribution
- Add screenshots or a workflow figure to this README
- Enable GitHub Releases and archive the release with Zenodo for a DOI

## Acknowledgments

Cellivista is built around the Seurat ecosystem for single-cell analysis and uses Shiny to provide an interactive web interface for non-programmatic use.
