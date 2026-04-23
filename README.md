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

Cellivista provides a guided, browser-based workflow for single-cell RNA-seq analysis, from raw 10x Genomics HDF5 input through quality control, doublet removal, multi-sample integration, clustering, differential expression, cluster annotation, and FeaturePlot-based gene visualization.

### Main workflow modules

1. Upload 10x Genomics `.h5` data and initialize a Seurat object
2. Compute quality-control metrics and generate QC plots
3. Apply post-QC filtering and export filtered objects
4. Detect and remove doublets with DoubletFinder
5. Integrate multiple samples using SCTransform-based Seurat integration
6. Perform downstream dimensionality reduction, clustering, and UMAP visualization
7. Identify cluster-wise marker genes
8. Annotate clusters with user-provided cell type labels
9. Visualize gene expression with FeaturePlots

## Live deployment

Production instance:

**https://au-cbgm-shiny.augusta.edu/Cellivista/**

## Repository structure

```text
Cellivista/
├── app.R
├── R/
│   ├── read_h5_to_seurat.R
│   ├── run_qc_metrics.R
│   ├── post_qc_filtering_and_plots.R
│   ├── remove_doublets.R
│   ├── integrate_seurat_samples.R
│   ├── run_downstream_analysis.R
│   ├── run_pairwise_dgea.R
│   ├── run_clusterwise_dgea.R
│   ├── annotate_clusters.R
│   └── plot_genes_feature.R
├── www/
│   └── logo.png
├── inst/
│   └── extdata/
│       └── sample_cluster_annotations.csv│       
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

- Primary input: 10x Genomics HDF5 (`.h5`) gene expression matrices
- Optional inputs: previously saved Seurat `.rds` objects and CSV-based cluster annotation maps

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

## Suggested citation text

If you cite the software in a manuscript, a concise software citation can follow this format:

> Cellivista: a Shiny-based interface for Seurat-driven single-cell RNA-seq analysis, including quality control, doublet removal, sample integration, clustering, differential expression, annotation, and gene-level visualization.

A machine-readable citation file is included as `CITATION.cff` and should be updated with the final author list, version, DOI, and GitHub URL before release.

## Reproducibility and publication notes

For manuscript use, the following are recommended before public release:

- Tag a release version matching the manuscript submission
- Add final author list and affiliations to `CITATION.cff`
- Add a project license appropriate for your intended distribution
- Add screenshots or a workflow figure to this README
- Enable GitHub Releases and archive the release with Zenodo for a DOI

## Acknowledgments

Cellivista is built around the Seurat ecosystem for single-cell analysis and uses Shiny to provide an interactive web interface for non-programmatic use.
