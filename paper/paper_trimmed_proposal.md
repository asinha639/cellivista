---
title: 'Cellivista: An interactive Seurat-based Shiny interface for single-cell RNA-seq analysis'
tags:
  - R
  - Shiny
  - single-cell RNA-seq
  - scRNA-seq
  - Seurat
  - bioinformatics
authors:
  - name: Arnoneel Sinha
    orcid: 0000-0003-1006-0683
    corresponding: true
    affiliation: 1
  - name: Kathryn Sinha
    orcid: 0009-0006-0978-1087
    affiliation: 2
  - name: Richard A. McIndoe
    orcid: 0000-0002-3040-3631
    affiliation: 1
affiliations:
  - index: 1
    name: Center for Biotechnology & Genomic Medicine, Augusta University, Augusta, GA, United States
  - index: 2
    name: Department of Biostatistics, Data Science, & Epidemiology, Augusta University, Augusta, GA, United States
date: 24 July 2026
# TODO(me): set the date to the release date.
bibliography: paper.bib
---

<!-- Trimmed proposal for paper.md. Not part of any build. Only deletions and merged sentences; no new claims or citations. -->

# Summary

Single-cell RNA sequencing (scRNA-seq) lets researchers measure gene activity in thousands of individual cells at once, revealing the different cell types and states that make up a tissue. Cellivista is a web-browser application that guides a user through a complete scRNA-seq analysis without programming. Built in R with the Shiny framework on top of the widely used Seurat toolkit [@Satija:2015], it walks the user step by step from loading data to labelling cell types and visualising the expression of genes of interest.

# Statement of need

Turning the raw measurements into biological insight usually requires writing and editing computer code, which is a barrier for many laboratory scientists. Although several graphical interfaces for scRNA-seq analysis exist, many emphasize data exploration or visualization rather than providing a complete, reproducible end-to-end workflow built directly on Seurat. By exposing a rigorous, reproducible Seurat workflow through a point-and-click interface, Cellivista makes modern single-cell analysis accessible to bench scientists while keeping the analytical steps transparent.

Rather than introducing new analytical methods, Cellivista builds directly upon the established Seurat ecosystem. Its primary contribution is providing a complete, end-to-end exploratory workflow within a single application, encompassing data import, quality control, filtering, doublet detection using DoubletFinder [@McGinnis:2019], SCTransform-based multi-sample integration, dimensionality reduction, clustering, marker identification, cluster annotation, and gene-level visualization. Unlike many hosted platforms, Cellivista can be deployed either locally or as a hosted web application, allowing researchers to analyze data within their own computing environments when required.

# State of the field

ASAP [@Gardeux:2017] is a web-based platform that provides an integrated environment for post-alignment scRNA-seq analysis, including filtering, normalization, clustering, differential expression analysis, functional enrichment, and interactive visualization. The Single Cell Portal (SCP) [@Tarhan:2023] enables researchers to upload, share, visualize, and interactively explore single-cell datasets through a web interface, with an emphasis on data dissemination and collaborative exploration. In contrast, Seurat [@Satija:2015] provides one of the most comprehensive and flexible frameworks for scRNA-seq analysis, but its use generally requires writing R code.

Table: Comparison of Cellivista with representative graphical scRNA-seq analysis platforms.

| Feature | Seurat | ASAP | Single Cell Portal | Cellivista |
|:---|:---:|:---:|:---:|:---:|
| Browser-based graphical interface | ✗ | ✓ | ✓ | ✓ |
| Local deployment | ✓ | ✗ | ✗ | ✓ |
| Complete downstream scRNA-seq workflow | ✓ | ✓ | Partial | ✓ |
| Quality control and filtering | ✓ | ✓ | Partial | ✓ |
| Doublet detection workflow | Scripted | ✗ | ✗ | ✓ |
| SCTransform integration | ✓ | ✗ | ✗ | ✓ |
| Intermediate analysis export | ✓ | ✗ | ✗ | ✓ |
| Native Seurat object compatibility | ✓ | Limited | ✗ | ✓ |
| No programming required | ✗ | ✓ | ✓ | ✓ |

<!-- TODO(me): check other R/Shiny front-ends to Seurat for this comparison table (no tools or citations were added here); add any you verify, with citations. -->

# Software design

Cellivista's architecture separates the user interface from the analysis logic. The reactive Shiny front end (`app.R`) presents the workflow as a sequence of tabs, while the computation lives in modular, UI-independent functions in the `R/` directory, one for each major stage of the analysis. This division was a deliberate trade-off: it makes the code easier to test, reuse, and maintain separately from the interface, at the cost of a small amount of extra plumbing to wire functions to reactive inputs.

![Overview of the Cellivista workflow.](figure1.png)

The workflow begins with data import from 10x Genomics HDF5 files, Seurat objects, or matrix/barcode/feature files, followed by quality control, normalization, dimensionality reduction, clustering, marker identification, and cell-type annotation. An optional Visualization Options tab adds colorblind-safe UMAP palettes (Okabe-Ito and viridis), a heatmap of the top 3, 5, 10 or 15 marker genes per cluster drawn from the cluster-wise marker results, and a stacked bar chart of cluster composition per sample. For comparisons between two groups, the package also provides a pairwise differential expression function, `run_pairwise_dgea()`, which is called from R rather than from the interface.

A second design decision was to make intermediate Seurat objects first-class, exportable outputs. Users can save the Seurat object after major steps (e.g., post-filtering, post-integration) and reload it in a later tab, so clustering parameters can be revisited without rerunning earlier stages, and the saved objects can be inspected or re-used in scripted Seurat workflows outside the app. A production instance is maintained at <https://au-cbgm-shiny.augusta.edu/cellivista/>, and the source code is available at <https://github.com/asinha639/cellivista> under the MIT License.

# Research impact statement

Cellivista was developed to support single-cell RNA-seq analysis within the Center for Biotechnology & Genomic Medicine (CBGM) Bioinformatics Core at Augusta University. To demonstrate the software, the application was used to analyze a publicly available benchmark dataset (human peripheral blood mononuclear cells, PBMCs) from the Gene Expression Omnibus under accession GSE132044 [@Ding:2020]: after quality control, filtering, dimensionality reduction, and graph-based clustering, marker detection recovered well-established immune markers including *CD4*, *CD8A*, and *MS4A1*, consistent with expected T-cell and B-cell populations.

The software is openly developed as a modular R/Shiny application with comprehensive documentation, automated unit testing, continuous integration through GitHub Actions, versioned releases archived with Zenodo, and example datasets that enable reviewers and users to reproduce analyses and execute the complete workflow locally.

# AI usage disclosure

OpenAI ChatGPT (GPT-5 family, 2026) was used during the development of Cellivista to assist with portions of R/Shiny code refinement, documentation, manuscript drafting and editing, and creation of the workflow schematic figure. AI assistance was also used to improve software documentation, including the README, vignette, and JOSS manuscript. Anthropic Claude (Claude Code) was used for code refinement, test development, documentation, and manuscript editing.

<!-- TODO(me): confirm the wording of the Claude (Claude Code) disclosure above, including whether to state a model name or version. -->

All AI-assisted code, figures, and text were reviewed, tested, validated, and edited by the authors before inclusion in the software or manuscript. All software architecture, analytical workflow design, scientific interpretation, and validation decisions were made by the authors, who take full responsibility for the accuracy, originality, licensing compliance, and integrity of the submitted software and manuscript.

# Acknowledgements

We thank the CBGM Bioinformatics Core and the Augusta University research computing team for their feedback and support during development. This work was supported internally by the CBGM Bioinformatics Core at Augusta University; the sponsor had no involvement in the design, analysis, or preparation of this work. The authors declare no conflicts of interest.

# References
