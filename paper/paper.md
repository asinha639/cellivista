---
title: 'cellivista: An interactive Seurat-based Shiny interface for single-cell RNA-seq analysis'
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
bibliography: paper.bib
---

# Summary

Single-cell RNA sequencing (scRNA-seq) lets researchers measure gene activity in thousands of individual cells at once, revealing the different cell types and states that make up a tissue. Turning the raw measurements into biological insight, however, usually requires writing and editing computer code, which is a barrier for many laboratory scientists. Although several graphical interfaces for scRNA-seq analysis exist, many emphasize data exploration or visualization rather than providing a complete, reproducible end-to-end workflow built directly on Seurat. `cellivista` is a web-browser application that guides a user through a complete scRNA-seq analysis without programming. Built in R with the Shiny framework on top of the widely used Seurat toolkit [@Satija:2015], it walks the user step by step from loading data, through quality control, removal of technical artefacts, combining multiple samples, grouping cells into clusters, identifying the genes that define each group, labelling those groups with cell-type names, and finally visualising the expression of genes of interest. At every stage, intermediate results can be saved and reloaded, so a long analysis can be paused and resumed. By exposing a rigorous, reproducible Seurat workflow through a point-and-click interface, `cellivista` makes modern single-cell analysis accessible to bench scientists while keeping the analytical steps transparent. The application is freely available as open-source software and can be run either through a hosted deployment or locally from the source repository. 

# State of the field

A number of graphical tools have been developed to make single-cell RNA sequencing (scRNA-seq) analysis more accessible to researchers without extensive programming experience. ASAP [@Gardeux:2017] is a web-based platform that provides an integrated environment for post-alignment scRNA-seq analysis, including filtering, normalization, clustering, differential expression analysis, functional enrichment, and interactive visualization. The Single Cell Portal (SCP) [@Tarhan:2023] enables researchers to upload, share, visualize, and interactively explore single-cell datasets through a web interface, with an emphasis on data dissemination and collaborative exploration. In contrast, Seurat [@Satija:2015] provides one of the most comprehensive and flexible frameworks for scRNA-seq analysis, but its use generally requires writing and modifying R code.

Table: Comparison of `cellivista` with representative graphical scRNA-seq analysis platforms.

| Feature | Seurat | ASAP | Single Cell Portal | `cellivista` |
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

Rather than introducing new analytical methods, `cellivista` builds directly upon the established Seurat ecosystem and exposes its standard workflow through an intuitive browser-based interface. Its primary contribution is providing a complete, end-to-end exploratory workflow within a single application, encompassing data import, quality control, filtering, doublet detection using DoubletFinder [@McGinnis:2019], SCTransform-based multi-sample integration, dimensionality reduction, clustering, marker identification, cluster annotation, and gene-level visualization. Intermediate Seurat objects remain fully exportable throughout the workflow, allowing analyses to be paused, resumed, audited, or continued outside the graphical interface using standard Seurat workflows.

Unlike many hosted platforms, `cellivista` is designed to be deployed either locally or as a hosted web application, allowing researchers to analyze data within their own computing environments when required. In addition, the analytical workflow is implemented as a collection of modular R functions that are separate from the reactive Shiny interface. This separation enables individual computational components to be independently tested, maintained, and reused while preserving an interactive user experience. Consequently, `cellivista` functions not only as a graphical application for non-programming users but also as a reusable and extensible software framework for bioinformatics core facilities and developers.

# Software design

`cellivista`'s architecture separates the user interface from the analysis logic. The reactive Shiny front end (`app.R`) presents the workflow as a sequence of tabs, while the actual computation lives in a set of modular helper functions in the `R/` directory, each corresponding to a major stage of the analysis workflow, including data import, quality control, filtering, doublet detection, integration, downstream analysis, marker identification, annotation, and visualization. This division was a deliberate trade-off: keeping each analytical step in a standalone, UI-independent function makes the code easier to test, reuse, and maintain separately from the interface, at the cost of a small amount of extra plumbing to wire functions to reactive inputs. 

![Overview of the cellivista workflow.](figure1.png)

**Figure 1:** Overview of the cellivista workflow.

The workflow begins with data import from 10x Genomics HDF5 files, Seurat objects, or matrix/barcode/feature files. Users are then guided through quality control, normalization, dimensionality reduction, clustering, marker identification, and cell-type annotation. At each major stage, results, publication-quality figures, and intermediate Seurat objects can be exported to support reproducible and iterative analyses.

Although cellivista is presented through a web-browser interface, the analytical workflow is implemented entirely in R using Shiny, Seurat, and related packages. The user interface serves as a front end to modular analysis functions rather than embedding computational logic directly in the reactive application code. The hosted deployment provides convenient access for users, while the complete source code remains openly available so the application can be executed, inspected, modified, and reproduced locally.

A second design decision was to make intermediate Seurat objects first-class, exportable outputs. Rather than treating an analysis as a single end-to-end run, `cellivista` lets users save the Seurat object after major steps (e.g., post-filtering, post-integration) and reload it in a later tab. This supports the iterative, exploratory nature of single-cell analysis — users can revisit clustering parameters without rerunning earlier stages — and preserves reproducibility because the saved objects can be inspected or re-used in scripted Seurat workflows outside the app. The application can be run locally from the repository or deployed on a hosted Shiny server, and a production instance is maintained at <https://au-cbgm-shiny.augusta.edu/cellivista/>. 

# Research impact statement

`cellivista` was developed to support single-cell RNA-seq analysis within the Center for Biotechnology & Genomic Medicine (CBGM) Bioinformatics Core at Augusta University, where it lowers the barrier for investigators to explore their own datasets without depending on bespoke scripting. To demonstrate the software, the application was used to analyze a publicly available benchmark dataset (human peripheral blood mononuclear cells, PBMCs) from the Gene Expression Omnibus under accession GSE132044 [@Ding:2020]: after quality control, filtering, dimensionality reduction, and graph-based clustering, marker detection recovered well-established immune markers including *CD4*, *CD8A*, and *MS4A1*, consistent with expected T-cell and B-cell populations. The analysis recovered expected biological populations and marker gene expression patterns from a well-characterized public dataset, demonstrating that the guided workflow faithfully reproduces a standard Seurat-based scRNA-seq analysis.

The software is openly developed as a modular R/Shiny application with comprehensive documentation, automated unit testing, continuous integration through GitHub Actions, versioned releases archived with Zenodo, and example datasets that enable reviewers and users to reproduce analyses and execute the complete workflow locally.

# AI usage disclosure

OpenAI ChatGPT (GPT-5 family, 2026) was used during the development of `cellivista` to assist with portions of R/Shiny code refinement, documentation, manuscript drafting and editing, and creation of the workflow schematic figure. AI assistance was also used to improve software documentation, including the README, vignette, and JOSS manuscript.

All AI-assisted code, figures, and text were reviewed, tested, validated, and edited by the authors before inclusion in the software or manuscript. All software architecture, analytical workflow design, scientific interpretation, and validation decisions were made by the authors, who take full responsibility for the accuracy, originality, licensing compliance, and integrity of the submitted software and manuscript. 

# Acknowledgements

We thank the CBGM Bioinformatics Core and the Augusta University research computing team for their feedback and support during development. This work was supported internally by the CBGM Bioinformatics Core at Augusta University; the sponsor had no involvement in the design, analysis, or preparation of this work. The authors declare no conflicts of interest. 
# References
