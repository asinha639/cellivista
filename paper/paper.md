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
    affiliation: 1
  - name: Kathryn Sinha
    orcid: 0000-0000-0000-0000 
    affiliation: 2
  - name: Richard A. McIndoe
    orcid: 0000-0002-3040-3631 
    corresponding: true
    affiliation: 1
affiliations:
  - index: 1
    name: Center for Biotechnology & Genomic Medicine, Augusta University, Augusta, GA, United States
  - index: 2
    name: Department of Biostatistics, Data Science, & Epidemiology, Augusta University, Augusta, GA, United States
date: 23 July 2026
bibliography: paper.bib
---

# Summary

Single-cell RNA sequencing (scRNA-seq) lets researchers measure gene activity in thousands of individual cells at once, revealing the different cell types and states that make up a tissue. Turning the raw measurements into biological insight, however, usually requires writing and editing computer code, which is a barrier for many laboratory scientists. `cellivista` is a web-browser application that guides a user through a complete scRNA-seq analysis without programming. Built in R with the Shiny framework on top of the widely used Seurat toolkit [@Satija:2015], it walks the user step by step from loading data, through quality control, removal of technical artefacts, combining multiple samples, grouping cells into clusters, identifying the genes that define each group, labelling those groups with cell-type names, and finally visualising the expression of genes of interest. At every stage, intermediate results can be saved and reloaded, so a long analysis can be paused and resumed. By exposing a rigorous, reproducible Seurat workflow through a point-and-click interface, `cellivista` makes modern single-cell analysis accessible to bench scientists while keeping the analytical steps transparent. The application is freely available as open-source software and can be run either through a hosted deployment or locally from the source repository. The archived software release is available at https://doi.org/10.5281/zenodo.21511557.

# Statement of need

Seurat [@Satija:2015] is one of the most widely adopted frameworks for scRNA-seq analysis, offering a flexible foundation for preprocessing, clustering, marker discovery, and visualisation. Despite its broad adoption, Seurat is used primarily by investigators who are comfortable writing and modifying R code. This requirement limits how easily many bench scientists can explore their own data independently, and it slows the iterative, parameter-tuning style of exploration (adjusting filtering thresholds, the number of dimensions, or the clustering resolution) that single-cell analysis demands. 

cellivista is designed for experimental biologists, bioinformatics core facilities, and other researchers who need to perform standard scRNA-seq analyses without developing custom R scripts. It packages the routine exploratory workflow into a guided tabbed interface. This includes input preparation and Seurat object creation, quality control metrics and filtering, doublet removal using DoubletFinder [@McGinnis:2019], SCTransform based multi sample integration, dimensionality reduction, clustering, marker identification, cluster annotation, and gene level visualization. It accepts 10x Genomics HDF5 (`.h5`) files, existing Seurat `.rds` objects, and matrix/barcode/feature files that it can convert to Seurat `.rds` format in-app. Because each step is backed by standard Seurat functions and intermediate objects are exportable, the analysis remains reproducible and auditable rather than hidden behind a black-box interface. 

# State of the field

Several graphical tools aim to make single-cell analysis more approachable. ASAP [@Gardeux:2017] is a web-based platform for the analysis and visualisation of scRNA-seq data, and the Broad Institute's Single Cell Portal [@Tarhan:2023] provides an interactive home for sharing and exploring single-cell genomics datasets. These platforms are valuable, but they can constrain customisation, depend on specific input structures or predefined workflows, or are oriented toward data sharing and visualisation of already-processed results rather than running the full upstream processing pipeline. At the other end of the spectrum, command-line Seurat offers complete flexibility but requires R programming proficiency. 

Rather than re-implementing single-cell algorithms, `cellivista` deliberately builds upon the Seurat ecosystem and exposes it through a browser interface. Its distinct contribution is delivering the complete exploratory Seurat workflow, including steps that hosted portals typically leave to the user's own scripts, such as doublet removal and SCT based integration, all within a single, no-code, locally runnable application. Throughout every stage, the underlying Seurat objects remain fully exportable. `cellivista` can be used either through a hosted deployment or by running the application locally from its source repository. Its computational workflow is implemented as modular helper functions that are separate from the reactive Shiny interface, allowing individual analytical steps to be inspected, tested, and reused independently while preserving the interactive user experience. 

# Software design

`cellivista`'s architecture separates the user interface from the analysis logic. The reactive Shiny front end (`app.R`) presents the workflow as a sequence of tabs, while the actual computation lives in a set of modular helper functions in the `R/` directory — for example `read_h5_to_seurat.R`, `convert_matrix_to_rds.R`, `run_qc_metrics.R`, `post_qc_filtering_and_plots.R`, `remove_doublets.R`, `integrate_seurat_samples.R`, `run_downstream_analysis.R`, `run_clusterwise_dgea.R`, `annotate_clusters.R`, and `plot_genes_feature.R`. This division was a deliberate trade-off: keeping each analytical step in a standalone, UI-independent function makes the code easier to test, reuse, and maintain separately from the interface, at the cost of a small amount of extra plumbing to wire functions to reactive inputs. 

![Overview of the cellivista workflow. Multiple input formats are accepted and processed through a guided Seurat-based analysis pipeline comprising quality control, normalization, dimensionality reduction, clustering, marker identification, and cell type annotation. Results, plots, and intermediate Seurat objects can be exported to support reproducible analyses.](figure1.png)

**Figure 1:** Overview of the cellivista workflow.

Although cellivista is presented through a web-browser interface, the analytical workflow is implemented entirely in R using Shiny, Seurat, and related packages. The user interface serves as a front end to modular analysis functions rather than embedding computational logic directly in the reactive application code. The hosted deployment provides convenient access for users, while the complete source code remains openly available so the application can be executed, inspected, modified, and reproduced locally.

A second design decision was to make intermediate Seurat objects first-class, exportable outputs. Rather than treating an analysis as a single end-to-end run, `cellivista` lets users save the Seurat object after major steps (e.g., post-filtering, post-integration) and reload it in a later tab. This supports the iterative, exploratory nature of single-cell analysis — users can revisit clustering parameters without rerunning earlier stages — and preserves reproducibility because the saved objects can be inspected or re-used in scripted Seurat workflows outside the app. The application can be run locally from the repository or deployed on a hosted Shiny server, and a production instance is maintained at <https://au-cbgm-shiny.augusta.edu/cellivista/>. 

# Research impact statement

`cellivista` was developed to support single-cell RNA-seq analysis within the Center for Biotechnology & Genomic Medicine (CBGM) Bioinformatics Core at Augusta University, where it lowers the barrier for investigators to explore their own datasets without depending on bespoke scripting. As a functionality demonstration, the application was used to reproduce a canonical analysis of a publicly available benchmark dataset (human peripheral blood mononuclear cells, PBMCs) from the Gene Expression Omnibus under accession GSE132044 [@Ding:2020]: after quality control, filtering, dimensionality reduction, and graph-based clustering, marker detection recovered well-established immune markers including *CD4*, *CD8A*, and *MS4A1*, consistent with expected T-cell and B-cell populations. This functionality demonstration shows that the guided workflow can recover expected biological populations and marker gene patterns from a well-characterized public dataset using the implemented Seurat-based analysis pipeline.

The software is openly developed as a modular R/Shiny application. The repository includes example data, local execution instructions, conversion utilities, and a detailed vignette to support reproducible analyses and allow reviewers and users to run the complete workflow locally. 

# AI usage disclosure

OpenAI ChatGPT (GPT-5 family, 2026) was used during the development of `cellivista` to assist with portions of R/Shiny code refinement, documentation, manuscript drafting and editing, and creation of the workflow schematic figure. AI assistance was also used to improve software documentation, including the README, vignette, and JOSS manuscript.

All AI-assisted code, figures, and text were reviewed, tested, validated, and edited by the authors before inclusion in the software or manuscript. All software architecture, analytical workflow design, scientific interpretation, and validation decisions were made by the authors, who take full responsibility for the accuracy, originality, licensing compliance, and integrity of the submitted software and manuscript. 

# Acknowledgements

We thank the CBGM Bioinformatics Core and the Augusta University research computing team for their feedback and support during development. This work was supported internally by the CBGM Bioinformatics Core at Augusta University; the sponsor had no involvement in the design, analysis, or preparation of this work. The authors declare no conflicts of interest. 
# References
