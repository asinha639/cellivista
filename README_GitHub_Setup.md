# GitHub setup checklist

## Suggested repository name

`cellivista`

## Suggested repository description

A Shiny-based application for end-to-end single-cell RNA-seq analysis using Seurat, including data conversion, quality control, doublet removal, integration, clustering, differential expression, annotation, and visualization.

## Suggested topics

`scrna-seq`, `single-cell`, `shiny`, `seurat`, `rstats`, `bioinformatics`, `transcriptomics`, `data-visualization`

## Suggested first release tag

`v0.1.0`

## After creating the remote repository

```bash
git init
git branch -M main
git add .
git commit -m "Initial public Cellivista release with RDS conversion and full analysis workflow"
git remote add origin https://github.com/asinha639/cellivista.git
git push -u origin main