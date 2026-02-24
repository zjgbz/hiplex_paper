# Hiplex Paper - Differential Analysis

## Overview
This pipeline performs differential analysis on whole genome correlation (WGC) data between treated and untreated conditions (V1/V2 vs T1/T2). It uses limma-voom to identify differentially co-bound genomic regions within each column cluster, generates scatter plots comparing Hi-Plex log2FC with gene expression log2FC, produces heatmaps of significant regions, and summarizes the number of significant regions per cluster as bar plots.

## Prerequisites
Before running the analysis, ensure you have the following R packages installed:


- **arrow**
- **edgeR** (Bioconductor)
- **limma** (Bioconductor)
- **GenomicRanges** (Bioconductor)
- **ComplexHeatmap** (Bioconductor)
- **matrixStats**
- **tibble**
- **glue**
- **latex2exp**
- **ggplot2**
- **ggpointdensity**
- **ggpubr**
- **ggrastr**
- **viridis**
- **svglite**
- **circlize**
- **data.table**
- **dplyr**
- **tidyr**

Note: The script will automatically install missing packages when run.

## Directory Structure
Make sure you are in the `hiplex_paper`directory:

`cd /path/to/your/hiplex_paper`

## How to Run
Run in terminal:

`bash differential_analysis/run_differential_analysis.sh` (approximately less than 1 hour)

Suggests request 256G

## Output/Results

The expected output directory is: 

`../results/FigureS5/column_cluster_fig`

Within this directory, you will find 

- `column_cluster_fig/` — Figures:
  - `limma_post-limmanorm_..._column_cluster-{cluster_id}_MA_plot.pdf` — MA plots per cluster.
  - `limma_post-limmanorm_..._column_cluster-{cluster_id}_significance_plot.pdf` — P-value/FDR histograms per cluster.
  - `scatter_{cluster_id}_5000_l2fc-0.5.pdf` — Scatter plots of Hi-Plex log2FC vs gene expression log2FC.
  - `publish_heatmap_..._column_cluster-{cluster_id}_limma_FDR-0.25_logFC-0.5_legend-on_colnames-on.pdf` — Heatmaps of significant regions.
  - `summary_sig_num_0.25_0.25_l2fc-0.5.pdf` — Bar plot of significant region counts per cluster.

Note: Refers to Figure S5 in Paper