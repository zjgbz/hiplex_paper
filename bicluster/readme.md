# Hiplex Paper - Bicluster

## Overview
This script performs biclustering analysis on whole genome correlation (WGC) data from CUT&Tag experiments. It generates heatmaps with k-means clustering on both rows (genomic regions) and columns (target pairs), and outputs cluster assignments for downstream analysis.

## Prerequisites
Before running the analysis, ensure you have the following R packages installed:

- **ComplexHeatmap**
- **rtracklayer**
- **arrow**
- **tibble**
- **latex2exp**
- **glue**

Note: The script will automatically install missing packages when run.

## Directory Structure
Make sure you are in the `hiplex_paper`directory:

`cd /path/to/your/hiplex_paper`

## How to Run
Run in terminal:

`bash bicluster/run_bicluster.sh` (approximately 10 minutes)

## Output/Results

The expected output directory is: 

`../results/Figure3/`

Within this directory, you will find `bicluster_V_mixed_all-qc_kmeans_euclidean_row_num-15_column_num-16_heatmap_row-LETTERS_dend-off_manually_reorganized_2026-02-19_02-36-32-EST.pdf`

Note: Refers to Figure3A in Paper
