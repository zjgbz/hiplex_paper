# Hiplex Paper - Peak Annotation

## Overview
This pipeline performs comprehensive genomic peak annotation for Hi-Plex CUT&Tag data. It annotates co-bound genomic regions using multiple reference databases — CCRE (cis-regulatory elements), ChromHMM, and RepeatMasker — and integrates with RNA-seq expression data and DNA methylation levels. The pipeline generates stacked bar plots of annotation distributions, hierarchically clustered dendrograms of target pairs, and combined summary figures showing CCRE, RepeatMasker, and methylation profiles side by side.

## Prerequisites
Before running the analysis, ensure you have the following R packages installed:

- **ggplot2**
- **cowplot**
- **GenomicRanges** (Bioconductor)
- **ChIPseeker** (Bioconductor)
- **ComplexHeatmap** (Bioconductor)
- **TxDb.Hsapiens.UCSC.hg38.knownGene** (Bioconductor)
- **rtracklayer** (Bioconductor)
- **plyranges** (Bioconductor)
- **glue**
- **RJSONIO**
- **tidyverse**
- **reshape2**
- **dplyr**
- **data.table**
- **arrow**
- **grid**
- **gridExtra**
- **stats**
- **ggrepel**
- **latex2exp**
- **dendextend**
- **stringr**
- **tools**
- **svglite**

Note: The script will automatically install missing packages when run.

## Directory Structure
Make sure you are in the `hiplex_paper` directory:

`cd /path/to/your/hiplex_paper`

## How to Run
Run in terminal:

`bash peak_annotation/run_peak_annotation.sh` (approximately 1 hour)

Suggests request 64G

This executes the following steps:

1. **Biclustering annotation** — Annotates biclustered genomic regions with CCRE, ChromHMM, and RepeatMasker, and integrates RNA-seq expression data per cluster.
2. **ChIPseeker + CCRE annotation** — Annotates peaks for all pairwise histone mark / TF combinations using CCRE cis-regulatory element overlaps and saves results as RDS.
3. **CCRE visualization** — Performs hierarchical clustering on CCRE annotation profiles and generates clustered stacked bar plots and dendrograms.
4. **RepeatMasker annotation** — Annotates the same peaks with RepeatMasker repetitive element classes and saves results as RDS.
5. **RepeatMasker visualization** — Generates clustered stacked bar plots of repetitive element distributions using the CCRE-derived clustering order.
6. **Annotation summary** — Combines CCRE annotations, RepeatMasker annotations, and DNA methylation M-value distributions into a single summary figure.

## Output/Results

The expected output directory is:

`../results/Figure2C_3B`

Within this directory, you will find:

- `bicluster_....pdf` — Combined annotation bar plots (CCRE, ChromHMM) and RNA expression boxplots per bicluster. (Figure 3B)
- `summary.pdf` — Combined figure with CCRE, RepeatMasker, and DNA methylation violin plots. (Figure 2C)

Note: Refers to Figures 2C and 3B in the paper.