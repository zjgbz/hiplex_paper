# Hiplex Paper - Fragment Decomposition

## Overview
This script is used to analyzes fragment length distributions from CUT&Tag data, performing differential analysis between conditions and generating visualizations including density plots, volcano plots, and barplots for different nucleosomal fragment categories (sub-nucleosomal, mono-nucleosomal, di-nucleosomal).

## Prerequisites
Before running the analysis, ensure you have the following R packages installed:

- **ggplot2**
- **ggrepel**
- **glue**
- **latex2exp**
- **cowplot**
- **dplyr**
- **sfsmisc**

Note: The script will automatically install missing packages when run.

## Directory Structure
Make sure you are in the `hiplex_paper`directory:

`cd /path/to/your/hiplex_paper`

## How to Run
Run in terminal:

1. To reproduce Figure2D: `bash fragment_decomposition/run_fragment_decomposition.sh` (approximately 2 minutes)

2. To reproduce all results: `bash fragment_decomposition/run_all_plots.sh` (approximately 5 minutes)

## Output/Results

The expected output directory is: 

`../results/Figure2/`

Within this directory, you will find: `barplot_frag_split_fastq-demux_per_V_valley-all-qc_summary.pdf`

Note: Refers to Figure2D in Paper
