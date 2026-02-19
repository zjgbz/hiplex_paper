# Hiplex Paper - Peak Count Heatmap

## Overview
This script is used to generate peak count heatmaps from genomic data. The heatmap visualizes the number of peaks for each target pair combination, with targets grouped by categories (Histone Modification, Writer, Transcription Factor).

## Prerequisites
Before running the analysis, ensure you have the following R packages installed:

- **ChIPseeker**
- **ComplexHeatmap**
- **glue**
- **latex2exp**
- **circlize**

Note: The script will automatically install missing packages when run.

## Directory Structure
Make sure you are in the `hiplex_paper`directory:

`cd /path/to/your/hiplex_paper`

## How to Run
Run in terminal:

`bash peak_count_heatmap/run_peak_count_heatmap.sh` (approximately 5 minutes)

## Output/Results

The expected output directory is: 

`../results/Figure2/`

Within this directory, you will find `V.pdf` - Heatmap showing log2-transformed peak counts for all target pair combinations

Note: Refers to Figure2B in Paper

