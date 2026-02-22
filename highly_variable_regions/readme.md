# Hiplex Paper - Highly Variable Regions

## Overview
This script identifies highly variable genomic regions from whole genome correlation (WGC) data. It fits a mean-variance relationship model (GAM or LOESS) to normalized WGC signals, computes hypervariance scores for each genomic bin, and generates mean-variance and density plots for visualization.

## Prerequisites
Before running the analysis, ensure you have the following R packages installed:

- **arrow**
- **matrixStats**
- **mgcv**
- **ggplot2**
- **cowplot**
- **ggExtra**
- **ggpointdensity**

Note: The script will automatically install missing packages when run.

## Directory Structure
Make sure you are in the `hiplex_paper`directory:

`cd /path/to/your/hiplex_paper`

## How to Run
Run in terminal:

`bash highly_variable_regions/run_highly_variable_regions.sh` (approximately 10 minutes)

Note: Suggests request 128G

## Output/Results

The expected output directory is: 

`../results/Figure3A/`

Within this directory, you will find 

- `V1V2_mixed_800_colQC-all-qc_libnorm_noAllZero_log2_qnorm_gam-40.feather` — Table of per-bin statistics including mean, variance, expected variance, hypervariance, and clipped hypervariance scores.

- `V1V2_mixed_800_colQC-all-qc_libnorm_noAllZero_log2_qnorm_gam-40.png` — Combined plot of (A) log2 mean vs log2 variance with fitted curve, and (B) log2 mean vs hypervariance.

- `V1V2_mixed_800_colQC-all-qc_libnorm_noAllZero_log2_qnorm_gam-40_smooth-300000_seed-42.png` — Density scatter plot of log2 mean vs log2 variance with fitted curve.

- `result_gam_21_40_42_300000.out` — Log file with stdout/stderr from the run.

Note: Refers to Figure 3A in Paper