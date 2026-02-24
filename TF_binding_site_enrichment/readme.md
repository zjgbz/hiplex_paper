# Hiplex Paper - TF Binding Site Enrichment

## Overview
This pipeline tests whether specific transcription factor (TF) binding motifs are enriched in genomic regions of interest compared to matched background control regions. It uses cisGenome to generate gene-distance-matched control regions, computes motif enrichment odds ratios and FDRs via Fisher's exact test, and visualizes the results as heatmaps.

The pipeline consists of four steps:
1. **Compute background regions** — For each input BED file, generate matched genomic control regions using cisGenome's `refgene_getmatchedcontrol`.
2. **Merge background regions** — Combine and deduplicate all control regions into a single unified background BED file (`ctrl.bed`).
3. **Calculate TF enrichment** — For each input BED file, compute motif enrichment (odds ratio, FDR) against the merged background.
4. **Plot heatmaps** — Generate heatmaps of log2(odds ratio) and z-scores across all clusters.

## Prerequisites

### Software
- **cisGenome 2.0** — Download and compile from [https://www.biostat.jhsph.edu/~hji/cisgenome/index.htm](https://www.biostat.jhsph.edu/~hji/cisgenome/index.htm). See installation instructions in the main README.

Installation Instruction:
`cd /path/to/your/paper_materials`
`wget http://jilab.biostat.jhsph.edu/software/cisgenome/executables/cisgenome_v2.0_linux.tar.gz`
`gzip -d cisgenome_v2.0_linux.tar.gz`
`tar xvf cisgenome_v2.0_linux.tar`
`cd cisgenome_project`
`./makefile`
`ls bin/`
`export PATH=/path/to/your/paper_materials/cisgenome_project/bin:$PATH`

### R Packages
The scripts will automatically install missing packages. Required packages include:

- **rtracklayer** (Bioconductor)
- **GenomicRanges** (Bioconductor)
- **plyranges** (Bioconductor)
- **ChIPseeker** (Bioconductor)
- **ComplexHeatmap** (Bioconductor)
- **matrixTests**
- **latex2exp**
- **glue**
- **circlize**

### Reference Data
The following reference files are expected in `../data/TF_binding_site_enrichment/`:

- `refFlat_sorted.txt` — Gene annotation file (refFlat format)
- `chrlen.txt` — Chromosome lengths file
- `bed_files/*.bed` — Input BED files (one per cluster/region set)

## Directory Structure
Make sure you are in the `hiplex_paper` directory:

`cd /path/to/your/hiplex_paper`

## How to Run
Make sure you install cisGenome First.

Run in terminal:

`bash TF_binding_site_enrichment/run_TF_bindingsite_enrichment.sh`

## Output/Results

The expected output directory is:

`../results/Figure3C/`

Within this directory, you will find:

-`reorganized_heatmap_filtered_rna_subtraction_5_t_zscore_clustered`

- `peak_no_perturb.pdf` — Full heatmap of log2(odds ratio) for all motifs passing filters (odds ratio ≥ 2, FDR ≤ 0.05).
- `peak_no_perturb.csv` — CSV of log2(odds ratio) values in clustered row order.
- `peak_no_perturb_fdr.csv` — CSV of FDR values in the same row order.
- `reorganized_heatmap_filtered_rna_subtraction_5_t.pdf` — Heatmap of top 5 enriched motifs per cluster (log2 odds ratio).
- `reorganized_heatmap_filtered_rna_subtraction_5_t_zscore.pdf` — Same heatmap with z-score normalization.
- `reorganized_heatmap_filtered_rna_subtraction_5_t_zscore_clustered.pdf` — Z-score heatmap with hierarchical clustering on columns.
- `reorganized_heatmap_filtered_rna_subtraction_5_t_cluster.pdf` — Log2(odds ratio) heatmap with column clustering.
- `{cluster_name}/peak_result/result.tsv` — Per-cluster motif enrichment results (odds ratio, FDR, target/background hit counts).

Intermediate data files are saved in `../data/TF_binding_site_enrichment/gw/`:

- `{bed_name}/control/{bed_name}_output_matchcontrol.bed` — Matched background regions per input BED file.
- `ctrl.bed` — Merged background control regions used for all enrichment tests.

Note: Refers to Figure 3C in Paper