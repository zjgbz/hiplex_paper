# Hiplex Paper - Whole Genome Count

## Overview
This script generates whole genome count matrices from BAM alignment files. It bins the genome into fixed-size windows and counts fragment overlaps for each target pair combination, then applies various normalization methods (library normalization, log2 transformation, min-max scaling, square root transformation) for downstream analysis.

## Prerequisites
Before running the analysis, ensure you have the following R packages installed:

- **data.table**
- **arrow**
- **GenomicAlignments**
- **GenomicRanges**
- **Biostrings**
- **BSgenome.Hsapiens.NCBI.GRCh38**
- **plyranges**
- **preprocessCore**
- **tibble**

Note: The script will automatically install missing packages when run.

## Directory Structure
Make sure you are in the `hiplex_paper`directory:

`cd /path/to/your/hiplex_paper`

## How to Run
Run in terminal:

`bash whole_genome_count/run_whole_genome_count.sh` (approximately 18:40 minutes)

## Output/Results

The expected output directory is: 

`../results/Figure3/`

Note: Refers to Figure3A in Paper

