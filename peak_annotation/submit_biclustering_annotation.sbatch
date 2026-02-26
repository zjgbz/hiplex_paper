#!/bin/bash

#SBATCH --job-name=motif_analysis
#SBATCH -p shared
##SBATCH -p pi-sysadmin
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1-02:00:00
#SBATCH --mem=120GB
#SBATCH --output=./err/%x_%A_%a.out
#SBATCH --error=./err/%x_%A_%a.err
#SBATCH --mail-type=END,ALL
#SBATCH --mail-user=faai.jyu@gmail.com
# source activate /users/kyu1/miniconda/envs/snakemake

# tsv_dir=/dcs05/hongkai/data/next_cutntag/bulk/wgc/mixed/800/
# tsv_list=($(ls $tsv_dir/bicluster_V_mixed_all-qc_kmeans_euclidean_row_num-15_column_num-16_heatmap_row_clusters_extend*.tsv))
tsv_file=../data/peak_annotation/bicluster_V_mixed_all-qc_kmeans_euclidean_row_num-15_column_num-16_heatmap_row_clusters_extend_sparsity-10_cor_quantile-0.75_no_sparse.tsv
# /dcs05/hongkai/data/next_cutntag/bulk/wgc/mixed/800//bicluster_V_mixed_all-qc_kmeans_euclidean_row_num-15_column_num-16_heatmap_row_clusters_extend_sparsity-10_cor_quantile-0.75_w_sparse.tsv

# tsv_file=${tsv_list[${SLURM_ARRAY_TASK_ID}]}
basename_without_ext="$(basename "${tsv_file%.*}")"
# input_file=/dcs05/hongkai/data/next_cutntag/bulk/wgc/mixed/800/A.bed
out_dir="../data_peak_annotation/${basename_without_ext}.pdf"
# out_dir="/dcs05/hongkai/data/next_cutntag/bulk/wgc/mixed/800/annotations/${basename_without_ext}.pdf"

module load conda_R
Rscript /dcs05/hongkai/data/next_cutntag/script/peak_annotation/scripts_for_manuscript_v2/biclustering_annotation.R \
    $tsv_file \
    $out_dir