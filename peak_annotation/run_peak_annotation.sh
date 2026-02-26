
module load conda_R


peak_dir="/dcs07/hongkai/data/mjiang/hiplex/paper_materials/data/peak/data_peak_auc_0.05_extend_window_narrow_wide_adaptive_binomial_0.05"
out_dir="../results/Figrue2C_3B"
mkdir -p $out_dir


tsv_file=../data/peak_annotation/bicluster_V_mixed_all-qc_kmeans_euclidean_row_num-15_column_num-16_heatmap_row_clusters_extend_sparsity-10_cor_quantile-0.75_no_sparse.tsv
basename_without_ext="$(basename "${tsv_file%.*}")"

output_dir="../data_peak_annotation/${basename_without_ext}.pdf"

Rscript ../hiplex_paper/peak_annotation/scripts_for_manuscript_v2/biclustering_annotation.R \
    $tsv_file \
    $output_dir

# peak annotation - Chipseeker
output_dir="${out_dir}/Chipseeker"
mkdir -p $output_dir
ccre_rds_file="${output_dir}/annotate.rds"

Rscript ../hiplex_paper/peak_annotation/peakAnnotationChIPseekerCCREForLoop_different_peak_type.R \
    $peak_dir \
    $output_dir

Rscript ../hiplex_paper/peak_annotation/plot_result_all_chipseeker_ccre_bash.R \
    $output_dir \
    $output_dir 

output_dir="${out_dir}/RepeatMasker"
mkdir -p $output_dir
repeatmasker_rds_file="${output_dir}/annotate.rds"

Rscript ../hiplex_paper/peak_annotation/peakAnnotationRepeatMasker_different_peak_types.R \
    $peak_dir \
    $output_dir

new_order_dir="${output_dir}/CCRE_anno_V_order.rds"
Rscript ../hiplex_paper/peak_annotation/plot_result_all_repeatmasker_bash.R \
    $output_dir \
    $new_order_dir \
    $output_dir

Rscript ../hiplex_paper/peak_annotation/annotation_summary_bash.R \
    $peak_dir \
    $repeatmasker_rds_file \
    $ccre_rds_file \
    $out_dir



