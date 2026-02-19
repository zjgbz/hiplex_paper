


peak_dir="../data/peak/data_peak_auc_0.05_extend_window_narrow_wide_adaptive_binomial_0.05"
out_dir="../results/Figure2"
mkdir -p $out_dir

module load conda_R
Rscript /dcs07/hongkai/data/mjiang/hiplex/paper_materials/hiplex_paper/peak_count_heatmap/peak_heatmap_bash.R \
    $peak_dir \
    $out_dir