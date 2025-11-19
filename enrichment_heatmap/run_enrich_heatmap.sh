bw_dir="/Users/kyu/Documents/thesis/data/bigwig_normalized/V"
peak_type=data_peak_auc_0.05_extend_window_narrow_wide_adaptive_binomial_0.05

echo "peak_type: $peak_type"
peak_dir="/Users/kyu/Documents/thesis/data/peak/${peak_type}/V/"
out_dir=/Users/kyu/Documents/thesis/scripts/enriched_heatmap/result_pintersect_subsample/${peak_type}/
mkdir -p $out_dir
Rscript /Users/kyu/Documents/thesis/scripts/enriched_heatmap/code_row.R ${bw_dir} ${peak_dir} "H3K9me3" "H3K9ac" "#97BE5A" "#f5ac8e" "#2878B5" "1" $out_dir 0.1
Rscript /Users/kyu/Documents/thesis/scripts/enriched_heatmap/code_row.R ${bw_dir} ${peak_dir} "H3K27me3" "H3K27ac" "#97BE5A" "#f5ac8e" "#2878B5" "1" $out_dir 0.1
Rscript /Users/kyu/Documents/thesis/scripts/enriched_heatmap/code_row.R ${bw_dir} ${peak_dir} "POLR2AphosphoS2" "H3K4me3" "#97BE5A" "#f5ac8e" "#2878B5" "1" $out_dir
