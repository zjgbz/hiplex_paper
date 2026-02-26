
read_dir="../data/differential_analysis"
out_dir="../results/FigureS5"
mkdir -p $out_dir

module load conda_R

Rscript ../hiplex_paper/differential_analysis/limma_column_cluster_orig_post-filter.R \
    ${read_dir} \
    ${out_dir}

Rscript ../hiplex_paper/differential_analysis/df_region_scatter_all_points.R \
    ${read_dir} \
    ${out_dir}

Rscript ../hiplex_paper/differential_analysis/heatmap_differential.R \
    ${read_dir} \
    ${out_dir}

result_dir="../results/FigureS5/column_cluster_result"
summary_file="${result_dir}/summary_sig_num_0.25_0.25_l2fc-0.5.tsv"

echo -e "column_cluster\tsig_num" > ${summary_file}
for i in $(seq 1 16); do
    file="${result_dir}/result_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-0.25_column_cluster-${i}_limma_FDR-0.25_logFC-0.5.tsv"
    if [ -f "$file" ]; then
        # subtract 1 for header
        count=$(($(wc -l < "$file") - 1))
    else
        count=0
    fi
    echo -e "${i}\t${count}" >> ${summary_file}
done

Rscript ../hiplex_paper/differential_analysis/summary_sig_num_barplot_single_case.R \
    ${read_dir} \
    ${out_dir}

Rscript ../hiplex_paper/differential_analysis/diff_heatmap_colcluster_equal_width.R \
    ${read_dir} \
    ${out_dir}

