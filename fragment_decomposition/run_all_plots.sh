
peak_dir="../data/frag_decomposition"
out_dir="../results/Figure2"
mkdir -p $out_dir


module load conda_R

Rscript ../hiplex_paper/fragment_decomposition/deco_density_valley.R \
    $peak_dir \
    $out_dir

Rscript ../hiplex_paper/fragment_decomposition/diff_frag_type.R \
    $peak_dir \
    $out_dir

Rscript ../hiplex_paper/fragment_decomposition/frag_hist_denstiy.R \
    $peak_dir \
    $out_dir

Rscript ../hiplex_paper/fragment_decomposition/barplot.R \
    $peak_dir \
    $out_dir