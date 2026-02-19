
read_dir="../data/frag_decomposition"
out_dir="../results/Figure2"
mkdir -p $out_dir


module load conda_R

Rscript ../hiplex_paper/fragment_decomposition/barplot.R \
    $read_dir \
    $out_dir