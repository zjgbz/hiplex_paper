

read_dir="../data/bicluster"
out_dir="../results/Figure3"
mkdir -p $out_dir


module load conda_R

Rscript ../hiplex_paper/bicluster/clusterheatmap_split_built-in.R \
    $read_dir \
    $out_dir