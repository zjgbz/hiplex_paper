# Pipeline for TF enrichment V0

bed_dir=/dcs05/hongkai/data/next_cutntag/bulk/wgc/mixed/800/bed/reorganized_wgc_final/
bed_list=($(ls $bed_dir/*.bed))
for bed_file in "${bed_list[@]}"; do
    # compute background region
    out_dir="/dcs05/hongkai/data/next_cutntag/bulk/motif_analysis/gw/"
    mkdir -p $out_dir
    bash /dcs05/hongkai/data/next_cutntag/script/motif_analysis/compute_background.sh $bed_file $out_dir
done

# merge background regions
module load conda_R
Rscript /dcs05/hongkai/data/next_cutntag/script/motif_analysis/merge_ctrl_seq.R /dcs05/hongkai/data/next_cutntag/bulk/motif_analysis/gw/ /dcs05/hongkai/data/next_cutntag/bulk/motif_analysis/gw/

# calculate TF enrichment
for bed_file in "${bed_list[@]}"; do
    echo $bed_file
    file_name=$(basename "${bed_file%.*}")
    out_dir="/dcs05/hongkai/data/next_cutntag/bulk/motif_analysis/gw/"
    control_bed_file=/dcs05/hongkai/data/next_cutntag/bulk/motif_analysis/gw/ctrl.bed
    Rscript /dcs05/hongkai/data/next_cutntag/script/motif_analysis/peak_analysis.R $bed_file $control_bed_file $peak_out_dir
done

# plot
mkdir -p ${out_dir}/${file_name}/peak_result/
peak_out_dir=${out_dir}/${file_name}/peak_result/result.tsv
module load conda_R
Rscript /dcs05/hongkai/data/next_cutntag/script/motif_analysis/peak_heatmap_bash.R /dcs05/hongkai/data/next_cutntag/bulk/motif_analysis/gw/ /dcs05/hongkai/data/next_cutntag/bulk/motif_analysis/hm/gw/peak_analysis/