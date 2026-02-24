
# require cisgenome from https://www.biostat.jhsph.edu/~hji/cisgenome/index.htm

cd /dcs07/hongkai/data/mjiang/hiplex/paper_materials/hiplex_paper

cisgenome_bin="../cisgenome_project/bin"

module load conda_R

read_dir="../data/TF_binding_site_enrichment/gw"
data_out_dir="../data/TF_binding_site_enrichment/gw"
mkdir -p $data_out_dir
out_dir="../results/Figure3C"
mkdir -p $out_dir

# # Pipeline for TF enrichment V0
# bed_dir=../data/TF_binding_site_enrichment/bed_files

# bed_list=($(ls $bed_dir/*.bed))
# for bed_file in "${bed_list[@]}"; do
#     # compute background region
#     bash ../hiplex_paper/TF_binding_site_enrichment/compute_background_refine.sh $bed_file $data_out_dir $cisgenome_bin
# done

# # merge background regions
# Rscript ../hiplex_paper/TF_binding_site_enrichment/merge_ctrl_seq.R \
#     ${data_out_dir} \
#     ${data_out_dir}

# # calculate TF enrichment
# for bed_file in "${bed_list[@]}"; do
#     echo $bed_file
#     file_name=$(basename "${bed_file%.*}")
#     control_bed_file=${data_out_dir}/ctrl.bed
#     mkdir -p ${data_out_dir}/${file_name}/peak_result/
#     peak_out_dir=${data_out_dir}/${file_name}/peak_result/result.tsv

#     Rscript ../hiplex_paper/TF_binding_site_enrichment/peak_analysis.R $bed_file $control_bed_file $peak_out_dir 1
# done

# # rm -rf ../data/TF_binding_site_enrichment/gw/*/peak_result/

# plot
Rscript ../hiplex_paper/TF_binding_site_enrichment/peak_heatmap_bash.R \
    ${data_out_dir} \
    ${out_dir}


