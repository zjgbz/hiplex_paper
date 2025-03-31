#!/bin/bash

tissue_name="K562"
frag_type="mixed"
bin_size="800"
process="orig"
region_selection="bulk_cluster_rows"
scen_list=( "region" "cell" )
# scen_list=( "cell" )
coexp_num_thres="1"
cell_reads_thres="0"

script_dir="/projects/foundation_model_for_single_cell_multiomics_data/hiplex/script/cluster"
sc_root_dir="/projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell"
wgc_root_dir="${sc_root_dir}/wgc"
cell_list_dir="${sc_root_dir}/${tissue_name}/cell_list"
cell_list_filename="cell_list_36x_qc-${cell_reads_thres}.txt"
cell_list_dir_filename="${cell_list_dir}/${cell_list_filename}"
region_cell_dir="${wgc_root_dir}/${tissue_name}/${frag_type}/region_cell/${bin_size}"
work_dir="${wgc_root_dir}/${tissue_name}/${frag_type}/cooccur/observation"
observation_root_dir="${work_dir}"
observation_all_root_dir="${observation_root_dir}/target_pair_all"
observation_avail_root_dir="${observation_root_dir}/target_pair_avail"
observation_summary_root_dir="${observation_root_dir}/target_pair_summary"

for scen in "${scen_list[@]}"
do
    observation_all_dir="${observation_all_root_dir}/${scen}"
    observation_avail_dir="${observation_avail_root_dir}/${scen}"
    observation_summary_dir="${observation_summary_root_dir}/${scen}"
    mkdir -p ${observation_all_dir}
    mkdir -p ${observation_avail_dir}
    mkdir -p ${observation_summary_dir}
    
    sbatch --job-name=cooccur_observation_${scen} \
    		--export=tissue_name=${tissue_name},frag_type=${frag_type},bin_size=${bin_size},process=${process},region_selection=${region_selection},cell_list_dir_filename=${cell_list_dir_filename},cell_reads_thres=${cell_reads_thres},region_cell_dir=${region_cell_dir},scen=${scen},coexp_num_thres=${coexp_num_thres},observation_all_dir=${observation_all_dir},observation_avail_dir=${observation_avail_dir},observation_summary_dir=${observation_summary_dir},work_dir=${work_dir},script_dir=${script_dir} \
    		submit_cooccur_observation.sbatch
done