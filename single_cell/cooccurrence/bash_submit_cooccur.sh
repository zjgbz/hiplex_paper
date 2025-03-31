#!/bin/bash

tissue_name="K562"
frag_type="mixed"
bin_size="800"
process="orig"
region_selection="bulk_cluster_rows"
perm_num="10000"
# scen_list=( "region" "cell" )
scen_list=( "region" )
row_cluster_sel_list=( "all" )
# row_cluster_sel_list=( "all" {A..O} )
coexp_num_thres="1"
cell_reads_thres_list=( "0" "1" )

script_dir="/projects/foundation_model_for_single_cell_multiomics_data/hiplex/script/cluster"
sc_root_dir="/projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell"
wgc_root_dir="${sc_root_dir}/wgc"
region_cell_dir="${wgc_root_dir}/${tissue_name}/${frag_type}/region_cell/${bin_size}"
work_dir="${wgc_root_dir}/${tissue_name}/${frag_type}/cooccur"
test_result_root_dir="${work_dir}/test_result"
region_root_dir="${work_dir}/region"
figure_root_dir="${work_dir}/figure"
observation_root_dir="${work_dir}/observation"
target_pair_avail_root_dir="${observation_root_dir}/target_pair_avail"
cell_list_dir="${sc_root_dir}/${tissue_name}/cell_list"

for cell_reads_thres in "${cell_reads_thres_list[@]}"
do
	cell_list_filename="cell_list_36x_qc-${cell_reads_thres}.txt"
	cell_list_dir_filename="${cell_list_dir}/${cell_list_filename}"
	for scen in "${scen_list[@]}"
	do
		test_result_scen_root_dir="${test_result_root_dir}/${scen}"
		region_scen_root_dir="${region_root_dir}/${scen}"
		figure_scen_root_dir="${figure_root_dir}/${scen}"
		mkdir -p ${test_result_scen_root_dir}
		mkdir -p ${region_scen_root_dir}
		mkdir -p ${figure_scen_root_dir}
		target_pair_avail_dir="${target_pair_avail_root_dir}/${scen}"
		for row_cluster_sel in "${row_cluster_sel_list[@]}"
		do
			test_result_dir="${test_result_scen_root_dir}/${row_cluster_sel}"
			region_dir="${region_scen_root_dir}/${row_cluster_sel}"
			figure_dir="${figure_scen_root_dir}/${row_cluster_sel}"
			mkdir -p ${test_result_dir}
			mkdir -p ${region_dir}
			mkdir -p ${figure_dir}

			target_pair_avail_filename="target_pair_avail_cell_reads_thres-${cell_reads_thres}_${scen}_${row_cluster_sel}_coexp_num_thres-${coexp_num_thres}.tsv"
			target_pair_avail_dir_filename="${target_pair_avail_dir}/${target_pair_avail_filename}"
			target_pair_avail_num=$(wc -l < ${target_pair_avail_dir_filename})
			# job_array_range="2-${target_pair_avail_num}"
			job_array_range="2-9,199-207,305-312"
			sbatch --job-name=cooccur_${scen}_${row_cluster_sel} \
					--array=${job_array_range} \
					--export=tissue_name=${tissue_name},frag_type=${frag_type},bin_size=${bin_size},process=${process},region_selection=${region_selection},cell_list_dir_filename=${cell_list_dir_filename},cell_reads_thres=${cell_reads_thres},region_cell_dir=${region_cell_dir},perm_num=${perm_num},scen=${scen},row_cluster_sel=${row_cluster_sel},test_result_dir=${test_result_dir},region_dir=${region_dir},figure_dir=${figure_dir},work_dir=${work_dir},script_dir=${script_dir},target_pair_avail_dir_filename=${target_pair_avail_dir_filename} \
					submit_cooccur.sbatch
		done
	done
done