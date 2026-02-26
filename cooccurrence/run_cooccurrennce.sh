#!/bin/bash

# module load conda_R  # 本地不需要
# module load anaconda3/2023.09  # 本地不需要，确保本地有 python 环境

cell_reads_thres_list=( "0" "1" )

test_result_root_dir="/projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell/wgc/K562/mixed/cooccur/test_result"
region_root_dir="/projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell/wgc/K562/mixed/cooccur/region"
figure_root_dir="/projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell/wgc/K562/mixed/cooccur/figure"
target_pair_avail_root_dir="/projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell/wgc/K562/mixed/cooccur/observation/target_pair_avail"
cell_list_dir="/projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell/K562/cell_list"

# test_result_root_dir="/projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell/wgc/K562/mixed/cooccur/test_result"
# region_root_dir="/projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell/wgc/K562/mixed/cooccur/region"
# figure_root_dir="/projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell/wgc/K562/mixed/cooccur/figure"
# target_pair_avail_root_dir="/projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell/wgc/K562/mixed/cooccur/observation/target_pair_avail"
# cell_list_dir="/projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell/K562/cell_list"

for cell_reads_thres in "${cell_reads_thres_list[@]}"
do
	cell_list_filename="cell_list_36x_qc-${cell_reads_thres}.txt"
	cell_list_dir_filename="${cell_list_dir}/${cell_list_filename}"

	test_result_scen_root_dir="${test_result_root_dir}/region"
	region_scen_root_dir="${region_root_dir}/region"
	figure_scen_root_dir="${figure_root_dir}/region"
	mkdir -p ${test_result_scen_root_dir}
	mkdir -p ${region_scen_root_dir}
	mkdir -p ${figure_scen_root_dir}
	target_pair_avail_dir="${target_pair_avail_root_dir}/region"

	test_result_dir="${test_result_scen_root_dir}/all"
	region_dir="${region_scen_root_dir}/all"
	figure_dir="${figure_scen_root_dir}/all"
	mkdir -p ${test_result_dir}
	mkdir -p ${region_dir}
	mkdir -p ${figure_dir}

	target_pair_avail_filename="target_pair_avail_cell_reads_thres-${cell_reads_thres}_region_all_coexp_num_thres-1.tsv"
	target_pair_avail_dir_filename="${target_pair_avail_dir}/${target_pair_avail_filename}"

	# 本地运行：用 for 循环代替 sbatch --array=2-9,199-207,305-312
	for task_id in $(seq 2 9) $(seq 199 207) $(seq 305 312)
	do
		echo "Running task_id=${task_id}, cell_reads_thres=${cell_reads_thres}"

		target_pair_1=$(awk -v awk_variable="${task_id}" 'NR==awk_variable {print $1}' ${target_pair_avail_dir_filename})
		target_pair_2=$(awk -v awk_variable="${task_id}" 'NR==awk_variable {print $2}' ${target_pair_avail_dir_filename})

		cd /projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell/wgc/K562/mixed/cooccur

		python /projects/foundation_model_for_single_cell_multiomics_data/hiplex/script/cluster/cooccur.py \
			K562 \
			mixed \
			800 \
			orig \
			bulk_cluster_rows \
			${cell_list_dir_filename} \
			${cell_reads_thres} \
			/projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell/wgc/K562/mixed/region_cell/800 \
			10000 \
			region \
			all \
			${target_pair_1} \
			${target_pair_2} \
			${test_result_dir} \
			${region_dir} \
			${figure_dir} \
			> result_cooccur_region_all_${cell_reads_thres}_${task_id}.out 2>&1

		echo "Finished task_id=${task_id}"
	done
done