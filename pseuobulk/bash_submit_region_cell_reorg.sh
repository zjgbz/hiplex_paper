#!/bin/bash

# process_list=( "libnorm" "orig" "libnorm_log2" )
process_list=( "orig" )
# tissue_name_list=( "brain" "iPSC" )
tissue_name_list=( "K562" )
# target_qc_type_list=( "all-36x-qc-3000" "all-36x-qc-3000" "all-36x-qc-3000" )
target_qc_type_list=( "all-17x" )
# target_qc_type_list=( "all" "all" )
# bin_size_list=( "promoter_-2000-2000" )
# bin_size_list=( "5000" "10000" "50000" "100000" )
# bin_size_list=( "BIRD_1000" "brain_iPSC_H3K27me3-H3K27me3_peak_combine" "brain_iPSC_H3K4me3-H3K4me3_peak_combine" )
# bin_size_list=( "brain_iPSC_H3K27me3-H3K4me3_peak_combine" )
# bin_size_list=( "brain_iPSC_tpc-01_peak_union" )
# bin_size_list=( "K562_tpc-01_peak_union" )
# bin_size_list=( "K562_H3K4me3-H3K4me3_peak" )
# bin_size_list=( "K562_H3K27me3-H3K4me3_peak" )
bin_size_list=( "K562_H3K27me3-H3K27me3_peak" )
# bin_size_list=( "800" )
script_dir="/projects/foundation_model_for_single_cell_multiomics_data/hiplex/script/cluster"
sc_root_dir="/projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell"
target_reads_thres="3000"

data_summary_dir="/projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell/data_summary"
target_pair_list_filename="sc_brain_iPSC_36x_target_pair_qc-${target_reads_thres}.txt"
# target_pair_list_filename="sc_brain_iPSC_K562_17x_target_pair_qc-${target_reads_thres}.txt"
# target_pair_list_filename="sc_K562_17x_target_pair_top_3.txt"
# target_pair_list_filename="K562_bulk_sc_qc_target_pairs.txt"
target_pair_list_dir_filename="${data_summary_dir}/${target_pair_list_filename}"
target_pair_num=$(wc -l < ${target_pair_list_dir_filename})
# job_array_range="2-${target_pair_num}"
# job_array_range="18,20,28"
job_array_range="18"
echo "target_pair_num: ${target_pair_num}, job_array_range: ${job_array_range}"

# for i in {0..1}
for i in {0..0}
do
    tissue_name=${tissue_name_list[$i]}
	target_qc_type=${target_qc_type_list[$i]}
	echo "$i, $tissue_name, $target_qc_type"

	for bin_size in "${bin_size_list[@]}"
	do
		for process in "${process_list[@]}"
		do
			sbatch --job-name=region_cell_reorg_${bin_size}_${tissue_name}_${process} --array=${job_array_range} --export=process=${process},bin_size=${bin_size},tissue_name=${tissue_name},target_qc_type=${target_qc_type},target_reads_thres=${target_reads_thres},sc_root_dir=${sc_root_dir},script_dir=${script_dir},target_pair_list_dir_filename=${target_pair_list_dir_filename} submit_region_cell_reorg.sbatch
		done
	done
done

# for tissue_name in "${tissue_name_list[@]}"
# do
# 	for bin_size in "${bin_size_list[@]}"
# 	do
# 		for process in "${process_list[@]}"
# 		do
# 			sbatch --job-name=region_cell_reorg_${tissue_name}_${process} --array=${job_array_range} --export=process=${process},bin_size=${bin_size},tissue_name=${tissue_name},target_reads_thres=${target_reads_thres},sc_root_dir=${sc_root_dir},script_dir=${script_dir},target_pair_list_dir_filename=${target_pair_list_dir_filename} submit_region_cell_reorg.sbatch
# 		done
# 	done
# done