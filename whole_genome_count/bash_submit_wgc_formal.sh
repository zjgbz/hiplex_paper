#!/bin/bash

hpc="jhpce"
# frag_type_list=( "mixed" "sub" "mono" "di" )
# frag_type_list=( "sub" "mono" "di" )
frag_type_list=( "mixed" )
bin_size_list=( "800" )
# bin_size_list=( "200" "300" "500" )
target_qc_type_list=( "all-qc" )
# target_qc_type_list=( "all" "hm-only" "hm-writer" )
libnorm_type_list=( "libnorm-mean" "libnorm-median" )

if [[ "${hpc}" == "rockfish" ]]; then
	align_dir="/data/hji7/minzhi/multitag/data_align"
	work_root_dir="/data/hji7/minzhi/multitag/wgc"
	script_dir="/data/hji7/minzhi/multitag/multitag_script"
elif [[ "${hpc}" == "jhpce" ]]; then
	# align_dir="/dcs05/hongkai/data/next_cutntag/bulk/homotone_heterotone_merged/data_align"
	align_dir="/dcs05/hongkai/data/next_cutntag/bulk/homotone_heterotone_merged/data_align/frag_decon/valley-all-qc"
	work_root_dir="/dcs05/hongkai/data/next_cutntag/bulk/wgc"
	script_dir="/dcs05/hongkai/data/next_cutntag/script/cluster"
fi

for target_qc_type in "${target_qc_type_list[@]}"
do
	for frag_type in "${frag_type_list[@]}"
	do
		for bin_size in "${bin_size_list[@]}"
		do
			for libnorm_type in "${libnorm_type_list[@]}"
			do
				if [[ ${bin_size} == "200" ]]; then
					mem="400GB"
				elif [[ ${bin_size} == "300" ]]; then
					mem="250GB"
				elif [[ ${bin_size} == "500" ]]; then
					mem="180GB"
				elif [[ ${bin_size} == "800" ]] && [[ ${target_qc_type} == "all" ]]; then
					mem="900GB" # 600G for all; 
					walltime="3-00:00:00" # 12h for all;
				elif [[ ${bin_size} == "800" ]] && [[ ${target_qc_type} == "all-qc" ]]; then
					mem="600GB" # 600G for all; 
					walltime="1-00:00:00" # 12h for all;
				elif [[ ${bin_size} == "800" ]] && [[ ${target_qc_type} == "hm-only" ]]; then
					mem="60GB"
					walltime="04:00:00"
				elif [[ ${bin_size} == "800" ]] && [[ ${target_qc_type} == "hm-writer" ]]; then
					mem="180GB"
					walltime="08:00:00"
				elif [[ ${bin_size} == "3200" ]]; then
					mem="50GB"
				fi
				sbatch --job-name=wgc_formal_${frag_type}_${bin_size}_${target_qc_type} --mem=${mem} --time=${walltime} --export=frag_type=${frag_type},bin_size=${bin_size},target_qc_type=${target_qc_type},libnorm_type=${libnorm_type},align_dir=${align_dir},work_root_dir=${work_root_dir},script_dir=${script_dir} submit_wgc_formal.sbatch
			done
		done
	done
done