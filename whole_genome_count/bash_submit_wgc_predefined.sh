#!/bin/bash

hpc="jhpce"
# frag_type_list=( "mixed" "sub" "di" "mono" "tri" )
frag_type_list=( "mixed" )
# bin_size_list=( "tzo-promoter_-1000-1000" "tzos-promoter_-1000-1000" "vzo-promoter_-1000-1000" "vzos-promoter_-1000-1000" )
# bin_size_list=( "promoter_-2000-0" "promoter_-2000-1000" "promoter_-2000-500" "promoter_-5000-0" "promoter_-5000-1000" "promoter_-5000-500" "naz-promoter_-2000-0" "naz-promoter_-2000-1000" "naz-promoter_-2000-500" "naz-promoter_-5000-0" "naz-promoter_-5000-1000" "naz-promoter_-5000-500" )
# bin_size_list=( "nas1-promoter_-1000-0" "nas1-promoter_-1000-1000" "nas1-promoter_-1000-500" "nas1-promoter_-2000-0" "nas1-promoter_-2000-1000" "nas1-promoter_-2000-500" "nas1-promoter_-5000-0" "nas1-promoter_-5000-1000" "nas1-promoter_-5000-500" )
# bin_size_list=( "promoter_-1000-1000" "promoter_-1000-0" "promoter_-1000-200" "promoter_-1000-500" "promoter_-2000-0" "promoter_-200-0" "promoter_-500-0" "naz-promoter_-1000-0" "naz-promoter_-1000-1000" "naz-promoter_-1000-200" "naz-promoter_-1000-500" "naz-promoter_-2000-0" "naz-promoter_-200-0" "naz-promoter_-500-0" )
# bin_size_list=( "promoter_-101000-101000" "promoter_-201000-201000" "promoter_-501000-501000" "promoter_-1001000-1001000" )
# bin_size_list=( "promoter_-50-300" "genebody_0-300" )
# bin_size_list=( "genebody_-300-3000" )
bin_size_list=( "promoter_-5000-5000" )
target_qc_type_list=( "all-qc" )

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
			if [[ ${bin_size} == "200" ]]; then
				mem="400GB"
			elif [[ ${bin_size} == "300" ]]; then
				mem="250GB"
			elif [[ ${bin_size} == "500" ]]; then
				mem="180GB"
			elif [[ ${bin_size} == "800" ]]; then
				mem="8GB"
			elif [[ ${bin_size} == "promoter_-1000-1000" ]]; then
				mem="10GB" # 1h
			elif [[ ${bin_size} == "ELS" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "genebody" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "genebody_0-300" ]]; then
				mem="100GB"
			elif [[ ${bin_size} == "genebody_-300-3000" ]]; then
				mem="200GB"
			elif [[ ${bin_size} == "promoter_-50-300" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "promoter_-1000-0" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "promoter_-1000-200" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "promoter_-2000-0" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "promoter_-200-0" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "promoter_-500-0" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "promoter_-1000-1000" ]]; then
				mem="8GB"
			elif [[ ${bin_size} == "promoter_-1000-500" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "naz-promoter_-1000-0" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "naz-promoter_-1000-200" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "naz-promoter_-2000-0" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "naz-promoter_-200-0" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "naz-promoter_-500-0" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "naz-promoter_-1000-1000" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "naz-promoter_-1000-500" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "tzo-promoter_-1000-1000" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "tzos-promoter_-1000-1000" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "vzo-promoter_-1000-1000" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "vzos-promoter_-1000-1000" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "naz-promoter_-5000-0" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "naz-promoter_-5000-1000" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "naz-promoter_-5000-500" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "naz-promoter_-2000-0" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "naz-promoter_-2000-1000" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "naz-promoter_-2000-500" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "promoter_-5000-0" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "promoter_-5000-1000" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "promoter_-5000-5000" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "promoter_-5000-500" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "promoter_-2000-0" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "promoter_-2000-1000" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "promoter_-2000-500" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "nas1-promoter_-1000-0" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "nas1-promoter_-1000-1000" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "nas1-promoter_-1000-500" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "nas1-promoter_-5000-0" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "nas1-promoter_-5000-1000" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "nas1-promoter_-5000-500" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "nas1-promoter_-2000-0" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "nas1-promoter_-2000-1000" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "nas1-promoter_-2000-500" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "promoter_-51000-51000" ]]; then
				mem="10GB"
			elif [[ ${bin_size} == "all_peak_set" ]]; then
				mem="50GB"
			elif [[ ${bin_size} == "promoter_-101000-101000" ]]; then
				mem="50GB"
			elif [[ ${bin_size} == "promoter_-201000-201000" ]]; then
				mem="50GB"
			elif [[ ${bin_size} == "promoter_-501000-501000" ]]; then
				mem="50GB"
			elif [[ ${bin_size} == "promoter_-1001000-1001000" ]]; then
				mem="50GB"
			elif [[ ${bin_size} == "promoter_-2001000-2001000" ]]; then
				mem="50GB"
			elif [[ ${bin_size} == "promoter_-5001000-5001000" ]]; then
				mem="50GB"
			elif [[ ${bin_size} == "promoter_-10001000-10001000" ]]; then
				mem="50GB"
			fi
			sbatch --job-name=wgc_predefined_${frag_type}_${bin_size}_${target_qc_type} --mem=${mem} --export=frag_type=${frag_type},bin_size=${bin_size},target_qc_type=${target_qc_type},align_dir=${align_dir},work_root_dir=${work_root_dir},script_dir=${script_dir} submit_wgc_predefined.sbatch
		done
	done
done