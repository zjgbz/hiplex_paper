#!/bin/bash

frag_type_list=( "mixed" )
bin_size_list=( "800" )
# bin_size_list=( "promoter_-200-200" )
# target_qc_type_list=( "orig" "hm-only" )
# target_qc_type_list=( "hm-writer" "writer-sptf" "hm-only" "hm-tf" )
target_qc_type_list=( "all-qc" )
# post_process_list=( "log10" "log10_noallzero" )
# post_process_list=( "libnorm_noAllZero_sqrt" "libnorm_noAllZero_log2" )
# post_process_list=( "libnorm_noAllZero_sqrt_qnorm" )
# post_process_list=( "libnorm_noAllZero_log2_qnorm" "libnorm_noAllZero_log2" "libnorm_noAllZero" )
post_process_list=( "libnorm_noAllZero_log2_qnorm" )
scen_list=( "V1V2" "T1T2" "V1V2-T1T2" )
# scen_list=( "V1V2-T1T2" )

for target_qc_type in "${target_qc_type_list[@]}"
do
	for frag_type in "${frag_type_list[@]}"
	do
		for bin_size in "${bin_size_list[@]}"
		do
			for post_process in "${post_process_list[@]}"
			do
				for scen in "${scen_list[@]}"
				do
					sbatch --job-name=mav_screen_${scen}_${frag_type}_${bin_size}_${target_qc_type}_${post_process} --export=frag_type=${frag_type},bin_size=${bin_size},qc_type=${target_qc_type},post_process=${post_process},scen=${scen} submit_mav_screen.sbatch
				done
			done
		done
	done
done