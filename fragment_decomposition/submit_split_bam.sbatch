#!/bin/bash

#SBATCH --job-name=split_bam
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=20:00:00
#SBATCH --mem=5GB
#SBATCH --output=/dcl02/hongkai/data/mjiang/output_error/multitag/%x_%A_%a.out
#SBATCH --error=/dcl02/hongkai/data/mjiang/output_error/multitag/%x_%A_%a.err
#SBATCH --mail-type=END,ALL
#SBATCH --mail-user=minzhi.hpc.status@gmail.com
##SBATCH --mail-user=mjiang26@jhu.edu
#SBATCH --array=0-1

module load samtools

# thres_name="valley-all-qc"
thres_name="valley-V-qc"
# scen_list=( "V" "T" "TV" )
# scen_list=( "NaB_36x_CT_T1_S1" "NaB_36x_CT_T2_S2" "NaB_36x_CT_V1_S3" "NaB_36x_CT_V2_S4" )
# scen_simple_name_list=( "T1" "T2" "V1" "V2" )
scen_list=( "T" "V" )
scen_simple_name_list=( "T" "V" )
scen=${scen_list[${SLURM_ARRAY_TASK_ID}]}
scen_simple_name=${scen_simple_name_list[${SLURM_ARRAY_TASK_ID}]}

target_list=( "IgG_control" "H3K36me3" "H3K4me1" "H3K27ac" "H3S10ph" "H2A_XS139ph" "H3K79me3" "H3K9me2" "H3K9me3" "H3K14ac" "H3K27me3" "H3K4me3" "SETD2" "MLL4_MLL2_KMT2B" "CBP_CREBBP" "EP300" "MSK1" "MSK2" "PIM1" "CDK8" "AURORA_Aurora_B" "EHMT2" "SuVar39_SUV39H1" "EHMT1" "EZH2" "MLL1_KMT2A" "CTCF" "POLR2AphosphoS2" "cJun" "cFos" "Max" "Myc" "USF1" "USF2" "NRF1" "YY1" "H3K9ac" "unknown" )
load_root_dir="/dcs05/hongkai/data/next_cutntag/bulk/homotone_heterotone_merged/data_align"
save_root_raw_dir="/dcs05/hongkai/data/next_cutntag/bulk/homotone_heterotone_merged/data_align/frag_decon"
save_root_dir=${save_root_raw_dir}/${thres_name}
mkdir -p ${save_root_dir}

load_dir=${load_root_dir}/${scen}
mixed_dir=${save_root_dir}/${scen_simple_name}_mixed
subneucleo_dir=${save_root_dir}/${scen_simple_name}_sub
mononeucleo_dir=${save_root_dir}/${scen_simple_name}_mono
dineucleo_dir=${save_root_dir}/${scen_simple_name}_di
# trineucleo_dir=${save_root_dir}/${scen}_tri

mkdir -p ${mixed_dir}
cp -r ${load_dir}/* ${mixed_dir}/
mkdir -p ${subneucleo_dir}
mkdir -p ${mononeucleo_dir}
mkdir -p ${dineucleo_dir}
# mkdir -p ${trineucleo_dir}

for target_1 in "${target_list[@]}"
do
	for target_2 in "${target_list[@]}"
	do
		if [ "${target_1}" \< "${target_2}" ]
		then
			mixed_filename="${target_1}-${target_2}.bam"
		elif [ "${target_1}" \> "${target_2}" ]
		then
			echo "${target_1} is larger than ${target_2}!"
		elif [ "${target_1}" = "${target_2}" ]
		then
			mixed_filename="${target_1}-${target_2}.bam"
		fi
		mixed_dir_filename=${load_dir}/${mixed_filename}

		# samtools view -h ${mixed_dir_filename} | awk 'substr($0,1,1)=="@" || ($9>=24 && $9<=120) || ($9<=-24 && $9>=-120)' | samtools view -b > ${subneucleo_dir}/${mixed_filename}
		# samtools view -h ${mixed_dir_filename} | awk 'substr($0,1,1)=="@" || ($9>=121 && $9<=298) || ($9<=-121 && $9>=-298)' | samtools view -b > ${mononeucleo_dir}/${mixed_filename}
		# samtools view -h ${mixed_dir_filename} | awk 'substr($0,1,1)=="@" || ($9>=299 && $9<=453) || ($9<=-299 && $9>=-453)' | samtools view -b > ${dineucleo_dir}/${mixed_filename}
		# samtools view -h ${mixed_dir_filename} | awk 'substr($0,1,1)=="@" || ($9>=454 && $9<=800) || ($9<=-454 && $9>=-800)' | samtools view -b > ${trineucleo_dir}/${mixed_filename}
		# samtools view -h ${mixed_dir_filename} | awk 'substr($0,1,1)=="@" || ($9>=24 && $9<=128) || ($9<=-24 && $9>=-128)' | samtools view -b > ${subneucleo_dir}/${mixed_filename}
		# samtools view -h ${mixed_dir_filename} | awk 'substr($0,1,1)=="@" || ($9>=129 && $9<=295) || ($9<=-129 && $9>=-295)' | samtools view -b > ${mononeucleo_dir}/${mixed_filename}
		# samtools view -h ${mixed_dir_filename} | awk 'substr($0,1,1)=="@" || ($9>=296 && $9<=800) || ($9<=-296 && $9>=-800)' | samtools view -b > ${dineucleo_dir}/${mixed_filename}
		samtools view -h ${mixed_dir_filename} | awk 'substr($0,1,1)=="@" || ($9>=24 && $9<=126) || ($9<=-24 && $9>=-126)' | samtools view -b > ${subneucleo_dir}/${mixed_filename}
		samtools view -h ${mixed_dir_filename} | awk 'substr($0,1,1)=="@" || ($9>=127 && $9<=297) || ($9<=-127 && $9>=-297)' | samtools view -b > ${mononeucleo_dir}/${mixed_filename}
		samtools view -h ${mixed_dir_filename} | awk 'substr($0,1,1)=="@" || ($9>=298 && $9<=800) || ($9<=-298 && $9>=-800)' | samtools view -b > ${dineucleo_dir}/${mixed_filename}
	done
done

exit
