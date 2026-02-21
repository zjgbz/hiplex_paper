#!/bin/bash

##SBATCH -p shared
#SBATCH -p cancergen
##SBATCH -p hongkai
##SBATCH -p pi-sysadmin
##SBATCH -p gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=/dcl02/hongkai/data/mjiang/output_error/multitag/%x_%A_%a.out
#SBATCH --error=/dcl02/hongkai/data/mjiang/output_error/multitag/%x_%A_%a.err
#SBATCH --mail-type=END,ALL
#SBATCH --mail-user=minzhi.hpc.status@gmail.com
#SBATCH --array=2-5

module load R/4.3

scen_list=( "V" "T" "T1" "T2" "V1" "V2" )
# scen_list=( "V" "T" "NaB_36x_CT_T1_S1" "NaB_36x_CT_T2_S2" "NaB_36x_CT_V1_S3" "NaB_36x_CT_V2_S4" )
# scen_list=( "pseudobulk" )

scen_raw=${scen_list[${SLURM_ARRAY_TASK_ID}]}
# if [[ "$frag_type" == "mixed" ]]; then
# 	scen=${scen_raw}
# else
# 	scen=${scen_raw}_${frag_type}
# fi
scen=${scen_raw}_${frag_type}
work_dir=${work_root_dir}/${frag_type}/${bin_size}

mkdir -p ${work_root_dir}
mkdir -p ${work_root_dir}/${frag_type}
mkdir -p ${work_dir}

cd ${work_dir}

Rscript ${script_dir}/wgc_formal.R ${scen} ${frag_type} ${bin_size} ${target_qc_type} ${libnorm_type} ${work_dir} ${align_dir} > result_${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out 2>&1

exit