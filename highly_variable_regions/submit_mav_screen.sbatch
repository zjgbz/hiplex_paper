#!/bin/bash

##SBATCH -p sas
#SBATCH -p cancergen
##SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=03:00:00
#SBATCH --mem=500GB
#SBATCH --output=/dcl02/hongkai/data/mjiang/output_error/multitag/%x_%A_%a.out
#SBATCH --error=/dcl02/hongkai/data/mjiang/output_error/multitag/%x_%A_%a.err
#SBATCH --mail-type=END,ALL
#SBATCH --mail-user=minzhi.hpc.status@gmail.com
##SBATCH --array=0-59
#SBATCH --array=0-20

module load R/4.3

# fitting_model="loess"
fitting_model="gam"

if [ $fitting_model == "loess" ]
then
	para_list=( "0.041" "0.042" "0.043" "0.044" "0.045" "0.046" "0.047" "0.048" "0.049" )
elif [ $fitting_model == "gam" ]
then
	# para_list=( "51" "52" "53" "54" "55" "56" "57" "58" "59" "60" "61" "62" "63" "64" "65" "66" "67" "68" "69" "70" "71" "72" "73" "74" "75" "76" "77" "78" "79" "80" "81" "82" "83" "84" "85" "86" "87" "88" "89" "90" "91" "92" "93" "94" "95" "96" "97" "98" "99" "100" "101" "102" "103" "104" "105" "106" "107" "108" "109" "110" )
	para_list=( "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "111" "112" "113" "114" "115" "116" "117" "118" "119" "120" )
fi

robust_check="21"
parameter=${para_list[${SLURM_ARRAY_TASK_ID}]}
seed="42"
nrow_sample="300000"

cd /dcs05/hongkai/data/next_cutntag/bulk/wgc/${frag_type}/${bin_size}

Rscript /dcs05/hongkai/data/next_cutntag/script/cluster/mav_screen.R ${fitting_model} ${robust_check} ${parameter} ${seed} ${nrow_sample} ${frag_type} ${bin_size} ${scen} ${qc_type} ${post_process} > result_${SLURM_JOB_NAME}_${fitting_model}_${robust_check}_${parameter}_${seed}_${nrow_sample}_${SLURM_ARRAY_JOB_ID}.out 2>&1

exit