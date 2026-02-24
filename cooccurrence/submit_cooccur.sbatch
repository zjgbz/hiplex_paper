#!/bin/bash

#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10:00:00
#SBATCH --mem=20GB
#SBATCH --output=/projects/foundation_model_for_single_cell_multiomics_data/output_error/hiplex/%x_%A_%a.out
#SBATCH --error=/projects/foundation_model_for_single_cell_multiomics_data/output_error/hiplex/%x_%A_%a.err
#SBATCH --mail-type=END,ALL
#SBATCH --mail-user=minzhi.hpc.status@gmail.com
##SBATCH --array=312

module load anaconda3/2023.09

target_pair_1=$(awk -vawk_variable="${SLURM_ARRAY_TASK_ID}" 'NR==awk_variable {print $1}' ${target_pair_avail_dir_filename})
target_pair_2=$(awk -vawk_variable="${SLURM_ARRAY_TASK_ID}" 'NR==awk_variable {print $2}' ${target_pair_avail_dir_filename})

cd ${work_dir}

python ${script_dir}/cooccur.py ${tissue_name} ${frag_type} ${bin_size} ${process} ${region_selection} ${cell_list_dir_filename} ${cell_reads_thres} ${region_cell_dir} ${perm_num} ${scen} ${row_cluster_sel} ${target_pair_1} ${target_pair_2} ${test_result_dir} ${region_dir} ${figure_dir} > result_${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out 2>&1

exit
