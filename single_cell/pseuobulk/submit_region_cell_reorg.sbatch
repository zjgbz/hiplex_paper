#!/bin/bash

#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=500GB
#SBATCH --output=/projects/foundation_model_for_single_cell_multiomics_data/output_error/hiplex/%x_%A_%a.out
#SBATCH --error=/projects/foundation_model_for_single_cell_multiomics_data/output_error/hiplex/%x_%A_%a.err
#SBATCH --mail-type=END,ALL
#SBATCH --mail-user=minzhi.hpc.status@gmail.com

# source activate /projects/foundation_model_for_single_cell_multiomics_data/software/env/py_310
# echo $(which python)
module load anaconda3/2023.09

target_pair_name=$(awk -vawk_variable="${SLURM_ARRAY_TASK_ID}" 'NR==awk_variable' ${target_pair_list_dir_filename})

mkdir -p ${sc_root_dir}/wgc/${tissue_name}/mixed/region_cell/${bin_size}
cd ${sc_root_dir}/wgc/${tissue_name}/mixed/region_cell/${bin_size}

python ${script_dir}/region_cell_reorg.py ${tissue_name} ${target_pair_name} ${process} ${sc_root_dir} ${target_reads_thres} ${bin_size} ${target_qc_type} > result_${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out 2>&1

exit