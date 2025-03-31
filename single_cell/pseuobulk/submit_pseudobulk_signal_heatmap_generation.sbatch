#!/bin/bash

#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3-00:00:00
#SBATCH --mem=150GB
#SBATCH --output=/projects/foundation_model_for_single_cell_multiomics_data/output_error/hiplex/%x_%A_%a.out
#SBATCH --error=/projects/foundation_model_for_single_cell_multiomics_data/output_error/hiplex/%x_%A_%a.err
#SBATCH --mail-type=END,ALL
#SBATCH --mail-user=minzhi.hpc.status@gmail.com
#SBATCH --array=1-24

module load R/4.4.0+Bioconductor
module load gcc/13.1.0

script_dir="/projects/foundation_model_for_single_cell_multiomics_data/hiplex/script/signal_visual"
work_dir="/projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell/wgc/K562/mixed"
cd ${wgc_dir}

Rscript ${script_dir}/pseudobulk_signal_heatmap_generation.R ${SLURM_ARRAY_TASK_ID} ${target_pair_selected} > result_${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out 2>&1
# Rscript ${script_dir}/pseudobulk_signal_heatmap_generation.R ${SLURM_ARRAY_TASK_ID} ${target_pair_selected}

exit