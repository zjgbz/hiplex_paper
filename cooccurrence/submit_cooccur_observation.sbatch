#!/bin/bash

#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=05:00:00
#SBATCH --mem=30GB
#SBATCH --output=/projects/foundation_model_for_single_cell_multiomics_data/output_error/hiplex/%x_%A_%a.out
#SBATCH --error=/projects/foundation_model_for_single_cell_multiomics_data/output_error/hiplex/%x_%A_%a.err
#SBATCH --mail-type=END,ALL
#SBATCH --mail-user=minzhi.hpc.status@gmail.com
#SBATCH --array=0-15

module load anaconda3/2023.09

row_cluster_sel_list=( "all" {A..O} )
row_cluster_sel=${row_cluster_sel_list[${SLURM_ARRAY_TASK_ID}]}

cd ${work_dir}

python ${script_dir}/cooccur_observation.py ${tissue_name} ${frag_type} ${bin_size} ${process} ${region_selection} ${cell_list_dir_filename} ${cell_reads_thres} ${region_cell_dir} ${scen} ${row_cluster_sel} ${coexp_num_thres} ${observation_all_dir} ${observation_avail_dir} ${observation_summary_dir} > result_${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out 2>&1

exit
