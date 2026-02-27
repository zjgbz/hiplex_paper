#!/bin/bash
#SBATCH --job-name=hiplex
# #SBATCH --partition=gpu     
# #SBATCH --gres=gpu:1
#SBATCH --time=1-00:00:00
#SBATCH --mem=512G
# #SBATCH --cpus-per-task=8               
#SBATCH --output=/dcs07/hongkai/data/mjiang/hiplex/paper_materials/hiplex_paper/submit_all_steps.out
#SBATCH --error=/dcs07/hongkai/data/mjiang/hiplex/paper_materials/hiplex_paper/submit_all_steps.err
#SBATCH --mail-user=yhu157@jh.edu
#SBATCH --mail-type=END,FAIL

cd /dcs07/hongkai/data/mjiang/hiplex/paper_materials/hiplex_paper

module load conda_R
bash /dcs07/hongkai/data/mjiang/hiplex/paper_materials/hiplex_paper/peak_annotation/run_peak_annotation.sh
# Rscript ../hiplex_paper/peak_annotation/peakAnnotationChIPseekerCCREForLoop_different_peak_type.R \
#     ../data/peak/data_peak_auc_0.05_extend_window_narrow_wide_adaptive_binomial_0.05/ \
#     ../results/Figrue2C_3B/Chipseeker_test 2>&1 | tail -50