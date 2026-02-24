#!/bin/bash
#SBATCH --job-name=hiple2
# #SBATCH --partition=gpu     
# #SBATCH --gres=gpu:1
#SBATCH --time=1-00:00:00
#SBATCH --mem=512G
# #SBATCH --cpus-per-task=8               
#SBATCH --output=/dcs07/hongkai/data/mjiang/hiplex/paper_materials/hiplex_paper/submit_all_steps_2_0.out
#SBATCH --error=/dcs07/hongkai/data/mjiang/hiplex/paper_materials/hiplex_paper/submit_all_steps_2_0.err
#SBATCH --mail-user=yhu157@jh.edu
#SBATCH --mail-type=END,FAIL

cd /dcs07/hongkai/data/mjiang/hiplex/paper_materials/hiplex_paper

bash /dcs07/hongkai/data/mjiang/hiplex/paper_materials/hiplex_paper/TF_binding_site_enrichment/run_TF_binding_site_enrichment.sh