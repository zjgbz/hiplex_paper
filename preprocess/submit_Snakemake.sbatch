#!/bin/bash
#SBATCH --job-name=Snakemake_Hiplex
#SBATCH --partition=cpu
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --mem=32g
#SBATCH --time=5-20:00:00
##SBATCH --mail-type=END,ALL
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.out

conda activate snakemake
cd /workdir/
SLURM_ARGS="-p cpu -n {resources.tasks} -c {resources.cpus_per_task} --mem={resources.mem_mb} -t {resources.runtime} --output=/dev/null --error=/dev/null"
which snakemake
snakemake --jobs 10 --configfile ./config.yaml --cluster "sbatch $SLURM_ARGS" --rerun-incomplete --latency-wait 150 --max-jobs-per-second 5 --max-status-checks-per-second 5 --retries 5
