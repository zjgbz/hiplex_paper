#########  for developer  #######
## goal: do alignment
## code used: "/dcs05/hongkai/data/next_cutntag/workflows/modules/analysis_pipeline_cut_tag.sh"
## input: index fastq1 fastq2
## output: bam
## todo: 
#########   be careful    #######



#!/bin/bash

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG

cleanup() {
    exit_code=$?
    if [ ${exit_code} == 0 ]
    then
	echo "Completed execution"
    else
	echo "\"${last_command}\" failed with exit code ${exit_code}."
    fi
}

# echo an error message before exiting
trap 'cleanup' EXIT INT TERM

##usage: bash analysis_pipeline.sh read_1.fastq read_2.fastq outname index

# warning 1 jhce03 has no bowtie2, but jhpce03 has bowtie/2.5.1, can be load using `module load bowtie`
# warning 2 add samtools
module load bowtie2
module load bedtools
module load samtools

##trim adaptor
# java -mx4000M -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar PE $1 $2 $3.R1.paired $3.R1.unpaired $3.R2.paired $3.R2.unpaired ILLUMINACLIP:/users/wzhou14/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
# rm $3.R1.unpaired
# rm $3.R2.unpaired

##alignment
bowtie2 -x $4 -I 10 -X 800 -5 45 -p 8 -q --local --very-sensitive-local --no-unal --no-mixed --no-discordant -1 $1 -2 $2 -S $3.sam
samtools view -bS -o $3.bam $3.sam
rm $3.sam
# rm $3.R1.paired
# rm $3.R2.paired

##filter reads with map quality and filter unwanted chromosome
samtools view -b -q 10 $3.bam > $3.filtered.bam
samtools view -h $3.filtered.bam | sed '/chrM/d;/random/d;/chrUn/d' | samtools view -Sb - > $3.filtered.clean.bam
rm $3.filtered.bam
rm $3.bam

##remove PCR duplicates
# warning 3 jar need to be installed
# warning 4 add -o option in sort
samtools sort $3.filtered.clean.bam -o $3.sort.bam
java -Xmx16g -jar /projects/foundation_model_for_single_cell_multiomics_data/software/picard-tools-1.137/picard.jar MarkDuplicates INPUT=$3.sort.bam OUTPUT=$3.filtered.clean.dup.bam METRICS_FILE=$3_log_file.txt REMOVE_DUPLICATES=true
#java -Xmx16g -jar ~/picard-tools-1.137/picard.jar MarkDuplicates INPUT=$3.sort.bam OUTPUT=$3.filtered.clean.dup.bam METRICS_FILE=$3_log_file.txt REMOVE_DUPLICATES=true
samtools sort -n $3.filtered.clean.dup.bam -o $3.filtered.clean.dup.sort.bam
rm $3.sort.bam
rm $3.filtered.clean.bam
rm $3.filtered.clean.dup.bam
mv $3.filtered.clean.dup.sort.bam $3.bam
