#!/usr/bin/bash 

# require cisgenome from https://www.biostat.jhsph.edu/~hji/cisgenome/index.htm
input_file=$1
out_dir=$2


file_name=$(basename "${input_file%.*}")

temp_dir=${out_dir}/${file_name}/temp/
control_dir=${out_dir}/${file_name}/control/
homer_dir=${out_dir}/${file_name}/homer_result/
mkdir -p ${temp_dir}
mkdir -p ${control_dir}
mkdir -p ${homer_dir}


input_cod_file="${temp_dir}/${file_name}.cod"
ouput_cod_file="${control_dir}/${file_name}_output_matchcontrol.cod"
ouput_bed_file="${control_dir}/${file_name}_output_matchcontrol.bed"
# Convert .bed to .cod
/users/kyu1/softwares/cisGenome-2.0/bin/file_bed2cod -i ${input_file} -o ${input_cod_file} 

# get matched genomic control regions
/users/kyu1/softwares/cisGenome-2.0/bin/refgene_getmatchedcontrol -d /dcs05/hongkai/data/wzhou14/tools/cisgenome/hg38/refFlat_sorted.txt -dt 1 -s human -c /dcs05/hongkai/data/wzhou14/tools/cisgenome/hg38/chrlen.txt -i $input_cod_file -o $ouput_cod_file -n 1 -l 200 -nr 1
# /users/kyu1/softwares/cisGenome-2.0/bin/refgene_getmatchedcontrol -d /dcs05/hongkai/data/wzhou14/tools/cisgenome/hg38/refFlat_sorted.txt -dt 1 -s human -c /dcs05/hongkai/data/wzhou14/tools/cisgenome/hg38/chrlen.txt -i /dcs05/hongkai/data/next_cutntag/script/motif_analysis/test_run_result/background/temp/A.cod -o "/dcs05/hongkai/data/next_cutntag/script/motif_analysis/test_run_result/background/control/A_output_matchcontrol.cod" -n 1 -l 1000 -nr 1

# Convert background.cod to background.bed
/users/kyu1/softwares/cisGenome-2.0/bin/file_cod2bed -i $ouput_cod_file -o $ouput_bed_file
