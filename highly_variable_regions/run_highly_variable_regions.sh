

module load conda_R

fitting_model="gam"
# para_list=("9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "111" "112" "113" "114" "115" "116" "117" "118" "119" "120" )
parameter="40"
robust_check="21"
seed="42"
nrow_sample="300000"
frag_type="mixed"
bin_size="800"
scen=V1V2
qc_type="all-qc" 
post_process="libnorm_noAllZero_log2_qnorm"

read_dir="../data/highly_variable_regions"
out_dir="../results/Figure3A"
mkdir -p ${out_dir}

Rscript ../hiplex_paper/highly_variable_regions/mav_screen_refine.R \
    ${fitting_model} \
    ${robust_check} \
    ${parameter} \
    ${seed} \
    ${nrow_sample} \
    ${frag_type} \
    ${bin_size} \
    ${scen} \
    ${qc_type} \
    ${post_process} \
    ${out_dir} \
    ${read_dir} > ${out_dir}/result_${fitting_model}_${robust_check}_${parameter}_${seed}_${nrow_sample}.out 2>&1
