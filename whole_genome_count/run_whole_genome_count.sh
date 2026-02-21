


module load conda_R

align_dir="../data/valley-V-qc"

bin_size="800"
target_qc_type="all-qc"
frag_type="mixed"
libnorm_type_list=( "libnorm-mean" "libnorm-median" )
scen_raw="V"

scen=${scen_raw}_${frag_type}
work_dir="../results/Figure3"
mkdir -p ${work_dir}

for libnorm_type in "${libnorm_type_list[@]}"; do
    Rscript ../hiplex_paper/whole_genome_count/wgc_formal.R \
        ${scen} \
        ${frag_type} \
        ${bin_size} \
        ${target_qc_type} \
        ${libnorm_type} \
        ${work_dir} \
        ${align_dir} > ${work_dir}/result_${libnorm_type}.out 2>&1 &
done
wait
