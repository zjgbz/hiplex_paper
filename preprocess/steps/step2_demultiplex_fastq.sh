
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


barcode_fwd_fasta=`echo $1`
barcode_rev_fasta=`echo $2`
R1_fastq=`echo $3`
R2_fastq=`echo $4`
dir_demux_fastq=`echo $5`

base_name=$(basename ${dir_demux_fastq})
mkdir -p $TMPDIR/data/${base_name}
ulimit -n 65536

cutadapt \
		-e 2 \
		-j 8 \
		--no-indels \
		--action none \
    -g ^file:$barcode_fwd_fasta \
    -G ^file:$barcode_rev_fasta \
    -o $dir_demux_fastq/{name1}-{name2}.R1.fastq.gz \
    -p $dir_demux_fastq/{name1}-{name2}.R2.fastq.gz \
    $R1_fastq $R2_fastq

# cp -R $TMPDIR/data/${base_name}/. $dir_demux_fastq
