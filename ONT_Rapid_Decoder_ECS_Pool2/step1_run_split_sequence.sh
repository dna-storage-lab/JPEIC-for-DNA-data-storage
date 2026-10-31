#!/bin/bash
echo "Begin Processing time : $(date '+%F %T')"
ulimit -s 655350

gcc -o split_fastq split_fastq.c

output_dir1=./ECS-Pool2_nanopore_data
if [ ! -d "$output_dir1" ]; then
    mkdir -p "$output_dir1"
fi

fasta_file=../ECS-Pool2_nanopore.fastq
need_reads=200000

 ./split_fastq $fasta_file $output_dir1 sub $need_reads 30

