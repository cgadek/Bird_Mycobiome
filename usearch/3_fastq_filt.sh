#!/bin/bash

#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBARCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --time=48:00:00
#SBATCH --job-name=usearch

IFS=$'\t'

while read sample R1 R2; do
./usearch11 -fastq_filter ${sample}_stripped.fq -fastaout ${sample}_filt.fa -fastq_maxee 1 -fastq_minlen 150 
done < sample_tab.txt
