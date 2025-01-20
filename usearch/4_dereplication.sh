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
./usearch11 -fastx_uniques ${sample}_filt.fa -sizeout -relabel ${sample}_Uniq -fastaout ${sample}_uniques.fa
done < sample_tab.txt
