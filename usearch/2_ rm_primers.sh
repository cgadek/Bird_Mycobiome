#!/bin/bash

#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBARCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --time=48:00:00
#SBATCH --job-name=usearch

#ITS3F primer GCATCGATGAAGAACGCAGC
#ITS4R primer TCCTCCGCTTATTGATATGC

IFS=$'\t'

while read sample R1 R2; do
./usearch11 -fastx_truncate ${sample}_merged_relabeled.fq -stripleft 20 -stripright 20 -fastqout ${sample}_stripped.fq
done < sample_tab.txt
