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
./usearch11 -fastq_mergepairs ${R1} \
	-reverse ${R2} \
	-fastqout ${sample}_merged.fq \
	-fastq_minmergelen 200 \
	-fastq_maxdiffs 20 \
	-fastq_pctid 10 \
	-report ${sample}_merge_report.txt \
	-tabbedout ${sample}_merged_tabedbout.txt
done < sample_tab.txt
