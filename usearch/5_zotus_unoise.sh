#!/bin/bash

#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBARCH --cpus-per-task=8
#SBATCH --mem-per-cpu=30G
#SBATCH --time=48:00:00
#SBATCH --job-name=usearch

#STEP1# ./usearch11 -unoise3 uniques_all.fa -zotus unoise_zotus.fa -tabbedout unoise_tab.txt
#STEP2# ./usearch11 -fastx_relabel unoise_zotus.fa -prefix zOTU -fastaout unoise_zotus_relabeled.fa -keep_annots
./usearch11 -otutab uniques_all.fa -zotus unoise_zotus.fa -otutabout unoise_otu_tab.txt -mapout unoise_map.txt -notmatched unoise_notmatched.fasta -dbmatched dbmatches.fasta -sizeout
