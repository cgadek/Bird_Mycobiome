#!/bin/bash

#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBARCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --time=48:00:00
#SBATCH --job-name=usearch

./usearch11 -sintax zotus_v2.fa -db ../../unite_alleuk_04.04.20204_V10.fasta.gz -tabbedout zotus.sintax.unite10 -strand both -sintax_cutoff 0.8
