#!/bin/bash

#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBARCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=40G
#SBATCH --time=4:00:00
#SBATCH --job-name=phylogeny

module load mafft/7.481-g3ca 
module load trimal/1.4.1-wan5
module load iqtree/1.6.12

fasta=$1
aln=$(echo ${fasta} | sed 's/.fasta/.msa/')
trim=${aln}.trim

mafft --auto $fasta > $aln
trimal -automated1 -in $aln -out $trim
iqtree -s $1 -m MF
#iqtree -s $1 -m TIM2e+G4 -bb 10000
