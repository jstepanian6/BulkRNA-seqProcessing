#!/bin/bash

#SBATCH --job-name=Stringtie
#SBATCH -p medium
#SBATCH -N 1
#SBATCH -n 1 #Tasks paralelos, recomendado para MPI, Default=1
#SBATCH --cpus-per-task=1 #Cores requeridos por task, recomendado para multi-thread, Default=1
#SBATCH --time=4:00:00
#SBATCH --mem=8000
#SBATCH --mail-user=j.stepanian@uniandes.edu.co
#SBATCH --mail-type=ALL
#SBATCH -o Stringtie.o%j

module load anaconda/python3.9
source activate stringTie2.2.1 

#se corre con el mismo for de hisat :)
p=$1;
o=/hpcfs/home/ing_sistemas/j.stepanian/PostMaster/Fertilidad/Processing_hg38/stringTieOutput
#input files
f1=/hpcfs/home/ing_sistemas/j.stepanian/PostMaster/Fertilidad/Processing_hg38/${p}_hisat2_sorted.bam
#referencia
REFERENCE=/hpcfs/home/ing_sistemas/j.stepanian/remapeo/reference/UCSC_hg38.gff

stringtie -G $REFERENCE -e -b ${p} -l ${p} -A ${o}/${p}_stringtie_abundances.txt -o ${o}/${p}_stringtie.gtf ${f1} 

