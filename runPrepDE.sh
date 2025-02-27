#!/bin/bash

#SBATCH --job-name=prepDE
#SBATCH -p medium
#SBATCH -N 1
#SBATCH -n 1 #Tasks paralelos, recomendado para MPI, Default=1
#SBATCH --cpus-per-task=1 #Cores requeridos por task, recomendado para multi-thread, Default=1
#SBATCH --time=48:00:00
#SBATCH --mem=64000
#SBATCH --mail-user=j.stepanian@uniandes.edu.co
#SBATCH --mail-type=ALL
#SBATCH -o prepDE.o%j

python2 prepDE.py -i listaGTF.txt -g geneCountMatrix_hg38.csv -t transcriptsMatrix_hg38.csv
