#!/bin/bash

#SBATCH --job-name=hisatAlign
#SBATCH -p medium
#SBATCH -N 1
#SBATCH -n 1 #Tasks paralelos, recomendado para MPI, Default=1
#SBATCH --cpus-per-task=16 #Cores requeridos por task, recomendado para multi-thread, Default=1
#SBATCH --time=72:00:00
#SBATCH --mem=100G
#SBATCH --mail-user=j.stepanian@uniandes.edu.co
#SBATCH --mail-type=ALL
#SBATCH -o hisatAlign.o%j


module load jdk/11.0.14 
module load hisat/2-2.2.1
module load samtools/1.16.1

#Se corre con for i in $(cat listaPE.txt); do sbatch runHisatPEoldmode.sh $i; done  -> la lista solo tiene el id de las muestras :)
#listaPE.txt contiene solo los ID de las muestras 

p=$1;


 # input files para samples JZ PE
#f1=DatosJovannyB1PE/${p}_R1_001.fastq.gz;
#f2=DatosJovannyB1PE/${p}_R2_001.fastq.gz;

 #input files para PE 
#f1=pairedEnd/${p}_1.fastq.gz
#f2=pairedEnd/${p}_2.fastq.gz

 # input files para SE
#f1=singleEnd/${p}.fastq.gz;


 # input files para colombianas SE
f1=colombianas/${p}.fastq.gz

 #Poner path a la ref
REFERENCE=reference/hg38_v0_Homo_sapiens_Idx; 

 # software variables. Write paths only if you can not install the programs or can not use installed versions
HISAT2=hisat2;
SAMTOOLS=samtools
JAVA=java;

 # jars for java packages | poner los paths
PICARD=/hpcfs/home/ing_sistemas/j.stepanian/software/picard_2.27.4.jar;
NGSEP=/hpcfs/home/ing_sistemas/j.stepanian/software/NGSEPcore_4.2.1.jar

 # map the reads and sort the alignment
mkdir ${p}_tmpdir;

 #hisat para PE
#${HISAT2} -p 16 --rg-id ${p} --rg SM:${p} --rg PL:ILLUMINA --dta -x ${REFERENCE} -1 ${f1} -2 ${f2} 2> ${p}_hisat2.log | ${JAVA} -Xmx8g -jar ${PICARD} SortSam MAX_RECORDS_IN_RAM=1000000 SO=coordinate CREATE_INDEX=true TMP_DIR=${p}_tmpdir I=/dev/stdin O=${p}_hisat2_sorted.bam >& ${p}_hisat2_sort.log;
rm -rf ${p}_tmpdir;


 #hisat para SE
${HISAT2} -p 16 --rg-id ${p} --rg SM:${p} --rg PL:ILLUMINA --dta -x ${REFERENCE} -U ${f1} 2>${p}_hisat2.log | ${JAVA} -Xmx8g -jar ${PICARD} SortSam MAX_RECORDS_IN_RAM=1000000 SO=coordinate CREATE_INDEX=true TMP_DIR=${p}_tmpdir I=/dev/stdin O=${p}_hisat2_sorted.bam >& ${p}_hisat2_sort.log;
rm -rf ${p}_tmpdir;

 # calculate statistics from the alignments file
${JAVA} -Xmx3g -jar ${NGSEP} BasePairQualStats -r ${REFERENCE} -o ${p}_readpos.stats ${p}_hisat2_sorted.bam >& ${p}_readpos.log
${JAVA} -Xmx3g -jar ${NGSEP} CoverageStats -i ${p}_sorted.bam -o ${p}_coverage.stats >& ${p}_coverage.log;
${SAMTOOLS} view -F 268 ${p}_hisat2_sorted.bam | awk '{l=$9;if(l>=0){i=sprintf("%d",l/25)+1;if(i<100)a[i]++;else aM++}}END{for(i=1;i<100;i++)print (i-1)*25,a[i];print "More",aM}' > ${p}_insertLength.stats;
