#!bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -pe smp 1
#$ -t 1-33
#$ -tc 10

module load salmon/1.8.0

R1=`ls -1 ../Clean/*R1.fastq.gz | head -n $SGE_TASK_ID | tail -n1`
R2=${R1/R1.fastq.gz/R2.fastq.gz}

sample=${R1/R1.fastq.gz/}

salmon_index=/Storage/data1/riano/Sugarcane/Pantranscriptome/Assemblies/CD-HIT_48_genotypes_transcriptome_salmonInx/

/usr/bin/time -v salmon quant -p $NSLOTS -i $salmon_index -l A -1 $R1 -2 $R2 -o $sample
