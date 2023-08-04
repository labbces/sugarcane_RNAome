#!/bin/bash

#PBS -N MMseqs2Cluster
#PBS -q par128
#PBS -l nodes=1:ppn=128
#PBS -e submitCluster_par128.sh.err
#PBS -o submitCluster_par128.sh.out

source /home/lovelace/proj/proj832/fvperes/miniconda3/etc/profile.d/conda.sh
conda activate mmseqs2

cd $PBS_O_WORKDIR

mmseqs createdb putative_ncRNA_consensus.fa DB

#DB - transcriptomes database
#DB_clust - clusterized database
#/work/fvperes - empty directory for temporary files

#running without --split-memory-limit 40G

mmseqs cluster --threads 128 -s 5.7 --cov-mode 2 --cluster-mode 2 -c 0.8 --min-seq-id 0.8 DB DB_clust /work/fvperes
mmseqs createsubdb DB_clust DB DB_clust_rep
mmseqs convert2fasta DB_clust_rep DB_clust_rep.fasta
mmseqs createtsv DB DB DB_clust DB_clust.tsv
