#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 1

module load miniconda3

/usr/bin/time -v ./lncRNAsCVDistribution.py -lncRNAs Hoang_lncRNAs.tsv -CV /Storage/data1/felipe.peres/Sugarcane_ncRNA/9_Fiber_and_Sugar/co-expression/Hoang/code/updated_filters/CNC/Hoang2017_counts_filters_VST_CV.txt -out Hoang_lncRNAs_CV_Distribution.png

/usr/bin/time -v ./lncRNAsCVDistribution.py -lncRNAs Correr_lncRNAs.tsv -CV /Storage/data1/felipe.peres/Sugarcane_ncRNA/9_Fiber_and_Sugar/co-expression/Correr/code/updated_filters/CNC/Correr2020_counts_filters_VST_CV.txt  -out Correr_lncRNAs_CV_Distribution.png

/usr/bin/time -v ./lncRNAsCVDistribution.py -lncRNAs Perlo_lncRNAs.tsv -CV /Storage/data1/felipe.peres/Sugarcane_ncRNA/9_Fiber_and_Sugar/co-expression/Perlo/code/updated_filters/CNC/Perlo2022_counts_filters_VST_CV.txt -out Perlo_lncRNAs_CV_Distribution.png


