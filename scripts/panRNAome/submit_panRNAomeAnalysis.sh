#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 1

module load miniconda3
conda activate panRNAomeAnalysis

/usr/bin/time -v ./panRNAomeAnalysis.py
