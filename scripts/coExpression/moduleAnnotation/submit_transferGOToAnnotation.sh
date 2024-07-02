#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -pe smp 1

module load miniconda3

/usr/bin/time -v ./transferGOToAnnotation.py
