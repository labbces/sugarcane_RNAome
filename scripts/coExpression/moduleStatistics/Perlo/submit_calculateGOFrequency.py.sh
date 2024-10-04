#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 1

module load miniconda3

/usr/bin/time -v ./calculateGOFrequency.py -modules ../moduleCNCwithlncRNAs.tsv -GO ./ -obo go-basic.obo -out_freq GOFrequencies.tsv -out_modules moduleCNCwithlncRNAswithEnrichedGO.tsv
