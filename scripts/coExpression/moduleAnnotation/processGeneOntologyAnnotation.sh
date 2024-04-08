#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -pe smp 1

module load miniconda3

annotation=GO_universe_annotation_list
panrnaome=panTranscriptome_panRNAomeClassificationTable_hyphen_Class.tsv
out=GO_annotations
ontology=MF
threshold=(0.3 0.4 0.5 0.6 0.7 0.8 0.9)

for i in "${threshold[@]}"
do
	echo threshold = $i
	./processGeneOntologyAnnotation.py --annotation $annotation -p $panrnaome -o $out -og $ontology -ppv $i
done
