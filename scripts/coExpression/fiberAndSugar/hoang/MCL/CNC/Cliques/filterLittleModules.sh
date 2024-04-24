#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -pe smp 1

dataset=Hoang2017
modules=$(ls -1 clique* | wc -l)
min_module_size=5

# generate formated modules
awk '{print $1, $2, substr(FILENAME, 8, length(FILENAME)-11)}' clique* | sort -k3n >> ${dataset}_counts_filters_VST_topCV_mcl_formated_cliques.csv

# create sequence with modules from 1 to X modules - filtering modules with less than 5 genes
seq 1 $modules | grep -vw "$(cut -f3 -d " " ${dataset}_counts_filters_VST_topCV_mcl_formated_cliques.csv | uniq -c | awk "\$1 < $min_module_size {print \$2}")" >> Hoang2017_counts_filters_VST_topCV_mcl_formated_cliques_filtered.csv


