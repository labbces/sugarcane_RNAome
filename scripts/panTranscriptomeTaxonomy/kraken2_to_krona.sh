#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -pe smp 1

ktImportTaxonomy -t 5 -m 3 -o 48_transcriptomes.report.krona.html 48_transcriptomes.report.kraken
