#!/bin/bash
#$ -q all.q
#$ -V
#$ -cwd
#$ -pe smp 1

module load Python/3.7.2

/usr/bin/time -v python3 pcc.py Hoang2017_tpm_vst_without80zeros.txt Hoang2017_tpm_vst_without80zeros_network.txt Hoang2017_tpm_vst_without80zeros_mcl_out.txt
