#!/bin/bash

#PBS -N DeepPlnc_1parte
#PBS -q serial
#PBS -e DeepPlnc_1parte.sh.e
#PBS -o DeepPlnc_1parte.sh.o
#PBS -J 1-47:2
#PBS -r y

cd $PBS_O_WORKDIR

NONCODING1=`ls -1 *.cpc_ncrnas.fa | head -n $PBS_ARRAY_INDEX | tail -n1`

NONCODING2=`ls -1 *.cpc_ncrnas.fa | head -n$((1+${PBS_ARRAY_INDEX})) | tail -n1`

#change header
sed 's/\s.*$//g' $NONCODING1 > ${NONCODING1}_mod

sed 's/\s.*$//g' $NONCODING2 > ${NONCODING2}_mod

#change sequences
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' ${NONCODING1}_mod > ${NONCODING1}_mod_single

awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' ${NONCODING2}_mod > ${NONCODING2}_mod_single

source /home/lovelace/proj/proj832/fvperes/miniconda3/etc/profile.d/conda.sh
conda activate deepplnc_labis

echo ### >> jobs1.txt
echo "JOB " $PBS_ARRAY_INDEX >> jobs1.txt
echo "NONCODING1 =" ${NONCODING1}_mod_single >> jobs1.txt

echo "JOB " $((1+${PBS_ARRAY_INDEX})) >> jobs1.txt
echo "NONCODING2 =" ${NONCODING2}_mod_single >> jobs1.txt
echo ### >> jobs1.txt

python3 DeepPlnc_1parte.py ${NONCODING1}_mod_single

python3 DeepPlnc_1parte.py ${NONCODING2}_mod_single

rm ${NONCODING1}_mod
rm ${NONCODING2}_mod

rm ${NONCODING1}_mod_single
rm ${NONCODING2}_mod_single
