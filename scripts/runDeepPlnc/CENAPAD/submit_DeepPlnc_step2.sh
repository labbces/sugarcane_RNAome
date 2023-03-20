#!/bin/bash

#PBS -N DeepPlnc_2parte
#PBS -q serial
#PBS -e DeepPlnc_2parte.sh.e
#PBS -o DeepPlnc_2parte.sh.o
#PBS -J 1-47:2
#PBS -r y

cd $PBS_O_WORKDIR

source /home/lovelace/proj/proj832/fvperes/miniconda3/etc/profile.d/conda.sh
conda activate deepplnc_labis

#Path RNAFold bin
RNAfold=/home/lovelace/proj/proj832/fvperes/Sugarcane_ncRNA/ViennaRNA-2.5.1/bin/RNAfold

#Input transcriptome
#FILE_INPUT="test"

FILE_INPUT1=`ls -1 ../step1/*fa_mod_single_coding.fa | head -n $PBS_ARRAY_INDEX | tail -n1`

FILE_INPUT2=`ls -1 ../step1/*fa_mod_single_coding.fa | head -n$((1+${PBS_ARRAY_INDEX})) | tail -n1`

BASENAME_INPUT1=`basename ${FILE_INPUT1}`

BASENAME_INPUT2=`basename ${FILE_INPUT2}`

#Output from DeepPlnc
TEMP_FILE1=${BASENAME_INPUT1}_temp

TEMP_FILE2=${BASENAME_INPUT2}_temp

#Python
#os.system(path+"RNAfold -j --noPS  "+filename+"_coding.fa | paste - - - | awk \'{print $3}\' > "+filename+"_RNA_coding")

#Bash
$RNAfold --jobs=1 --noPS $FILE_INPUT1 | paste - - - | awk '{print $3}' > ${BASENAME_INPUT1}_RNA_coding

$RNAfold --jobs=1 --noPS $FILE_INPUT2 | paste - - - | awk '{print $3}' > ${BASENAME_INPUT2}_RNA_coding

#Python
#filename_2 = str("temp")
#os.system("cat "+filename+"_coding.fa | paste - - "+filename+"_RNA_coding >"+filename_2)

#Bash
cat $FILE_INPUT1 | paste - - ${BASENAME_INPUT1}_RNA_coding > $TEMP_FILE1

cat $FILE_INPUT2 | paste - - ${BASENAME_INPUT2}_RNA_coding > $TEMP_FILE2

echo "BASENAME_INPUT1 =" $BASENAME_INPUT1 >> jobs2.txt
echo "BASENAME_INPUT2 =" $BASENAME_INPUT2 >> jobs2.txt
echo "BASENAME_INPUT1_process =" ${FILE_INPUT1} >> jobs2.txt
echo "BASENAME_INPUT2_process =" ${FILE_INPUT2} >> jobs2.txt
echo "BASENAME_INPUT1_RNAfold =" ${BASENAME_INPUT1}_RNA_coding >> jobs2.txt
echo "BASENAME_INPUT2_RNAfold =" ${BASENAME_INPUT2}_RNA_coding >> jobs2.txt
echo "BASENAME_INPUT1_final =" $TEMP_FILE1 >> jobs2.txt
echo "BASENAME_INPUT2_final =" $TEMP_FILE2 >> jobs2.txt
