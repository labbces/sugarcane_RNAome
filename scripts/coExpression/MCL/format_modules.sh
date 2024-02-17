#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 1

IN_FILE=`ls -1 *_mcl.txt`
OUT_FILE=pearson_mcl/out
num=($(LC_NUMERIC="en_US.UTF-8". seq 1.3 0.5 6))

# Get counts in each cluster
for i in "${num[@]}"
do
	echo $i
	awk -F' ' '{print NF}' ${OUT_FILE}.${i} > ${OUT_FILE}.${i}.number
done

# Reformat modules output

dir=pearson_mcl
cd $dir
for file in $(ls $dir | grep -v number | grep -v formated)
do
	echo $file	
	rm -fv ${dir}/${file}.formated.csv
	count=1
	for line in $(cat $file)
	do
        	sed -n "${count}p" $file | awk -F" " -v c="$count" '{for(i=1; i<=NF; i++) {print $i,c}}' >> ${file}.formated.csv
        	let count++
	done
done
cd -
