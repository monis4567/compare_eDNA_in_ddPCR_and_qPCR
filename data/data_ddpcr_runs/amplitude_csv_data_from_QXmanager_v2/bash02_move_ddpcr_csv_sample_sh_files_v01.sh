#!/bin/bash
# -*- coding: utf-8 -*-

##get the present directory
WD=$(pwd)
INDIR1="/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2022/ddPCR_qPCR_MST/ddpcr_resultater"
# define directory to place the csv files from the QX manager v.1.2 software in
OUDIR1="csv_sample_sheets_from_QXmanager"
# remove the previous version of this output directory
rm -rf "$OUDIR1"
# make a new version of this output directory
mkdir "$OUDIR1"
# make a list that holds the prefix to the names of the ddpcr csv report files
LSSSHF=$(ls "$INDIR1" | grep 'ddpcr0' | grep '\.csv')
#echo $LSSSHF

#make the list of samples an array you can iterate over
declare -a FARRAY=("${LSSSHF}")

#iterate over prefix to the names of the ddpcr csv reports
# and make new directories and copy files in to these directories 
for F in ${FARRAY[@]}
do
	echo $F
	#
	cat "$INDIR1"/$F |  cut -d',' -f2 | sed -e 's/Sample description 1/Sample/g' > "$WD"/"$OUDIR1"/tmp_sample_column.txt
	paste -d',' "$INDIR1"/$F "$WD"/"$OUDIR1"/tmp_sample_column.txt | \
	sed 's/\r//' | \
	LC_ALL=C sed -e 's/\r$//' > "$WD"/"$OUDIR1"/tmpfl.txt
	
	cat "$WD"/"$OUDIR1"/tmpfl.txt | sed 's/ //g' > "$WD"/"$OUDIR1"/$F
	rm "$WD"/"$OUDIR1"/tmpfl.txt
	#cat "$F" | awk 'BEGIN {FS=","; OFS=","} {print $1,".",$2,",",$3,".", $4,",",$5}'
	#cp "$INDIR1"/$F "$WD"/"$OUDIR1"/.
done
# remove the temporary column file
rm -rf "$WD"/"$OUDIR1"/tmp_sample_column.txt

