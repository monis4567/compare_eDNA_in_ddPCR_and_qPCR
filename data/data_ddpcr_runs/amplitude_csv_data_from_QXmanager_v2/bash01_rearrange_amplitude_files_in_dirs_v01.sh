#!/bin/bash
# -*- coding: utf-8 -*-

##get the present directory
WD=$(pwd)

INDIR1="amplitude_data_from_QX_manager_v2"
TMPDIR="tmp_dir"

#remove the old versions of the in- and output directory
rm -rf ${TMPDIR}
#make new versions of the in- and output directory
#mkdir ${TMPDIR}

cd "$WD"/"$INDIR1"
# make a list that holds the prefix to the names of the ddpcr csv report files
LSSMPL=$(ls | cut -d'_' -f1-2 | uniq | grep 'dd')

#make the list of samples an array you can iterate over
declare -a SMPLARRAY=("${LSSMPL}")

#iterate over prefix to the names of the ddpcr csv reports
# and make new directories and copy files in to these directories 
for D in ${SMPLARRAY[@]}
do
	ND=$(echo "$D" | sed 's/Myaara/Myaare/g' | sed 's/ddPCR0/ddpcr0/g')

	echo "$D"
    #done
	rm -rf "$WD"/"$ND"
	mkdir "$WD"/"$ND"
	#
	F=$(ls | grep "${D}")
	# notice that $F needs to have no quotation marks around it
	# see this question: https://unix.stackexchange.com/questions/405527/i-get-message-file-name-too-long-when-running-for-in-and-touch
	for file in $F
		do
			Nfile=$(echo "$file" | sed 's/Myaara/Myaare/g' | sed 's/ddPCR0/ddpcr0/g')
		
		#echo $file
		# # this question : https://stackoverflow.com/questions/42004482/shell-script-replace-a-specified-column-with-sed
		# replace in lines in file to get point as delimeter between 1st and 2nd comma separated element
		# and btween 3rd and 4th comma separated element
		cat "$file" | awk 'BEGIN {FS=","; OFS=""} {print $1,".",$2,",",$3,".", $4,",",$5}' | \
		sed -e 's:1Amplitude\.Ch2Amplitude:1Amplitude,Ch2Amplitude:g' | \
		sed 's/\r//' | \
		LC_ALL=C sed -e 's/\r$//' | \
		sed 's/\.,//g' > "$WD"/"$ND"/"$Nfile"
		# To SKIP the first N lines:
		# https://stackoverflow.com/questions/604864/print-a-file-skipping-the-first-x-lines-in-bash
		# tail -n +5 > "$WD"/"$D"/"$file"
		# #sed ':a;N;$!ba;s/\n//g' 
		#cp $file "$WD"/"$D"/.
		#cat "$file" | head -3

		#check the last line of the output file
		lastline=$(tail -1 "$WD"/"$ND"/"$Nfile")
		if [[ "$lastline" = *Amplitude* ]]; then
			#echo "$file"
			#echo "$lastline"
			echo "0,0,0" >> "$Nfile"
		fi
		
		done 
done


#head -12 "$WD"/"ddpcr0023_Psefar/ddpcr0023_Psefar_20221025_123517_481_A06_Amplitude.csv"

