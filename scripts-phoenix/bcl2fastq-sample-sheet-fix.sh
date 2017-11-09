#!/bin/bash


##
## fix bcl2fastq demultiplexing sample sheet to get rid of problematic characters
##


# input
proj=$1

if [ -z "$1" ]
then
	echo "ERROR! NO ARGUMENT SUPPLIED."
	exit 1
fi

basecalls_dir="/ifs/data/sequence/Illumina/production/${proj}/Data/Intensities/BaseCalls"
ss=${basecalls_dir}/SampleSheet.csv

printf "\n\n FIX SAMPLE SHEET $ss \n\n"

# make the output group-writeable
umask 007

# fix and show sample sheet
if [ -s $ss ]
then

	# fix newlines
	dos2unix --quiet $ss
	mac2unix --quiet $ss

	# replace commas inside quoted fields with dashes
	awk -F '"' -v OFS='' '{ for (i=2; i<=NF; i+=2) gsub(",", "-", $i) } 1' $ss > ${ss}.tmp && mv ${ss}.tmp $ss
	# replace periods, parentheses, quotes, and blanks in sample names with dashes
	awk -F ',' 'BEGIN { OFS="," } { gsub(/\.|\(|\)|\#|\:|\/|\047|[[:blank:]]/, "-", $3); print }' $ss > ${ss}.tmp && mv ${ss}.tmp $ss
	# replace multiple dashes
	sed -i 's/--*/-/g' $ss
	# replaces dashes at the beginning of the field
	sed -i 's/,-/,/g' $ss
	# replaces dashes at the end of the field
	sed -i 's/-,/,/g' $ss
	# remove lines missing values
	sed -i '/^,,,,,/d' $ss
	# add newline to end of file if one does not exist (some scripts may complain)
	sed -i -e '$a\' $ss

	# check for extra columns in sample sheet
	max_comma_count=0
	while read i
	do
		comma_count=$(echo $i | tr -d -c "," | wc -c)
		if [ $comma_count -gt $max_comma_count ]
		then
			max_comma_count=$comma_count
		fi
	done < $ss

	# if too many commas, replace trailing commas
	if [ $max_comma_count -gt 9 ]
	then
		sed -i 's/,,*$//g' $ss
	fi

	# display sample sheet for easy review
	column -s "," -t $ss

else

	printf "\n\n NO SAMPLE SHEET FOUND AT $ss \n\n"
	sleep 5
	exit 1

fi



# end
