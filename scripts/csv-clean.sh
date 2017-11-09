#!/bin/bash


##
## make csv file more compatible with various tools by removing problematic characters
##


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")

# check for correct number of arguments
if [ ! $# == 1 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name file.csv \n" >&2
	exit 1
fi

# arguments
csv=$1

# check that input exists
if [ ! -s "$csv" ] ; then
	echo -e "\n $script_name ERROR: file $csv does not exist \n" >&2
	exit 1
fi

# fix newlines
dos2unix --quiet $csv
mac2unix --quiet $csv

# replace commas inside quoted fields with dashes
awk -F '"' -v OFS='' '{ for (i=2; i<=NF; i+=2) gsub(",", "-", $i) } 1' $csv > ${csv}.tmp && mv ${csv}.tmp $csv

# remove quotes
sed -i 's/\"//g' $csv

# remove lines missing any values (only commas present)
sed -i '/^,,*$/d' $csv

# add newline to end of file if one does not exist (some scripts may complain)
sed -i -e '$a\' $csv



# end
