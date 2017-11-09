#!/bin/bash


##
## merge any number of tab or comma-separated files (coreutils join can only do 2 at a time)
## for tab field separator, use $'\t'
##


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")

# check for correct number of arguments
if [ $# -lt 3 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name field_separator missing_field_char in1.txt [in2.txt in3.txt ...] > merged.txt \n" >&2
	exit 1
fi

# load recent coreutils ("-o auto" support added in release 8.12)
module load coreutils/8.24

# arguments
separator="$1"
shift
empty_char="$1"
shift

# check if at least the first file exists
if [ ! -s "$1" ] ; then
	echo -e "\n $script_name ERROR: file $1 does not exist \n" >&2
	exit 1
fi

# recursive join function
function rjoin {
	if [[ $# -gt 1 ]]; then
		LC_ALL=C join -t "$separator" -a1 -a2 -o auto -e "$empty_char" - <(LC_ALL=C sort "$1") | rjoin "${@:2}"
	else
		LC_ALL=C join -t "$separator" -a1 -a2 -o auto -e "$empty_char" - <(LC_ALL=C sort "$1")
	fi
}

rjoin "${@:2}" < "$1"



# end
