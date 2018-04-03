#!/bin/bash


##
## 10X Cell Ranger
## cellranger aggr - aggregates count data from multiple runs of the 'cellranger count'
##


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")

# check for correct number of arguments
if [ $# -lt 1 ] || [ $# -gt 2 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name sample_sheet [name] \n" >&2
	exit 1
fi

# arguments
sample_sheet=$1
analysis_name=$2

# settings (many sub-steps seem to be single-threaded, so threads are mostly irrelevant)
threads="4"
mem="32"

# output (add analysis name if provided)
sample_name="aggregated"
if [ -n "$analysis_name" ] ; then
	sample_name="${sample_name}-${analysis_name}"
fi
web_summary_html="${sample_name}/outs/web_summary.html"

# check that input exists
if [ ! -s "$sample_sheet" ] ; then
	echo -e "\n ERROR: sample sheet $sample_sheet does not exist \n" >&2
	exit 1
fi

# delete empty .po files to keep directory clean
rm -rf cellranger*.po*

echo -e "\n $(date) \n" >&2

# check if output already exits
if [ -s "$web_summary_html" ]; then
	echo -e "\n ERROR: summary $web_summary_html already exists \n" >&2
	exit 1
fi

# clean up sample sheet
dos2unix -q "$sample_sheet"
sed -i 's/"//g' "$sample_sheet"
sed -i -e '$a\' "$sample_sheet"

module unload gcc
module load cellranger/2.1.0

# display settings
echo " * cellranger: $(which cellranger) "
echo " * sample sheet: $sample_sheet "

# cellranger aggr command

# id      A unique run id, used to name output folder [a-zA-Z0-9_-]+.
# csv     Path of CSV file enumerating 'cellranger count' outputs.

cellranger_cmd="
cellranger aggr \
--jobmode local \
--localcores $threads \
--localmem $mem \
--id $sample_name \
--csv $sample_sheet
"
echo -e "\n CMD: $cellranger_cmd \n"
$cellranger_cmd

sleep 15

# check that output html summary (and probably everything else) exists
if [ ! -s "$web_summary_html" ] ; then
	echo -e "\n ERROR: summary $web_summary_html does not exist \n" >&2
	exit 1
fi

# copy html summary to top level for easy navigation
rsync -t "$web_summary_html" "./${sample_name}.html"

# clean up (temp files)
rm -rf "${sample_name}/SC_RNA_COUNTER_CS"

echo -e "\n $(date) \n"



# end
