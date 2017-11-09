#!/bin/bash


##
## 10X Cell Ranger - processes Chromium single cell RNA-seq output
##


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")

# check for correct number of arguments
if [ ! $# == 3 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name genome_name sample_name fastq_dir \n" >&2
	exit 1
fi

# arguments
genome_name=$1
sample_name_fastq=$2
sample_name_out="count-$sample_name_fastq"
fastq_dir=$(readlink -f "$3")

# settings
threads=$NSLOTS
mem=$(echo "$threads * 8" | bc)
transcriptome_dir="/ifs/data/cellranger-refdata/refdata-cellranger-${genome_name}-1.2.0"
alt_transcriptome_dir="/ifs/home/id460/ref/${genome_name}/cellranger"

# check that input exists
if [ ! -d "$fastq_dir" ] ; then
	echo -e "\n ERROR: fastq dir $fastq_dir does not exist \n" >&2
	exit 1
fi

if [ ! -d "$transcriptome_dir" ] ; then
	echo -e "\n WARNING: genome dir $transcriptome_dir does not exist \n" >&2
	echo -e "\n setting genome dir to $alt_transcriptome_dir \n" >&2
	transcriptome_dir="${alt_transcriptome_dir}"
fi

if [ ! -d "$transcriptome_dir" ] ; then
	echo -e "\n ERROR: genome dir $transcriptome_dir does not exist \n" >&2
	exit 1
fi

# delete empty .po files to keep directory clean
rm -rf cellranger*.po*

module unload gcc
module load cellranger/2.0.0

# display settings
echo " * cellranger: $(which cellranger) "
echo " * threads: $threads "
echo " * mem: $mem "
echo " * transcriptome dir: $transcriptome_dir "
echo " * fastq dir: $fastq_dir "
echo " * sample: $sample_name_fastq "
echo " * out dir: $sample_name_out "

echo -e "\n $(date) \n" >&2

# cellranger run command

# id             A unique run id, used to name output folder
# fastqs         Path of folder created by 10x demultiplexing or bcl2fastq
# sample         Prefix of the filenames of FASTQs to select
# transcriptome  Path of folder containing 10X-compatible transcriptome

cellranger_cmd="
cellranger count \
--localcores $threads \
--localmem $mem \
--transcriptome $transcriptome_dir \
--fastqs $fastq_dir \
--sample $sample_name_fastq \
--id $sample_name_out
"
echo -e "\n CMD: $cellranger_cmd \n"
$cellranger_cmd

sleep 15

web_summary_html="./${sample_name_out}/outs/web_summary.html"

# check that output html summary (and probably everything else) exists
if [ ! -s "$web_summary_html" ] ; then
	echo -e "\n ERROR: summary $web_summary_html does not exist \n" >&2
	exit 1
fi

# copy html summary to top level for easy navigation
rsync -tv $web_summary_html ./${sample_name_out}.html

# delete empty .po files to keep directory clean
rm -rf cellranger*.po*

echo -e "\n $(date) \n"



# end
