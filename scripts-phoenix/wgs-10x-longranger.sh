#!/bin/bash


##
## Whole Genome Phasing and SV Calling from 10x Genomics Chromium Linked-Reads using Long Ranger.
##
## Usage:
## qsub -N longranger -M ${USER}@nyumc.org -m ae -j y -cwd -pe threaded 16 -b y \
## -hard -l mem_free=128G -l mem_token=8G \
## bash /ifs/home/id460/public/genomics/scripts-phoenix/wgs-10x-longranger.sh <sample> <fastq_dir>
##


#########################


# check for correct number of arguments
if [ ! $# == 2 ] ; then
	echo -e "\n ERROR: wrong number of arguments supplied \n" >&2
	echo -e "\n USAGE: bash wgs-10x-longranger.sh sample fastq_dir \n" >&2
	exit 1
fi

# arguments
sample="$1"
fastq_dir=$(readlink -f "$2")

# check that input exists
if [ ! -d "$fastq_dir" ] ; then
	echo -e "\n ERROR: fastq dir $fastq_dir does not exist \n" >&2
	exit 1
fi


#########################


# system-specific settings

# Long Ranger directory
longranger_version="2.2.2"
longranger_dir="/ifs/home/id460/software/longranger/longranger-${longranger_version}"

# Long Ranger reference directory
longranger_ref_dir="/ifs/home/id460/ref/hg38/longranger-2.1.0"

# GATK path (Long Ranger 2.2 compatible with versions 3.3-4.0, excluding 3.6)
gatk_jar="/ifs/home/id460/software/GenomeAnalysisTK/gatk-4.0.4.0/gatk-package-4.0.4.0-local.jar"

# unload all loaded modulefiles
module purge
module load local

# load java (for GATK)
module load java/1.8


#########################


# system load settings
threads=$NSLOTS
mem=$(echo "$threads * 8" | bc)

# output settings
run_id="longranger-${sample}"

# display settings
echo
echo " * sample:               $sample "
echo " * FASTQ dir:            $fastq_dir "
echo " * Long Ranger bin dir:  $longranger_dir "
echo " * GATK jar file:        $gatk_jar "
echo " * threads:              $threads "
echo " * mem:                  $mem "
echo " * run output:           $run_id "
echo

echo -e "\n analysis started: $(date) \n" >&2

longranger_cmd="
${longranger_dir}/longranger wgs \
--fastqs ${fastq_dir} \
--sample ${sample} \
--id ${run_id} \
--reference ${longranger_ref_dir} \
--vcmode=gatk:${gatk_jar} \
--localcores ${threads} \
--localmem ${mem} \
"
echo -e "\n CMD: $longranger_cmd \n"
$longranger_cmd

longranger_out_dir=$(readlink -f "$(pwd)/${run_id}")

echo -e "\n analysis ended: $(date) \n" >&2


#########################


# check that output generated

if [ ! -e "${longranger_out_dir}/outs/summary.csv" ] ; then
	echo -e "\n ERROR: output ${longranger_out_dir}/outs/summary.csv does not exist \n" >&2
	exit 1
fi

if [ ! -e "${longranger_out_dir}/outs/loupe.loupe" ] ; then
	echo -e "\n ERROR: output ${longranger_out_dir}/outs/loupe.loupe does not exist \n" >&2
	exit 1
fi


#########################


# cleanup

# check file size before cleanup
du -sh "$run_id"

# delete temp files
rm -rf "${run_id}/PHASER_SVCALLER_CS"

# check file size after cleanup
du -sh "$run_id"


#########################



# end
