#!/bin/bash


##
## De novo assembly from 10x Genomics Chromium Linked-Reads using Supernova.
##
## Usage:
## qsub -N supernova -M ${USER}@nyumc.org -m e -j y -cwd -pe threaded 16 -b y \
## -hard -l mem_free=256G -l mem_token=16G \
## bash /ifs/home/id460/public/genomics/scripts-phoenix/assembly-10x-supernova.sh <fastq_dir>
##


#########################


# system-specific settings

# supernova directory
supernova_dir="/ifs/home/id460/software/supernova/supernova-1.2.1"


#########################


# check for correct number of arguments
if [ ! $# == 1 ] ; then
	echo -e "\n ERROR: wrong number of arguments supplied \n" >&2
	echo -e "\n USAGE: bash assembly-10x-supernova.sh fastq_dir \n" >&2
	exit 1
fi

# arguments
fastq_dir=$(readlink -f "$1")

# check that input exists
if [ ! -d "$fastq_dir" ] ; then
	echo -e "\n ERROR: fastq dir $fastq_dir does not exist \n" >&2
	exit 1
fi


#########################


# step 1: assembly

# run name (used to name output folder)
run_id="assembly-supernova"

# supernova settings
threads=$NSLOTS
mem=$(echo "$threads * 16" | bc)

# display settings
echo " * fastq dir: $fastq_dir "
echo " * supernova bin dir: $supernova_dir "
echo " * threads: $threads "
echo " * mem: $mem "

echo -e "\n assembly started: $(date) \n" >&2

# supernova assembly command
# adjust --bcfrac and --maxreads for smaller genomes
# "Sequence 800-1200 M reads for a 3.2 Gb genome ... (400-600 M) for a 1.6 Gb genome"
supernova_cmd="
${supernova_dir}/supernova run \
--localcores ${threads} \
--localmem ${mem} \
--id ${run_id} \
--fastqs ${fastq_dir}
"
echo -e "\n CMD: $supernova_cmd \n"
$supernova_cmd

echo -e "\n assembly ended: $(date) \n" >&2

# check that output generated
supernova_out_dir=$(readlink -f "$(pwd)/${run_id}")
if [ ! -e "${supernova_out_dir}/outs/report.txt" ] ; then
	echo -e "\n ERROR: output ${supernova_out_dir}/outs/report.txt does not exist \n" >&2
	exit 1
fi


#########################


# step 2: generate fasta file

# display settings
echo " * assembly dir: ${supernova_out_dir}/outs/assembly "
echo " * fasta prefix: ${supernova_out_dir}/assembly "

# generate different style fasta files
styles="raw megabubbles pseudohap pseudohap2"
for s in $styles; do

	# supernova mkoutput command
	${supernova_dir}/supernova mkoutput \
	--asmdir "${supernova_out_dir}/outs/assembly" \
	--outprefix "${supernova_out_dir}/assembly.${s}" \
	--style "${s}"

done

# check that output generated
if [ ! -e "${supernova_out_dir}/assembly.raw.fasta.gz" ] ; then
	echo -e "\n ERROR: output ${supernova_out_dir}/assembly.raw.fasta.gz does not exist \n" >&2
	exit 1
fi

# check that output generated
styles_out="raw megabubbles pseudohap pseudohap2.1 pseudohap2.2"
for s in $styles_out; do

	# check that output generated
	if [ ! -e "${supernova_out_dir}/assembly.${s}.fasta.gz" ] ; then
		echo -e "\n ERROR: output ${supernova_out_dir}/assembly.${s}.fasta.gz does not exist \n" >&2
		exit 1
	fi

done


#########################


# clean up (temp files)
rm -rf ${run_id}/ASSEMBLER_CS
# clean up (large assembly files)
rm -rf ${run_id}/outs/assembly/a*
rm -rf ${run_id}/outs/assembly/closures*
rm -rf ${run_id}/outs/assembly/data


#########################



# end
