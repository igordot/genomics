#!/bin/bash


##
## De novo assembly from 10x Genomics Chromium Linked-Reads using Supernova.
##
## Usage:
## sbatch --job-name=supernova --nodes=1 --ntasks=1 --cpus-per-task=17 --mem-per-cpu=32G --time=15-00 \
## --partition=fn_long --mail-user=${USER}@nyulangone.org --mail-type=END,FAIL,REQUEUE --export=NONE \
## --wrap="bash /gpfs/data/igorlab/public/genomics/scripts-bigpurple/assembly-10x-supernova.sh fastq_dir [max_reads]"
##


#########################


# system-specific settings

# supernova directory
supernova_version="2.1.1"
supernova_dir="/gpfs/data/igorlab/software/supernova/supernova-${supernova_version}"

# prepare environment
module purge
module add default-environment


#########################


# check for correct number of arguments
if [ $# -lt 1 ] ; then
	echo -e "\n ERROR: wrong number of arguments supplied \n" >&2
	echo -e "\n USAGE: bash assembly-10x-supernova.sh fastq_dir [max_reads] \n" >&2
	exit 1
fi

# arguments
fastq_dir=$(readlink -f "$1")
max_reads="$2"

# check that input exists
if [ ! -d "$fastq_dir" ] ; then
	echo -e "\n ERROR: fastq dir $fastq_dir does not exist \n" >&2
	exit 1
fi


#########################


# step 1: assembly (supernova run)

# million reads cutoff (default is 1200M)
# set the number of reads so as to achieve 56x raw coverage: (genome size) x 56 / 150, assuming 150bp reads
# coverage significantly greater than 56x can sometimes help but can also be deleterious, depending on the dataset
# default value is 1.2B, which only makes sense for ~3.2 Gb genomes
if [ -n "$max_reads" ] ; then
	max_reads_m="$max_reads"
else
	max_reads_m="1200"
fi

# system load settings (reserve an extra thread for overhead)
threads=$SLURM_CPUS_PER_TASK
threads=$(( threads - 1 ))
mem=$(echo "$threads * 32" | bc)

# run name (used to name output folder)
supernova_version_nodot=$(echo "$supernova_version" | sed 's/\.//g')
run_id="assembly-supernova-v${supernova_version_nodot}-reads${max_reads_m}M"

# display settingse
echo
echo " * fastq dir:               $fastq_dir "
echo " * supernova bin dir:       $supernova_dir "
echo " * reads cutoff (million):  $max_reads_m "
echo " * threads:                 $threads "
echo " * mem:                     $mem "
echo " * run name (output dir):   $run_id "
echo

echo -e "\n assembly started: $(date) \n" >&2

# supernova assembly command

supernova_cmd="
${supernova_dir}/supernova run \
--maxreads ${max_reads_m}000000 \
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


# step 2: generate fasta file (supernova mkoutput)

# display settings
echo
echo " * assembly dir:  ${supernova_out_dir}/outs/assembly "
echo " * fasta prefix:  ${supernova_out_dir}/assembly "
echo

# generate different style fasta files
styles="raw megabubbles pseudohap pseudohap2"
for s in $styles; do

	echo -e "\n generate fasta: style $s \n" >&2

	# supernova mkoutput command
	${supernova_dir}/supernova mkoutput \
	--asmdir "${supernova_out_dir}/outs/assembly" \
	--outprefix "${supernova_out_dir}/assembly.${s}" \
	--style "${s}"

done

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


# cleanup

# check file size before cleanup
echo
echo "file size before cleanup"
du -sh "$run_id"
echo

# delete large assembly files (keep small ones just in case)
rm -rf ${run_id}/outs/assembly/a*
rm -rf ${run_id}/outs/assembly/closures*
rm -rf ${run_id}/outs/assembly/data
# delete temp files
rm -rf ${run_id}/ASSEMBLER_CS

# check file size after cleanup
echo
echo "file size after cleanup"
du -sh "$run_id"
echo


#########################



# end
