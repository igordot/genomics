#!/bin/bash


##
## 10X Cell Ranger
## processes Chromium single cell RNA-seq output with expression and antibody libraries
##


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")

# check for correct number of arguments
if [ ! $# == 4 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name genome_name sample_name libraries_csv features_csv \n" >&2
	exit 1
fi

# arguments
genome_name=$1
sample_name=$2
sample_name_out="count-$sample_name"
libraries_csv=$(readlink -f "$3")
features_csv=$(readlink -f "$4")

# settings
# threads=$NSLOTS
# threads=$SLURM_NTASKS
threads=16
# mem=$(echo "$threads * 8" | bc)
mem=128

# make the output group-writeable
umask 007

if [[ "$genome_name" == "hg19" ]] ; then
	transcriptome_dir="/gpfs/data/sequence/cellranger-refdata/refdata-cellranger-hg19-3.0.0"
elif [[ "$genome_name" == "hg38" ]] ; then
	transcriptome_dir="/gpfs/data/sequence/cellranger-refdata/refdata-cellranger-GRCh38-3.0.0"
elif [[ "$genome_name" == "GRCh38" ]] ; then
	transcriptome_dir="/gpfs/data/sequence/cellranger-refdata/refdata-cellranger-GRCh38-3.0.0"
elif [[ "$genome_name" == "mm10" ]] ; then
	transcriptome_dir="/gpfs/data/sequence/cellranger-refdata/refdata-cellranger-mm10-3.0.0"
elif [[ "$genome_name" == "hg19_and_mm10" ]] ; then
	transcriptome_dir="/gpfs/data/sequence/cellranger-refdata/refdata-cellranger-hg19-and-mm10-3.0.0"
else
	transcriptome_dir="/gpfs/data/sequence/cellranger-refdata/ref/${genome_name}/cellranger"
	#transcriptome_dir="/gpfs/home/id460/ref/${genome_name}/cellranger"
fi

# unload all loaded modulefiles
module purge

# check that input exists
if [ ! -s "$libraries_csv" ] ; then
	echo -e "\n ERROR: libraries csv $libraries_csv does not exist \n" >&2
	exit 1
fi

if [ ! -s "$features_csv" ] ; then
	echo -e "\n ERROR: features csv $features_csv does not exist \n" >&2
	exit 1
fi

if [ ! -d "$transcriptome_dir" ] ; then
	echo -e "\n ERROR: genome dir $transcriptome_dir does not exist \n" >&2
	exit 1
fi

module load cellranger/3.1.0

# display settings
echo " * cellranger:        $(which cellranger) "
echo " * threads:           $threads "
echo " * mem:               $mem "
echo " * transcriptome dir: $transcriptome_dir "
echo " * libraries csv:     $libraries_csv "
echo " * features csv:      $features_csv "
echo " * sample:            $sample_name "
echo " * out dir:           $sample_name_out "

echo -e "\n $(date) \n" >&2

# cellranger run command

# id             A unique run id, used to name output folder
# fastqs         Path of folder created by 10x demultiplexing or bcl2fastq (must not be passed)
# sample         Prefix of the filenames of FASTQs to select
# transcriptome  Path of folder containing 10X-compatible transcriptome
# libraries      Path to a file declaring FASTQ paths and library types of input libraries
# feature-ref    Path to a file declaring the Feature Barcoding reagents

cellranger_cmd="
cellranger count \
--localmem $mem \
--localcores $threads \
--transcriptome $transcriptome_dir \
--libraries $libraries_csv \
--feature-ref $features_csv \
--id $sample_name_out \
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
rsync -tv "$web_summary_html" "./${sample_name_out}.html"

# delete temp files
rm -rf "./${sample_name_out}/SC_RNA_COUNTER_CS"

echo -e "\n $(date) \n"



# end
