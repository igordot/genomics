#!/bin/bash


##
## Processes 10x Genomics Chromium single-cell RNA-seq FASTQs with Cell Ranger (cellranger count).
## Provide a FASTQ directory for classic expression libraries.
## Provide tables of libraries and features for expression and antibody/hashtag libraries.
##
## Usage:
## sbatch --job-name=cellranger-${sample} --ntasks=1 --cpus-per-task=17 --mem=128G --time=10:00:00 \
##   --mail-user=${USER}@nyumc.org --mail-type=FAIL,END --export=NONE \
##   --wrap="bash ./scrna-10x-cellranger-count.sh module_version genome_build sample_name fastq_dir"
##


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")

# check for correct number of arguments
if [ $# -lt 4 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e " Usage:" >&2
	echo -e "" >&2
	echo -e "   RNA only:        ./${script_name} module_version genome_name sample_name fastq_dir" >&2
	echo -e "   RNA and ADT/HTO: ./${script_name} module_version genome_name sample_name libraries_csv features_csv" >&2
	echo -e "" >&2
	echo -e "   RNA example:         ./${script_name} 9.0.0 hg38 my_sample /gpfs/data/fastq" >&2
	echo -e "" >&2
	if [ $# -gt 0 ] ; then echo -e " Provided arguments: $* \n" >&2 ; fi
	exit 1
fi

# arguments
module_version=$1
genome_name=$2
sample_name_fastq=$3
sample_name_out="count-$sample_name_fastq"
if [ $# -eq 4 ] ; then
	fastq_dir=$(readlink -f "$4")
fi
if [ $# -eq 5 ] ; then
	libraries_csv=$(readlink -f "$4")
	features_csv=$(readlink -f "$5")
fi

# settings (16 threads and 64G does not finish with 10h time limit)
threads=16
mem=128

# make the output group-writeable
umask 007

if [[ "$genome_name" == "hg19" ]] ; then
	transcriptome_dir="/gpfs/data/sequence/cellranger-refdata/refdata-cellranger-hg19-3.0.0"
elif [[ "$genome_name" == "hg38" ]] ; then
	# transcriptome_dir="/gpfs/data/sequence/cellranger-refdata/refdata-gex-GRCh38-2020-A"
	transcriptome_dir="/gpfs/data/sequence/cellranger-refdata/refdata-gex-GRCh38-2024-A"
elif [[ "$genome_name" == "GRCh38" ]] ; then
	# transcriptome_dir="/gpfs/data/sequence/cellranger-refdata/refdata-gex-GRCh38-2020-A"
	transcriptome_dir="/gpfs/data/sequence/cellranger-refdata/refdata-gex-GRCh38-2024-A"
elif [[ "$genome_name" == "mm10" ]] ; then
	transcriptome_dir="/gpfs/data/sequence/cellranger-refdata/refdata-gex-mm10-2020-A"
elif [[ "$genome_name" == "mm39" ]] ; then
	transcriptome_dir="/gpfs/data/sequence/cellranger-refdata/refdata-gex-GRCm39-2024-A"
elif [[ "$genome_name" == "GRCm39" ]] ; then
	transcriptome_dir="/gpfs/data/sequence/cellranger-refdata/refdata-gex-GRCm39-2024-A"
elif [[ "$genome_name" == "GRCh38_and_mm10" ]] ; then
	transcriptome_dir="/gpfs/data/sequence/cellranger-refdata/refdata-gex-GRCh38-and-mm10-2020-A"
elif [[ "$genome_name" == "GRCh38_and_GRCm39" ]] ; then
	transcriptome_dir="/gpfs/data/sequence/cellranger-refdata/refdata-gex-GRCh38_and_GRCm39-2024-A"
else
	# additional custom genomes
	transcriptome_dir="/gpfs/data/igorlab/ref/${genome_name}/cellranger"
fi

# check that input exists

if [ ! -d "$transcriptome_dir" ] ; then
	echo -e "\n ERROR: genome dir $transcriptome_dir does not exist \n" >&2
	exit 1
fi

if [ -n "$fastq_dir" ] ; then
	# RNA-only command
	if [ ! -d "$fastq_dir" ] ; then
		echo -e "\n ERROR: fastq dir $fastq_dir does not exist \n" >&2
		exit 1
	fi
else
	# RNA and ADT command
	if [ ! -s "$libraries_csv" ] ; then
		echo -e "\n ERROR: libraries csv $libraries_csv does not exist \n" >&2
		exit 1
	fi
	if [ ! -s "$features_csv" ] ; then
		echo -e "\n ERROR: features csv $features_csv does not exist \n" >&2
		exit 1
	fi
fi

# clean up the environment
module purge
module add default-environment
module add cellranger/${module_version}

# display settings
echo " * cellranger:        $(which cellranger) "
echo " * threads:           $threads "
echo " * mem:               $mem "
echo " * transcriptome dir: $transcriptome_dir "
echo " * out dir:           $sample_name_out "
if [ -n "$fastq_dir" ] ; then
	echo " * sample name:       $sample_name_fastq "
	echo " * fastq dir:         $fastq_dir "
	extra_args="--sample $sample_name_fastq --fastqs $fastq_dir"
else
	echo " * libraries csv:     $libraries_csv "
	echo " * features csv:      $features_csv "
	extra_args="--libraries $libraries_csv --feature-ref $features_csv"
fi

echo -e "\n $(date) \n" >&2

# cellranger count command

# transcriptome  Path of folder containing 10x-compatible transcriptome reference
# id             A unique run id and output folder name [a-zA-Z0-9_-]+
# sample         Prefix of the filenames of FASTQs to select
# fastqs         Path of folder created by 10x demultiplexing or bcl2fastq
# libraries      CSV file declaring input library data sources
# feature-ref    Feature reference CSV file, declaring Feature Barcode constructs and associated barcodes

cellranger_cmd="
cellranger count \
--create-bam true \
--localmem $mem \
--localcores $threads \
--transcriptome $transcriptome_dir \
--id $sample_name_out \
$extra_args \
--disable-ui \
"
echo -e "\n CMD: $cellranger_cmd \n"
eval "$cellranger_cmd"

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
