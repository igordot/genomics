#!/bin/bash


##
## 10X Cell Ranger
## Processes Chromium single-cell V(D)J and Gene Expression output (cellranger multi)
## Enables the analysis of multiple library types together (compared to using cellranger vdj and cellranger count separately)
##
## Usage:
## sbatch --job-name=cellranger-${sample} --ntasks=1 --cpus-per-task=17 --mem=128G --time=10:00:00 \
##   --mail-user=${USER}@nyumc.org --mail-type=FAIL,END --export=NONE \
##   --wrap="bash ./scrna-10x-cellranger-multi.sh module_version sample_name config_csv"
##


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")

# check for correct number of arguments
if [ ! $# == 3 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name module_version sample_name \n" >&2
	exit 1
fi

# arguments
module_version=$1
sample_name=$2
config_csv=$(readlink -f "$3")

# settings (16 threads and 64G does not finish with 10h time limit)
threads=16
mem=128

# make the output group-writeable
umask 007

# check that input exists
if [ ! -s "$config_csv" ] ; then
	echo -e "\n ERROR: config $config_csv does not exist \n" >&2
	exit 1
fi

module purge
module add default-environment
module add cellranger/${module_version}

# display settings
echo " * cellranger:        $(which cellranger) "
echo " * threads:           $threads "
echo " * mem:               $mem "
echo " * out dir:           $sample_name "

echo -e "\n $(date) \n" >&2

# cellranger multi command

# id             A unique run id and output folder name [a-zA-Z0-9_-]+
# csv            Path of CSV file enumerating input libraries and analysis parameters
# sample         Prefix of the filenames of FASTQs to select

# The multi config CSV contains both the library definitions and experiment configuration variables.
# It is composed of up to four sections: [gene-expression], [feature], [vdj], and [libraries].
# Template: https://support.10xgenomics.com/multi-config-template.csv

cellranger_cmd="
cellranger multi \
--localmem $mem \
--localcores $threads \
--csv $config_csv \
--id $sample_name \
--disable-ui \
"
echo -e "\n CMD: $cellranger_cmd \n"
$cellranger_cmd

sleep 15

web_summary_html="./${sample_name}/outs/per_sample_outs/${sample_name}/web_summary.html"

# check that output html summary (and probably everything else) exists
if [ ! -s "$web_summary_html" ] ; then
	echo -e "\n ERROR: summary $web_summary_html does not exist \n" >&2
	exit 1
fi

# copy html summary to top level for easy navigation
rsync -tv "$web_summary_html" "./${sample_name}.html"

# delete temp files
rm -rf "./${sample_name}/SC_MULTI_CS"

echo -e "\n $(date) \n"



# end
