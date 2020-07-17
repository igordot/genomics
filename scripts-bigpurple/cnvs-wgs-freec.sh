#!/bin/bash


##
## WGS copy number variant analysis using Control-FREEC (with optional matched normal).
## Based on SNS WES version.
##
## Usage:
## sbatch --job-name=cnvs-wgs-${sample} --nodes=1 --ntasks=1 --cpus-per-task=5 --mem=100G --time=5:00:00 \
## --mail-user=${USER}@nyulangone.org --mail-type=FAIL,REQUEUE --export=NONE \
## --wrap="bash ./cnvs-wgs-freec.sh project_dir genome_build sample_name bam [control_bam] [window_size]"
##


#########################


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ $# -lt 4 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir genome_build sample_name bam [control_bam] [window_size] \n" >&2
	exit 1
fi

# arguments
proj_dir=$(readlink -f "$1")
genome_build="$2"
sample_t="$3"
bam_t=$(readlink -f "$4")
bam_n="$5"
win_size="$6"


#########################


# check if control sample and/or window size are specified

if [ -n "$win_size" ] ; then
	# both control sample and window size are specified
	bam_n=$(readlink -f "$bam_n")
	win_size_label="paired-${win_size}"
elif [ -n "$bam_n" ] ; then
	# either control sample or window size are specified
	if [ -e "$bam_n" ] ; then
		bam_n=$(readlink -f "$bam_n")
		win_size=""
		win_size_label="paired-auto"
	else
		win_size="$bam_n"
		bam_n=""
		win_size_label="$win_size"
	fi
else
	# no control sample or window size are specified
	win_size_label="auto"
fi

# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: PROJ DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$bam_t" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam_t DOES NOT EXIST \n" >&2
	exit 1
fi

if [ -n "$bam_n" ] && [ ! -s "$bam_n" ] ; then
	echo -e "\n $script_name ERROR: CONTROL BAM $bam_n DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# settings and files

sample="${sample_t}"

cnv_freec_dir="${proj_dir}/CNV-FREEC-${win_size_label}"
mkdir -p "$cnv_freec_dir"

# need a separate directory for each sample since some auto-generated files have identical filenames
sample_freec_logs_base_dir="${proj_dir}/logs-${segment_name}-${win_size_label}"
sample_freec_logs_dir="${sample_freec_logs_base_dir}/${sample_t}"
mkdir -p "$sample_freec_logs_dir"

summary_csv="${sample_freec_logs_base_dir}/${sample}.${segment_name}.csv"

config_txt="${sample_freec_logs_dir}/config.txt"

out_base_sample="${sample_freec_logs_dir}/$(basename $bam_t)"
# out_base_control="${sample_freec_logs_dir}/$(basename $bam_n)"
fixed_base="${cnv_freec_dir}/${sample_t}"

cpn_sample="${out_base_sample}_sample.cpn"
cpn_control="${out_base_control}_control.cpn"

cnvs_original="${out_base_sample}_CNVs"
ratio_original="${out_base_sample}_ratio.txt"
info_original="${out_base_sample}_info.txt"

# repeated later again based on resolution
cnvs_fixed="${fixed_base}.CNVs.txt"
ratio_fixed="${fixed_base}.ratio.txt"
graph_fixed="${fixed_base}.png"

# unload all loaded modulefiles
module purge
module add default-environment


#########################


# genome-specific settings

if [[ "$genome_build" == "hg19" ]] ; then
	chr_files_dir="/gpfs/data/igorlab/ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/"
	chr_len_file="/gpfs/data/igorlab/ref/hg19/genome.fa.fai"
	gem="/gpfs/data/igorlab/ref/hg19/FREEC/out100m2_hg19.gem"
elif [[ "$genome_build" == "mm10" ]] ; then
	chr_files_dir="/gpfs/data/igorlab/ref/iGenomes/Mus_musculus/UCSC/mm10/Sequence/Chromosomes/"
	chr_len_file="/gpfs/data/igorlab/ref/mm10/genome.fa.fai"
	gem="/gpfs/data/igorlab/ref/mm10/FREEC/out100m4_mm10.gem"
else
	echo -e "\n $script_name ERROR: UNSUPPORTED GENOME \n" >&2
	exit 1
fi


#########################


# skip if output exits already

if [ -s "$cpn_sample" ] || [ -s "$cnvs_original" ] || [ -s "$cnvs_fixed" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 0
fi


#########################


# create config

# either coefficientOfVariation or window must be specified for whole genome sequencing data
if [ -n "$win_size" ] ; then
	win_size_config="window = $win_size"
else
	win_size_config="coefficientOfVariation = 0.05"
fi

# config for control dataset

if [ -n "$bam_n" ] ; then
control_config="
[control]
mateFile = $bam_n
inputFormat = BAM
mateOrientation = FR
"
else
	control_config=""
fi

config_contents="

# whole genome sequencing config
# based on: https://github.com/BoevaLab/FREEC/blob/master/data/config_WGS.txt


[general]

# output directory
outputDir = .

# number of threads
maxThreads = 4

# path to sambamba (used only to read .BAM files)
sambamba = /gpfs/data/igorlab/software/sambamba/sambamba-0.7.0

# file with chromosome lengths (fa.fai accepted starting from v9.3)
chrLenFile = $chr_len_file

# path to the directory with chromosomes fasta files
# necessary to calculate a GC-content profile if a control dataset and GC-content profile are not available
chrFiles = $chr_files_dir

# information about mappable positions (GEM output)
gemMappabilityFile = $gem

# use a mappability profile to correct read counts (provided with gemMappabilityFile)
uniqueMatch = TRUE

# genome ploidy
# you can set different values and Control-FREEC will select the one that explains most observed CNAs
ploidy = 2

# sample sex
# sex=XY will not annotate one copy of chr X and Y as a losssex=XY
sex = XY

# either coefficientOfVariation or window must be specified for whole genome sequencing data
# for whole exome sequencing: window=0
$win_size_config

# set to 1 or 2 to correct the Read Count (RC) for GC-content bias and low mappability
# Default (WGS): 0
# Default (WES): 1 (≥ v9.5) and 0 (< v9.5)
# forceGCcontentNormalization = 1

# degree of polynomial
# Default: 3&4 (GC-content based normalization, WGS) or 1 (control-read-count-based normalization, WES)
# degree = 1

# desired behavior in the ambiguous regions
# 4: make a separate fragment of this unknown region and do not assign any copy number to this region at all
# breakPointType = 4

# positive value of threshold for segmentation of normalized profiles
# Default: 0.8 (for WGS)
breakPointThreshold = 0.8

# threshold on the minimal number of reads per window in the control sample
# recommended value >=50 for for exome data
# readCountThreshold = 10

# additional output in BedGraph format for the UCSC genome browser
BedGraphOutput = TRUE


[sample]

# file with mapped reads (can be single end reads, mate-pairs or paired-end reads)
mateFile = $bam_t

# format of reads (in mateFile)
# SAM, BAM, pileup, bowtie, eland, arachne, psl (BLAT), BED, Eland
inputFormat = BAM

# format of reads (in mateFile)
# 0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)
mateOrientation = FR

$control_config

"

echo "$config_contents" > "$config_txt"


#########################


# Control-FREEC

# FREEC compiled with GCC 6.1.0 (same GCC must be in the environment when running)
module add gcc/6.1.0
# bedtools to create .pileup files for WES data
module add bedtools/2.27.1
# samtools to create .pileup files (for BAF) (even with sambamba enabled)
module add samtools/1.3

cd "$sample_freec_logs_dir"

freec_dir="/gpfs/data/igorlab/software/FREEC/FREEC-11.5"
freec_bin="${freec_dir}/src/freec"

echo
echo " * FREEC: $(readlink -f $freec_bin) "
echo " * sample T : $sample_t "
echo " * BAM T : $bam_t "
echo " * BAM N : $bam_n "
echo " * window size : $win_size "
echo " * CNVs original: $cnvs_original "
echo " * CNVs fixed: $cnvs_fixed "
echo " * ratio original: $ratio_original "
echo " * ratio fixed: $ratio_fixed "
echo

freec_cmd="$freec_bin -conf $config_txt"
echo -e "\n CMD: $freec_cmd \n"
($freec_cmd)

sleep 30


#########################


# check that output generated

if [ ! -s "$cnvs_original" ] ; then
	echo -e "\n $script_name ERROR: $cnvs_original NOT GENERATED \n" >&2
	exit 1
fi

if [ ! -s "$ratio_original" ] ; then
	echo -e "\n $script_name ERROR: $ratio_original NOT GENERATED \n" >&2
	exit 1
fi


#########################


# clean up

# delete raw copy number profiles
rm -v "$cpn_sample"

if [ -s "$cpn_control" ] ; then
	rm -v "$cpn_control"
fi


#########################


# get resolution and add to output names

# get resolution
res=$(cat "$info_original" | grep "Window" | cut -f 2)

# adjust output file names
cnvs_fixed="${fixed_base}.${res}.CNVs.txt"
ratio_fixed="${fixed_base}.${res}.ratio.txt"
graph_fixed="${fixed_base}.${res}.png"

echo
echo " * res: $res "
echo " * CNVs fixed: $cnvs_fixed "
echo " * ratio fixed: $ratio_fixed "
echo " * graph fixed: $graph_fixed "
echo


#########################


# post-processing

module add r/3.6.1

echo
echo " * R: $(readlink -f $(which R)) "
echo " * R version: $(R --version | head -1) "
echo " * Rscript: $(readlink -f $(which Rscript)) "
echo " * Rscript version: $(Rscript --version 2>&1) "
echo

# required libraries: rtracklayer

# add p-value to the predicted CNVs
# add Wilcoxon test and Kolmogorov-Smirnov test p-values to _CNVs file (also add headers to the columns)
freec_asses_sig_cmd="cat ${freec_dir}/scripts/assess_significance.R | R --slave --args $cnvs_original $ratio_original"
echo -e "\n CMD: $freec_asses_sig_cmd \n"
eval "$freec_asses_sig_cmd"

# add "chr" to CNV table chromosomes
cat "${cnvs_original}.p.value.txt" | sed 's/^\([0-9XY]\)/chr\1/' | LC_ALL=C sort -k1,1 -k2,2n | uniq > "$cnvs_fixed"

# visualize normalized copy number profile with predicted CNAs as well as BAF profile by running makeGraph.R
freec_makegraph_cmd="cat ${freec_dir}/scripts/makeGraph.R | R --slave --args 2 $ratio_original"
echo -e "\n CMD: $freec_makegraph_cmd \n"
eval "$freec_makegraph_cmd"


#########################


# fix some of the names

mv -v "$ratio_original" "$ratio_fixed"
mv -v "${ratio_original}.png" "$graph_fixed"


#########################


# check that output generated

if [ ! -s "$cnvs_fixed" ] ; then
	echo -e "\n $script_name ERROR: CNVs $cnvs_fixed NOT GENERATED \n" >&2
	exit 1
fi


#########################


# summary

# ratios and predicted copy number alterations for each window
num_bins=$(cat "$ratio_fixed" | grep -v 'MedianRatio' | wc -l)
echo "num bins: $num_bins"

# predicted copy number alterations
num_cnas=$(cat "$cnvs_fixed" | grep -v 'uncertainty' | wc -l)
echo "num CNAs: $num_cnas"

# header for summary file
echo "#SAMPLE,res,bins,CNAs" > "$summary_csv"

# summarize log file
echo "${sample},${res},${num_bins},${num_cnas}" >> "$summary_csv"

sleep 30

# combine all sample summaries
cat ${sample_freec_logs_base_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq \
> "${proj_dir}/summary.${segment_name}-${win_size_label}.csv"


#########################


# annotate

annot_cmd="bash /gpfs/data/igorlab/public/sns/segments/annot-regions-annovar.sh $proj_dir $sample $cnvs_fixed"
echo -e "\n CMD: $annot_cmd \n"
($annot_cmd)


#########################



# end
