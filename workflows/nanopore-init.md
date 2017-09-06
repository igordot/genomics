# Processing of Oxford Nanopore Technologies (ONT) data


## Sequencing

There are multiple sequencing protocols. The basic one is 1D, which is analogous to Illumina's single-read option where each DNA fragment is sequenced once. For higher accuracy, there was a 2D workflow where each fragment will generate both a template and complement reads (sepated by hairpin). In May 2017, ONT replaced the 2D system (part of a legal dispute since Pacific Biosciences patented the hairpin approach) with 1D^2 or "1D squared". 

During a standard MinION run with a single unbarcoded sample, the MinKNOW software writes a FAST5 file for each DNA molecule in a local directory. These FAST5 files contain aggregated signal measurements and may be basecalled.

FAST5 files overview from [Simpson Lab blog](http://simpsonlab.github.io/2017/09/06/nanopolish-v0.8.0/):
> Oxford Nanopore’s sequencers measure the disruption in electric current caused by single-stranded DNA moving through the nanopore. The device samples the current six thousand times per second and writes the samples to a FAST5 file. We refer to these measurements as “the raw signal” or “the raw samples”, or simply “the raw”. For the past three years nanopore basecallers have converted the raw samples into segments called “events”, with the boundaries between events roughly corresponding to movements of DNA through the pore. After the samples are segmented into events, the basecaller predicts which k-mer was in the pore when the samples for each event were taken. The basecalling results are stored in a new FAST5 file that has a table containing every event and its k-mer label.

## Basecalling

MinKNOW can output basecalled and non-basecalled FAST5 files, but MinKNOW basecalled files may cause issues downstream. ONT also offers [Albacore](https://community.nanoporetech.com/protocols/albacore-offline-basecalli/v/abec_2003_v1_revx_29nov2016) for local offline basecalling.

Albacore needs Python 3.4+ and can by installed using `pip` (it's probably best to use `virtualenv`):
```bash

# create a new environment in your virtualenv directory
cd /virtualenv_path/
pyvenv albacore

# activate virtualenv
source /virtualenv_path/albacore/bin/activate

# upgrade pip
pip install --upgrade pip

# install albacore from .whl (from ONT website)
pip install /path/ont_albacore-x.x.x.whl

# deactivate virtualenv
deactivate
```

Perform basecalling:
```bash
# it may be necessary to unset PYTHONPATH if it was set
unset PYTHONPATH

# activate virtualenv
source /virtualenv_path/albacore/bin/activate

# check available flowcells and kits (must be specified for basecalling)
read_fast5_basecaller.py --list_workflows

# run basecaller (use & to run the process in the background)
read_fast5_basecaller.py \
--worker_threads 8 \
--recursive \
--save_path ./albacore-out \
--output_format fast5 \
--input ./fast5 \
--flowcell FLO-MIN000 \
--kit SQK-NSK000 \
&

# check the number of processed reads
cat ./albacore-out/pipeline.log | grep "Finished" | wc -l

# deactivate virtualenv
deactivate
```

## Data extraction

Basecalled FAST5 files can be converted to FASTQ format that is more compatible with verious downstream analysis tools.

[Poretools](http://poretools.readthedocs.io/) is a popular tool for extracting data and information from FAST5 files.

Create a variable for the basecalled reads directory:
```bash
fast5_dir="/albacore_path/workspace"
```

Calculate overall stats (number of reads, mean read length, etc.):
```bash
poretools stats --type 2D $fast5_dir > reads.2D.stats.txt
```

Determine nucleotide composition:
```bash
poretools nucdist $fast5_dir > reads.2D.nucdist.txt
```

Generate gzipped FASTQ file:
```bash
poretools fastq --min-length 500 --type 2D $fast5_dir | gzip > reads.2D.fastq.gz
```

Generate FASTA file:
```bash
poretools fasta --min-length 500 --type 2D $fast5_dir > reads.2D.fasta
```

Generate a histogram of read lengths:
```bash
poretools hist --theme-bw --min-length 0 --max-length 20000 --num-bins 39 --saveas reads.2D.hist.png $fast5_dir
```

Determine read lengths:
```bash
samtools faidx reads.2D.fasta
cat reads.2D.fasta.fai | grep "_Basecall_2D_2d" | cut -f 1,2 > reads.2D.length.txt
```

## Assembly

[Canu](https://canu.readthedocs.io/) is a common de novo assembler for Nanopore long reads. It performs error correction, but additional polishing is helpful. [Nanopolish](https://github.com/jts/nanopolish) can calculate an improved consensus sequence for a draft genome assembly.

If you have both Illumina and Nanopore data, then [SPAdes](http://bioinf.spbau.ru/spades) is a good option for hybrid assembly. SPAdes will use Nanopore reads for gap closure and repeat resolution.

## Additional info

PoreCamp is a training bootcamp based around Oxford Nanopore MinION sequencing that provides great [tutorials](https://porecamp.github.io/2017/) for basic processing of the data.
