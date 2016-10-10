# Basic QC and setup of an Oxford Nanopore run

During a standard MinION run with a single unbarcoded sample, the MinKNOW software writes a FAST5 file for each DNA 
molecule in a local directory `/minion_path/reads/`. These FAST5 files contain aggregated signal measurements and are 
not basecalled.

Basecalling of the FAST5 files is done by the Metrichor that uploads each pre-basecalled FAST5 file from the
`/minion_path/reads/` directory to the cloud-based basecalling service. Uploaded FAST5 files are moved to the 
`/metrichor_path/reads/uploads` directory. Basecalled FAST5 files are downloaded to the `/metrichor_path/reads/downloads` 
directory and then further separated into `pass` and `fail` directories. For a 2D workflow, a read will pass if
basecalling was successful and a 2D read (template, complement, and hairpin) was produced with a mean base quality 
score greater than 9.

Basecalled FAST5 files can be converted to FASTQ format that is more compatible with verious downstream analysis tools.

Create a variable for the pass reads directory:
```
fast5_dir="metrichor_path/downloads/pass"
```

Check the number of passing FAST5 files:
```
find $fast5_dir -name "*.fast5" | wc -l
```

Calculate overall stats (number of reads, mean read length, etc.):
```
poretools stats --type 2D $fast5_dir > reads.2D.stats.txt
```

Determine nucleotide composition:
```
poretools nucdist $fast5_dir > reads.2D.nucdist.txt
```

Generate FASTQ file
```
poretools fastq --min-length 500 --type 2D $fast5_dir | gzip > reads.2D.fastq.gz
```

Generate FASTA file:
```
poretools fasta --min-length 500 --type 2D $fast5_dir > reads.2D.fasta
```

Generate a histogram of read lengths:
```
poretools hist --theme-bw --min-length 0 --max-length 20000 --num-bins 39 --saveas reads.2D.hist.png $fast5_dir
```

Determine read lengths:
```
samtools faidx reads.2D.fasta
cat reads.2D.fasta.fai | grep "_Basecall_2D_2d" | cut -f 1,2 > reads.2D.length.txt
```
