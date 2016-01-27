# Working with FASTQ Files

Subset a FASTQ (where NNN is the desired number of reads):
```
seqtk sample -s 100 in.fastq.gz NNN | gzip > out.fastq.gz
```
Seqtk obtained from https://github.com/lh3/seqtk  
Alternatives: 
* fastq-tools fastq-sample (http://homes.cs.washington.edu/~dcjones/fastq-tools/fastq-sample.html)
* sample_fastq.py (https://github.com/mojones/random_scripts/blob/master/sample_fastq.py)
