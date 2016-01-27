# Working with FASTQ Files

Sample random reads from a FASTQ (where NNN is the number of reads):
```
seqtk sample -s 100 in.fastq.gz NNN | gzip > out.fastq.gz
```
Seqtk obtained from https://github.com/lh3/seqtk  
Alternative: sample_fastq.py from https://github.com/mojones/random_scripts/blob/master/sample_fastq.py
