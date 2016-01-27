# Working with chrom.sizes Files

Generate a chrom.sizes file (from an indexed FASTA file):
```
cut -f 1,2 genome.fa.fai > chrom.sizes
```

Alternatives:
* fetchChromSizes (http://hgdownload.soe.ucsc.edu/admin/exe/)
