# Creating ribosomal RNA reference sequence


RNA-seq libraries are typically prepared from total RNA using poly(A) enrichment of the mRNA to remove ribosomal RNAs,
but this method fails to capture non-poly(A) transcripts or partially degraded mRNAs.
As an alternative, there are total RNA-seq protocols that require a separate rRNA depletion step.
To test the effectiveness of rRNA depletion, it's a good idea to check rRNA levels in the final RNA-seq library.

Some rRNAs should be annotated in the GTF file and show up along other genes in the final output,
but they can also fall in the repetitive regions of genome, resulting in alignment-related issues.
Some tools may filter out multi-mapping reads, causing rRNA abundance to be underrepresented.
Thus, it may be useful to create a separate rRNA reference sequence to align against.
This would also be necessary for tools like
[FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
that allow screening a library against a set of sequence databases to check its composition.

A good resource for rRNA sequences is
[Rfam](http://rfam.xfam.org/), a database of non-coding RNA.

Download and merge Rfam FASTA files:
```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/12.0/fasta_files/*
cat RF*.gz > Rfam.all.fa.gz
```

Convert a multi-line FASTA file to a single-line FASTA file:
```bash
zcat Rfam.all.fa.gz \
| fasta_formatter -w 0 \
| gzip \
> Rfam.nowrap.fa.gz
```

Remove empty lines, replace spaces with underscores, and keep just "ribosomal" sequences:
```bash
zcat Rfam.all.nowrap.fa.gz \
| sed '/^$/d' \
| sed 's/\s/_/g' \
| awk 'BEGIN {RS=">"; ORS="";} tolower($1) ~ /ribosomal_rna/ {print ">"$0}' \
| gzip \
> Rfam.rRNA.nowrap.fa.gz
```

Set a variable for species of interest. For example:
```bash
species="homo_sapiens"
species="mus_musculus"
species="drosophila_melanogaster"
species="mycobacterium"
species="mycobacterium_tuberculosis"
```

Extract species-specific ribosomal sequences:
```bash
zcat Rfam.rRNA.nowrap.fa.gz \
| grep -A 1 -i "${species}" \
| grep -v '^--$' \
| fasta_formatter -w 80 \
> rRNA.${species}.fa
```

Index the species-specific FASTA file:
```bash
samtools faidx rRNA.${species}.fa
```

Build Bowtie2 index:
```bash
bowtie2-build rRNA.${species}.fa rRNA.${species}
```
