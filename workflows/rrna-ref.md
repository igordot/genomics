# Creating ribosomal RNA reference sequence


RNA-seq libraries are typically prepared from total RNA using poly(A) enrichment of the mRNA to remove ribosomal RNAs, but this method fails to capture non-poly(A) transcripts or partially degraded mRNAs.
As an alternative, there are total RNA-seq protocols that require a separate rRNA depletion step.
To test the effectiveness of rRNA depletion, it's a good idea to check rRNA levels in the final RNA-seq library.

Some rRNAs are annotated in the GTF file and show up along with other genes in the final output.
However, rRNA abundance may be substantially underrepresented, since those sequences can fall in the repetitive regions of genome and many tools filter out multi-mapping reads.
Thus, it may be useful to create a separate set of rRNA sequences to align against.
This would also be necessary for tools like [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/), which check the composition of a library by screening it against a set of sequence databases.

A good resource for rRNA sequences is [RNAcentral](https://rnacentral.org/), a database of non-coding RNA from multiple databases such as Rfam and RDP (Ribosomal Database Project).

Download and merge RNAcentral FASTA files:

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/12.0/sequences/rnacentral_species_specific_ids.fasta.gz
```

Convert multi-line sequences to single-line (using `fasta_formatter` from FASTX-Toolkit):

```bash
gzip -cd rnacentral_species_specific_ids.fasta.gz \
  | fasta_formatter -w 0 \
  | gzip \
  > rnacentral.nowrap.fasta.gz
```

Remove empty lines, replace spaces with underscores, and keep just ribosomal sequences:

```bash
gzip -cd rnacentral.nowrap.fasta.gz \
  | sed '/^$/d' \
  | sed 's/\s/_/g' \
  | grep -E -A 1 "ribosomal_RNA|rRNA" \
  | grep -v "^--$" \
  | gzip \
  > rnacentral.ribosomal.nowrap.fasta.gz
```

Set a variable for the species of interest. For example:

```bash
species="homo_sapiens"
species="mus_musculus"
species="drosophila_melanogaster"
```

Extract species-specific ribosomal sequences:

```bash
zcat rnacentral.ribosomal.nowrap.fasta.gz \
  | grep -A 1 -F -i "${species}" \
  | grep -v "^--$" \
  | fasta_formatter -w 80 \
  > rRNA.${species}.fa
```

Index the species-specific FASTA file (not necessary, but will confirm that the FASTA file is valid):

```bash
samtools faidx rRNA.${species}.fa
```

Check the number of sequences per species (there are currently around 6,000 for human and 600 for mouse in RNAcentral):

```bash
wc -l *fai
```

Build the species-specific Bowtie2 index for tools like FastQ Screen:

```bash
bowtie2-build rRNA.${species}.fa rRNA.${species}
```
