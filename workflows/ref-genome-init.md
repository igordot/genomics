# Initializing reference genome directory


Most bioinformatic tools require a reference, but the exact requirements can vary. These are my steps for setting up 
a reference genome directory that can be used for common analysis types (RNA-seq, ChIP-seq, WGS/WES). This workflow is 
inspired by Illumina iGenomes, a collection of sequence and annotation files for commonly analyzed genomes.

The assumption is that you already have two basic files:
 - `genome.fa` - genome sequence in FASTA format
 - `genes.gtf` - gene annotations in GTF format

It's important that sequence names (such as "chr1") are identical in the FASTA and GTF files.

Generate FASTA index `genome.fa.fai` using samtools:
```
samtools faidx genome.fa
```

Generate sequence dictionary `genome.dict` using Picard (used by Picard and GATK):
```
java -Xmx8G -jar ${PICARD_ROOT}/picard.jar CreateSequenceDictionary VERBOSITY=WARNING REFERENCE=genome.fa OUTPUT=genome.dict
```

Generate `chrom.sizes` (used by bedtools and UCSC utilities):
```
cut -f 1,2 genome.fa.fai > chrom.sizes
```

Generate genome BED file using bedtools:
```
bedtools makewindows -g chrom.sizes -w 999999999 | LC_ALL=C sort -k1,1 > genome.bed
```

Generate Bowtie 2 index:
```
mkdir Bowtie2 && cd Bowtie2 && ln -s ../genome.fa && bowtie2-build genome.fa genome && cd ..
```

Generate BWA index:
```
mkdir BWA && cd BWA && ln -s ../genome.fa && bwa index -a bwtsw genome.fa && cd ..
```

Generate gene predictions in refFlat format:
```
gtfToGenePred -genePredExt -geneNameAsName2 genes.gtf refFlat.tmp.txt
paste <(cut -f 12 refFlat.tmp.txt) <(cut -f 1-10 refFlat.tmp.txt) > refFlat.txt
rm -f refFlat.tmp.txt
gzip refFlat.txt
```

Generate STAR index:
```
mkdir STAR
STAR --runMode genomeGenerate --sjdbOverhang 100 --runThreadN 8 \
  --genomeDir STAR --genomeFastaFiles genome.fa --sjdbGTFfile genes.gtf
```
