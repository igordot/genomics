# Converting Files

## GFF to GTF
Convert from GFF to GTF:
```
gffread in.gff -T -o out.gtf
```
gffread program is part of the Cufflinks package

***

## GTF to refFlat
Convert gene annotations from GTF to genePred refFlat format:
```
gtfToGenePred -genePredExt -geneNameAsName2 genes.gtf refFlat.tmp.txt
paste <(cut -f 12 refFlat.tmp.txt) <(cut -f 1-10 refFlat.tmp.txt) > refFlat.txt
rm refFlat.tmp.txt
gzip refFlat.txt
```
Tested with Ensembl GTF file and used for Picard CollectRnaSeqMetrics (must be gzipped).  
gtfToGenePred obtained from http://hgdownload.cse.ucsc.edu/admin/exe/

***

## SRA to FASTQ
Convert SRA (Sequence Read Archives) files to FASTQs:
```
/path/sratoolkit.2.3.4/bin/fastq-dump -v --split-files --gzip file.sra
```
SRA run ID to FASTQ (will download the SRA file and put it in a temp directory):
```
/path/sratoolkit.2.3.4/bin/fastq-dump -v --split-files --gzip SRR0000000
```
NCBI SRA Toolkit obtained from http://eutils.ncbi.nih.gov/Traces/sra/?view=software

***

## DOT to SVG/PNG
Convert DOT database schema file to SVG or PNG:
```
dot -Tsvg file.dot > file.svg
dot -Tpng file.dot > file.png
```
