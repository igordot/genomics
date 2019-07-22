# Initializing draft genome directory


The requirement for [setting up a reference genome directory](https://github.com/igordot/genomics/blob/master/workflows/ref-genome-init.md) is having two basic files:

* `genome.fa` - genome sequence in FASTA format
* `genes.gtf` - gene annotations in GTF format

Properly formatted sequence and annotation files are readily available for commonly analyzed genomes. Less popular and 
draft genomes are less standardized and can be more difficult to work with.

## Saimiri boliviensis boliviensis (Bolivian squirrel monkey) example

As of September 2016, the best reference for the Bolivian squirrel monkey is the first preliminary assembly SaiBol1.0 
(GCA_000235385.1), provided by the Broad Institute in October 2011. The assembly comprises 2685 top level sequences, all 
of which are unplaced scaffolds (from 151,413 contigs).

When working with a new genome, Ensembl is usually a good place to start as it contains well-formatted reference files for 
many species. You can find the [SaiBol1 genome](http://pre.ensembl.org/Saimiri_boliviensis/Info/Index).

The FASTA and GTF files are available and can be downloaded:

```bash
wget -O genome.ensembl.pre.fa.gz ftp://ftp.ensembl.org/pub/pre/fasta/dna/saimiri_boliviensis/Saimiri_boliviensis.SaiBol1.0.dna_rm.toplevel.fa.gz
wget -O genes.ensembl.pre.gtf.gz ftp://ftp.ensembl.org/pub/pre/gtf/saimiri_boliviensis/SaiBol1.0.genes.gtf.gz
```

Once you have the reference files, it's a good idea to spot-check them.

```bash
zcat genome.ensembl.pre.fa.gz | grep ">" | head
```

Output:

```
>scaffold:SaiBol1.0:JH378105.1:1:72162052:1 scaffold JH378105.1
>scaffold:SaiBol1.0:JH378106.1:1:71252344:1 scaffold JH378106.1
>scaffold:SaiBol1.0:JH378107.1:1:58292249:1 scaffold JH378107.1
>scaffold:SaiBol1.0:JH378108.1:1:54856640:1 scaffold JH378108.1
>scaffold:SaiBol1.0:JH378109.1:1:50794693:1 scaffold JH378109.1
>scaffold:SaiBol1.0:JH378110.1:1:49021937:1 scaffold JH378110.1
>scaffold:SaiBol1.0:JH378111.1:1:46157118:1 scaffold JH378111.1
>scaffold:SaiBol1.0:JH378112.1:1:45331107:1 scaffold JH378112.1
>scaffold:SaiBol1.0:JH378113.1:1:44311105:1 scaffold JH378113.1
>scaffold:SaiBol1.0:JH378114.1:1:44255708:1 scaffold JH378114.1
```

```bash
zcat genes.ensembl.pre.gtf.gz | head
```

Output:

```
JH378796.1	protein_coding	exon	805	910	.	+	.	 gene_id "ENSP00000271139_1"; transcript_id "ENSP00000271139_1"; exon_number "1"; gene_biotype "protein_coding";
JH378796.1	protein_coding	CDS	805	910	.	+	0	 gene_id "ENSP00000271139_1"; transcript_id "ENSP00000271139_1"; exon_number "1"; gene_biotype "protein_coding"; protein_id "ENSP00000271139_1";
JH378796.1	protein_coding	exon	2580	3055	.	+	.	 gene_id "ENSP00000271139_1"; transcript_id "ENSP00000271139_1"; exon_number "2"; gene_biotype "protein_coding";
JH378796.1	protein_coding	CDS	2580	3055	.	+	2	 gene_id "ENSP00000271139_1"; transcript_id "ENSP00000271139_1"; exon_number "2"; gene_biotype "protein_coding"; protein_id "ENSP00000271139_1";
JH378584.1	protein_coding	exon	731	835	.	+	.	 gene_id "ENSP00000250416_1"; transcript_id "ENSP00000250416_1"; exon_number "1"; gene_biotype "protein_coding";
JH378584.1	protein_coding	CDS	731	835	.	+	0	 gene_id "ENSP00000250416_1"; transcript_id "ENSP00000250416_1"; exon_number "1"; gene_biotype "protein_coding"; protein_id "ENSP00000250416_1";
JH378584.1	protein_coding	exon	2441	2603	.	+	.	 gene_id "ENSP00000250416_1"; transcript_id "ENSP00000250416_1"; exon_number "2"; gene_biotype "protein_coding";
JH378584.1	protein_coding	CDS	2441	2603	.	+	0	 gene_id "ENSP00000250416_1"; transcript_id "ENSP00000250416_1"; exon_number "2"; gene_biotype "protein_coding"; protein_id "ENSP00000250416_1";
JH378584.1	protein_coding	exon	3155	3293	.	+	.	 gene_id "ENSP00000250416_1"; transcript_id "ENSP00000250416_1"; exon_number "3"; gene_biotype "protein_coding";
JH378584.1	protein_coding	CDS	3155	3293	.	+	2	 gene_id "ENSP00000250416_1"; transcript_id "ENSP00000250416_1"; exon_number "3"; gene_biotype "protein_coding"; protein_id "ENSP00000250416_1";
```

There are several issues:

* FASTA contig names have spaces. This is usually fine, but will cause errors with some tools.
* GTF contig names (first column) do not match FASTA contig names (the full line after `>`).
* GTF records have `gene_id`, but not `gene_name`. Most biologists will want to know the gene names.

The first two problems can be solved with a few command line commands. The last one is more complicated.

Let's try NCBI next.

NCBI has a lot of databases, so it can be difficult to navigate. It also has the 
[SaiBol1 genome](https://www.ncbi.nlm.nih.gov/genome/6907), which is also based GCA_000235385.1 like Ensembl.

There is no GTF, but there is a GFF.

```bash
wget -O genome.ncbi.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/Saimiri_boliviensis/CHR_Un/39432_ref_SaiBol1.0_chrUn.fa.gz
wget -O genes.ncbi.gff3.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/Saimiri_boliviensis/GFF/ref_SaiBol1.0_scaffolds.gff3.gz
```

There are actually two GFF files: `scaffolds.gff3.gz` and `top_level.gff3.gz`, but they are identical based on both file 
size and `md5sum`.

```
$ zcat genome.ncbi.fa.gz | grep ">" | head
>gi|395726353|ref|NW_003943604.1| Saimiri boliviensis boliviensis isolate 3227 unplaced genomic scaffold, SaiBol1.0 scaffold00001, whole genome shotgun sequence
>gi|395726233|ref|NW_003943605.1| Saimiri boliviensis boliviensis isolate 3227 unplaced genomic scaffold, SaiBol1.0 scaffold00002, whole genome shotgun sequence
>gi|395726111|ref|NW_003943606.1| Saimiri boliviensis boliviensis isolate 3227 unplaced genomic scaffold, SaiBol1.0 scaffold00003, whole genome shotgun sequence
>gi|395725977|ref|NW_003943607.1| Saimiri boliviensis boliviensis isolate 3227 unplaced genomic scaffold, SaiBol1.0 scaffold00004, whole genome shotgun sequence
>gi|395725897|ref|NW_003943608.1| Saimiri boliviensis boliviensis isolate 3227 unplaced genomic scaffold, SaiBol1.0 scaffold00005, whole genome shotgun sequence
>gi|395725895|ref|NW_003943609.1| Saimiri boliviensis boliviensis isolate 3227 unplaced genomic scaffold, SaiBol1.0 scaffold00006, whole genome shotgun sequence
>gi|395725818|ref|NW_003943610.1| Saimiri boliviensis boliviensis isolate 3227 unplaced genomic scaffold, SaiBol1.0 scaffold00007, whole genome shotgun sequence
>gi|395725816|ref|NW_003943611.1| Saimiri boliviensis boliviensis isolate 3227 unplaced genomic scaffold, SaiBol1.0 scaffold00008, whole genome shotgun sequence
>gi|395725734|ref|NW_003943612.1| Saimiri boliviensis boliviensis isolate 3227 unplaced genomic scaffold, SaiBol1.0 scaffold00009, whole genome shotgun sequence
>gi|395725652|ref|NW_003943613.1| Saimiri boliviensis boliviensis isolate 3227 unplaced genomic scaffold, SaiBol1.0 scaffold00010, whole genome shotgun sequence
```

```
$ zcat genes.ncbi.gff3.gz | grep -v "#" | head
NW_003943604.1	RefSeq	region	1	72162052	.	+	.	ID=id0;Dbxref=taxon:39432;Name=Unknown;chromosome=Unknown;gbkey=Src;genome=genomic;isolate=3227;mol_type=genomic DNA;sex=female;sub-species=boliviensis
NW_003943604.1	Gnomon	gene	8363	13782	.	-	.	ID=gene0;Dbxref=GeneID:101049931;Name=LOC101049931;gbkey=Gene;gene=LOC101049931
NW_003943604.1	Gnomon	mRNA	8363	13782	.	-	.	ID=rna0;Parent=gene0;Dbxref=GeneID:101049931,Genbank:XM_010336801.1;Name=XM_010336801.1;gbkey=mRNA;gene=LOC101049931;product=breast cancer type 2 susceptibility protein;transcript_id=XM_010336801.1
NW_003943604.1	Gnomon	exon	13673	13782	.	-	.	ID=id1;Parent=rna0;Dbxref=GeneID:101049931,Genbank:XM_010336801.1;gbkey=mRNA;gene=LOC101049931;product=breast cancer type 2 susceptibility protein;transcript_id=XM_010336801.1
NW_003943604.1	Gnomon	exon	11227	11475	.	-	.	ID=id2;Parent=rna0;Dbxref=GeneID:101049931,Genbank:XM_010336801.1;gbkey=mRNA;gene=LOC101049931;product=breast cancer type 2 susceptibility protein;transcript_id=XM_010336801.1
NW_003943604.1	Gnomon	exon	8363	8619	.	-	.	ID=id3;Parent=rna0;Dbxref=GeneID:101049931,Genbank:XM_010336801.1;gbkey=mRNA;gene=LOC101049931;product=breast cancer type 2 susceptibility protein;transcript_id=XM_010336801.1
NW_003943604.1	Gnomon	CDS	13673	13739	.	-	0	ID=cds0;Parent=rna0;Dbxref=GeneID:101049931,Genbank:XP_010335103.1;Name=XP_010335103.1;gbkey=CDS;gene=LOC101049931;product=breast cancer type 2 susceptibility protein;protein_id=XP_010335103.1
NW_003943604.1	Gnomon	CDS	11227	11475	.	-	2	ID=cds0;Parent=rna0;Dbxref=GeneID:101049931,Genbank:XP_010335103.1;Name=XP_010335103.1;gbkey=CDS;gene=LOC101049931;product=breast cancer type 2 susceptibility protein;protein_id=XP_010335103.1
NW_003943604.1	Gnomon	CDS	8363	8619	.	-	2	ID=cds0;Parent=rna0;Dbxref=GeneID:101049931,Genbank:XP_010335103.1;Name=XP_010335103.1;gbkey=CDS;gene=LOC101049931;product=breast cancer type 2 susceptibility protein;protein_id=XP_010335103.1
NW_003943604.1	Gnomon	gene	16666	24141	.	+	.	ID=gene1;Dbxref=GeneID:101050261;Name=ZAR1L;gbkey=Gene;gene=ZAR1L
```

These reference files aren't perfect either, but at least the gene names are present.

Fix the FASTA contig names so they match the GFF contig names (`NW_*` or `NC_*`):

```bash
zcat genome.ncbi.fa.gz | perl -pe 's/gi\|.*\|(N._.+?)\|.*/\1/g' > genome.fa
```

Convert GFF to GTF using `gffread` (part of Cufflinks suite):

```bash
zcat genes.ncbi.gff3.gz | gffread - -T -o genes.ncbi.gff3.gtf
```

A few (24 out of over 850,000) of the GTF entries do not contain `gene_id` or `gene_name`. Remove those:

```bash
cat genes.ncbi.gff3.gtf | grep "gene_name" > genes.gtf
```

This leaves us with a clean `genome.fa` and `genes.gtf` for [setting up a reference genome directory](https://github.com/igordot/genomics/blob/master/workflows/ref-genome-init.md).
