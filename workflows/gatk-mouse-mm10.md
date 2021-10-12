# Creating GATK mm10 resource bundle


The GATK resource bundle is a collection of standard files for working with human resequencing data.
It contains known SNPs and indels to be used for BaseRecalibrator, RealignerTargetCreator, and IndelRealigner.
This is an attempt to recreate a similar bundle for the mouse genome (UCSC build mm10).

For mouse SNPs, it's possible to use the dbSNP database, which should be comparable to the human version.

Download dbSNP GRCm38 VCF files (each chromosome is in a separate file):

```bash
wget --recursive --no-parent --no-directories \
--accept vcf*vcf.gz \
ftp://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/
```

Delete the non-primary chromosomes if they are not included in the reference FASTA file.

Add "chr" to each chromosome (convert from GRCm38 to mm10 format):

```bash
for vcf in $(ls -1 vcf_chr_*.vcf.gz) ; do
  vcf_new=${vcf/.vcf.gz/.vcf}
  echo $vcf
  zcat $vcf | sed 's/^\([0-9XY]\)/chr\1/' > $vcf_new
  rm -fv $vcf
done
```

Combine all dbSNP VCF files into one:

```bash
# generate parameter string containing all VCF files
vcf_file_string=""
for vcf in $(ls -1 vcf_chr_*.vcf) ; do
  vcf_file_string="$vcf_file_string -V $vcf"
done
echo $vcf_file_string
```

For the `MergeVcfs` script to work it needs a dictonary of the reference genom. This must be done with

```bash
picard CreateSequenceDictionary R=genome.fa O=mm10.dict
```

```bash
# concatenate VCF files
picard MergeVcfs $vcf_file_string O=output_variants.vcf.gz D=/home/yeroslaviz/projects/Guru/WES_b2m/reference/mm10/mm10.dict
```


More recent dbSNP releases include a merged `00-All.vcf.gz` file in addition to the separate chromosome files.
Although it will not need to be concatenated, it will likely still need to be adjusted for GATK compatibility.

For mouse indels, the Sanger Mouse Genetics Programme (Sanger MGP) is probably the best resource.

Download all MGP indels (6/2018 release):

```bash
wget ftp://ftp-mouse.sanger.ac.uk/REL-1807-SNPs_Indels/mgp.v6.merged.norm.snp.indels.sfiltered.vcf.gz \
-O mgp.v6.indels.vcf.gz
```

Filter for passing variants with chr added:

```bash
# adjust header
zcat mgp.v6.indels.vcf.gz | head -1000 | grep "^#" | cut -f 1-8 \
| grep -v "#contig" | grep -v "#source" \
> mgp.v5.indels.pass.chr.vcf
# keep only passing and adjust chromosome name
zcat mgp.v6.indels.vcf.gz | grep -v "^#" | cut -f 1-8 \
| grep -w "PASS" | sed 's/^\([0-9MXY]\)/chr\1/' \
>> mgp.v5.indels.pass.chr.vcf
```

Sort VCF (automatically generated index has to be deleted due to a known bug):

```bash
picard SortVcf VERBOSITY=WARNING \
# SD=/mm10.dict # is not needed anymore.\
I=mgp.v6.indels.pass.chr.vcf \
O=mgp.v6.indels.pass.chr.sort.vcf
rm -fv mgp.v5.indels.pass.chr.sort.vcf.idx
```

Additional info:

* [GATK Resource Bundle](https://software.broadinstitute.org/gatk/download/bundle)
* [What should I use as known variants/sites for running tool X?](http://gatkforums.broadinstitute.org/gatk/discussion/1247/what-should-i-use-as-known-variants-sites-for-running-tool-x)
