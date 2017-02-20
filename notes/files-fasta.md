# Working with FASTA Files


Index FASTA:
```bash
samtools faidx genome.fa
```

***

Remove empty records (description without sequence):
```bash
awk '$2{print RS}$2' FS='\n' RS=\> ORS= in.fasta > out.fasta
```

***

Remove blank lines:
```bash
sed -i '/^$/d' in.fasta
```

***

Remove problematic characters (they may cause issues with some tools):
```bash
sed -i -e "s/[ ,\(\)\.\/\|:=]/_/g" in.fasta
sed -i 's/___/__/g' in.fasta
```

***

Filter FASTA file by sequence length.

Using `awk`:
```bash
# if applicable, convert multi-line FASTA to single-line FASTA
# using awk:
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' file.fa > file.nowrap.fa
# using FASTX-Toolkit:
fasta_formatter -w 0 -i file.fa -o file.nowrap.fa

# filter by sequence size (1000 in this example)
awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 1000 {print ">"$0}' file.nowrap.fa > file.1000.fa
```
Using `faFilter`:
```bash
faFilter -minSize=N -maxSize=N in.fa out.fa
```
`faFilter` obtained from http://hgdownload.soe.ucsc.edu/admin/exe/
