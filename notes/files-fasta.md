# Working with FASTA Files

Index FASTA:
```
samtools faidx genome.fa
```

***

Remove empty records (description without sequence):
```
awk '$2{print RS}$2' FS='\n' RS=\> ORS= in.fasta > out.fasta
```

***

Remove blank lines:
```
sed -i '/^$/d' in.fasta
```

***

Remove problematic characters (may cause issues with some tools):
```
sed -i -e "s/[ ,\(\)\.\/\|:=]/_/g" in.fasta
sed -i 's/___/__/g' in.fasta
```

***

Filter FASTA file by sequence length:
```
# FASTA file must not be wrapped (each sequence should appear on a single line)

# if necessary, convert multi-line FASTA to single-line FASTA
# using awk:
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' file.fa > file.nowrap.fa
# using FASTX-Toolkit:
fasta_formatter -w 0 -i file.fa -o file.nowrap.fa

# filter by sequence size (1000 in this example)
awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 1000 {print ">"$0}' file.nowrap.fa > file.1000.fa
```
