# Creating a reference genome including GFP


You can use a standard reference genome for most sequencing-based studies, but the experiment may be more complicated.
Studies may involve knock-in or transgenic organisms where the genome sequence is altered.
Since a FASTA file can contain multiple sequences, it's trivial to create a combined one if you don't care about the
exact position in the genome where the foreign sequence is introduced.
However, for RNA-seq, you would need to also modify gene annotations, which is more involved.

This example uses green fluorescent protein (GFP) as the introduced sequence, since that is a common use case.
Searching for the exact GFP sequence yields many variants, since there are multiple source species of wild-type GFP and
various engineered derivatives.
Most of the mammalian expression vectors use the "enhanced" or "eukaryotic" GFP (EGFP).

This is the sequence that was used to create a GFP FASTA file `genome.GFP.fa`:
```
>GFP
ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAA
GTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGC
TGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAG
CAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTA
CAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGG
ACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAAC
GGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACAC
CCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACG
AGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAG
```

To create the corresponding GTF, it's necessary to find the sequence length:
```bash
cat genome.GFP.fa | grep -v "^>" | tr -d "\n" | wc -c
```

The length is 717. This was used to manually create `genes.GFP.gtf`:
```
GFP	unknown	gene	1	717	.	+	.	gene_id "GFP"; gene_name "GFP"; gene_biotype "protein_coding";
GFP	unknown	transcript	1	717	.	+	.	gene_id "GFP"; transcript_id "GFP"; gene_name "GFP"; gene_biotype "protein_coding";
GFP	unknown	exon	1	717	.	+	.	gene_id "GFP"; transcript_id "GFP"; gene_name "GFP"; gene_biotype "protein_coding";
```

Make sure that the sequence name (column 1) in the GTF matches the one in the FASTA file.

Different tools expect different content from a GTF file.
Using `gene`, `transcript`, and `exon` features seems to be sufficient.
The `gene_biotype` attribute was only added to help with filtering.

Finally, merge the standard genome FASTA and GTF files with the GFP ones:
```bash
cat genome.mm10.fa genome.GFP.fa > genome.fa
cat genes.mm10.gtf genes.GFP.gtf > genes.gtf
```

This produces the FASTA and GTF files for the combined reference genome that should be compatible with most tools.
