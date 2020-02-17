# Creating a reference genome with exogenous sequences such as GFP 


A standard reference genome is sufficient for most sequencing-based studies, but the experiment may be more complicated.
Studies may involve knock-in or transgenic organisms where the genome sequence is altered.
Since a FASTA file can contain multiple sequences, it's trivial to create a combined one if the exact position in the genome where the foreign sequence is introduced is not relevant.
However, for RNA-seq, you additionally need to modify gene annotations, which is more involved.

Green fluorescent protein (GFP) is a frequently introduced sequence.
Searching for the exact GFP sequence yields many variants, since there are multiple source species of wild-type GFP and
various engineered derivatives, such as yellow fluorescent protein (YFP) or TurboGFP.
This example uses the "enhanced" or "eukaryotic" GFP (EGFP), commonly used in the mammalian expression vectors.

This is the sequence that was used to create a GFP FASTA file `genome.EGFP.fa`:
```
>EGFP
ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAA
GTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGC
TGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAG
CAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTA
CAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGG
ACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAAC
GGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACAC
CCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACG
AGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA
```

To create the corresponding GTF, it's necessary to find the sequence length:
```bash
cat genome.EGFP.fa | grep -v "^>" | tr -d "\n" | wc -c
```

The length is 720. This was used to manually create `genes.EGFP.gtf`:
```
EGFP	unknown	gene	1	720	.	+	.	gene_id "EGFP"; gene_name "EGFP"; gene_biotype "protein_coding";
EGFP	unknown	transcript	1	720	.	+	.	gene_id "EGFP"; transcript_id "EGFP"; gene_name "EGFP"; gene_biotype "protein_coding";
EGFP	unknown	exon	1	720	.	+	.	gene_id "EGFP"; transcript_id "EGFP"; gene_name "EGFP"; gene_biotype "protein_coding";
```

Make sure that the sequence name (column 1) in the GTF matches the one in the FASTA file.

Different tools expect different content from a GTF file.
Using `gene`, `transcript`, and `exon` features seems to be sufficient.
The `gene_biotype` attribute was only added to help with filtering.

Finally, merge the standard genome FASTA and GTF files with the GFP ones:
```bash
cat genome.mm10.fa genome.EGFP.fa > genome.fa
cat genes.mm10.gtf genes.EGFP.gtf > genes.gtf
```

This produces the FASTA and GTF files for the combined reference genome that should be compatible with most GTF-based tools, such as STAR or Cell Ranger.
Pseudoaligners like Kallisto or Salmon build an index from a FASTA formatted file of target sequences, so you can simply append `genome.EGFP.fa` to the cDNA FASTA file.
