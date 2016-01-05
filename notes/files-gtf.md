# Working with GTF Files

Some viral and bacterial GTF/GFF files only contain "CDS" feature types (column 3), but certain programs require "exon" featues. Convert "CDS" to "exon", which would be equivalent for that purpose:
```
cat in.gtf | perl -pe 's/\tCDS\t/\texon\t/g' > out.gtf
```

***

Verify the GTF file format and that the genes specified by the GTF file do not violate the rules of gene structure:
```
validate_gtf.pl genes.gtf
```
validate_gtf.pl (by Evan Keibler) obtained from http://mblab.wustl.edu/software.html (there is also a different version included in the Eval package)
