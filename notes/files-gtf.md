# Working with GTF Files

Some viral and bacterial GTF/GFF files only contain "CDS" feature types (column 3), but certain programs require "exon" featues. Convert "CDS" to "exon", which would be equivalent for that purpose:
```
cat in.gtf | perl -pe 's/\tCDS\t/\texon\t/g' > out.gtf
```
