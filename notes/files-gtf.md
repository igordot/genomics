# Working with GTF Files


Verify the GTF file format and confirm that the genes specified do not violate the rules of gene structure:
```bash
validate_gtf.pl genes.gtf
```
`validate_gtf.pl` (by Evan Keibler) obtained from http://mblab.wustl.edu/software.html (there is also a different version included in the Eval package)

***

Add gene names to GTF gene IDs to make them more readable (merge `gene_id` and `gene_name`):
```bash
cat genes.gtf \
| grep "transcript_id" \
| perl -pe 's/(gene_id "(.+?)"; )(.*)(gene_name "(.+?)"; )/gene_id "\5:\2"; \3 \4/g' \
> genes.name-id.gtf
```

***

Some viral and bacterial GTF files only contain `CDS` features (column 3), but some tools require `exon` features.
Convert `CDS` to `exon`, which would be equivalent for that purpose:
```bash
cat in.gtf | perl -pe 's/\tCDS\t/\texon\t/g' > out.gtf
```
