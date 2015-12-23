# R

Add chr to chromosome names RangedData data structure (from NCBI/Ensembl to UCSC style).
```
ann
# RangedData with 38293 rows and 2 value columns across 51 spaces
#                       space               ranges   |    strand
#                    <factor>            <IRanges>   | <integer>
# ENSMUSG00000090025        1   [3054233, 3054733]   |         1
# ENSMUSG00000064842        1   [3102016, 3102125]   |         1
# ENSMUSG00000051951        1   [3205901, 3671498]   |        -1
 
names(ann) = paste("chr", names(ann), sep="")
 
ann
# RangedData with 38293 rows and 2 value columns across 51 spaces
#                       space               ranges   |    strand
#                    <factor>            <IRanges>   | <integer>
# ENSMUSG00000090025     chr1   [3054233, 3054733]   |         1
# ENSMUSG00000064842     chr1   [3102016, 3102125]   |         1
# ENSMUSG00000051951     chr1   [3205901, 3671498]   |        -1
```
