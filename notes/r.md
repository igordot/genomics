# R

## General R

Vignette:
```r
# show vignettes for a package
browseVignettes(package = "package")
# get vignette
vignette("topic")
```

Get version of package:
```r
packageVersion("package")
```

Working with methods:
```r
# prints source code for method
getMethod(method, "class")
# find method
selectMethod(method, "class")
# show all methods for class
showMethods(classes = "class")
methods(class = "class")
# method help for S3 objects
?"method.class"
```

***

## GenomicRanges

The GenomicRanges package defines general purpose containers for storing and manipulating genomic intervals and variables defined along a genome.

Basics:
```r
library(GenomicRanges)
# create a new GRanges object with one range
g = GRanges("chr1", IRanges(10001, 10100), strand = "+")
# get basic info for ranges
start(g)
end(g)
width(g)
strand(g)
# get metadata columns (additional optional information for ranges)
mcols(g)
# get IRanges
ranges(g)
# get chromosomes for each range
seqnames(g)
# get all chromosomes
seqlevels(g)
```

Intra-range methods (modify each range independently):
* `shift`:	move the ranges by a specific number of base pairs
* `resize`:	resizes to width, keeping start for + and end for -
* `narrow`:	narrows by relative position within range
* `flank`:	returns flanking ranges upstream
* `promoters`:	similar to flank
* `restrict`:	restricts ranges to a start and end position
* `trim`:	trims out of bound ranges
* `+/-`:	add or subtract a fixed amount
* `?"intra-range-methods"`: summarize all intra-range methods

Inter-range methods (comparisons between ranges):
* `reduce`:	merge overlapping ranges to produce a simplified set
* `gaps`:	get gaps between the ranges
* `disjoin`:	break into discrete non-overlapping ranges based on original starts/ends
* `?"inter-range-methods"`: summarize all inter-range methods

Distance methods (compare each range in `x` to `subject`):
* `nearest`:	get an integer vector containing the index of the nearest neighbor range in subject
* `precede`:	get the index of the range in subject that is directly preceded by the range in x
* `follow`:	get the index of the range in subject that is directly followed by the range in x
* `distanceToNearest`:	get the distances to the nearest neighbor in subject (Hits object)
* `distance`: get the distances to the nearest neighbor (integer vector)

Overlaps:
```r
# vector of which x ranges overlap y ranges
x %over% y
# overlaps Hits object
o = findOverlaps(x, y)
# relative to x
queryHits(o)
# relative to y
subjectHits(o)
```

***

Add chr to chromosome names RangedData data structure (from NCBI/Ensembl to UCSC style).
```r
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
