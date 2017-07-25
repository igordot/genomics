# RNA-seq differential expression analysis

Basic RNA-seq differential expression analysis in R.

## DESeq2

Load library:
```r
library(DESeq2)
```

Import data:
```r
# import raw counts matrix
dds = DESeqDataSetFromMatrix(counts, coldata, ~ group)
# import SummarizedExperiment
dds = DESeqDataSet(se, ~ group)
# import tximport
dds = DESeqDataSetFromTximport(txi, coldata, ~ group)
```

Analysis:
```r
dds = DESeq(dds)
res = results(dds)
```

***

## edgeR

Load library:
```r
library(edgeR) 
```

Import data (assuming four RNA-Seq libraries in two groups with counts are stored in a tab-delimited text file and gene symbols in a column "gene"):
```r
# import data 
counts = read.delim("counts.txt", row.names = "gene")
group = factor(c(1,1,2,2))
# edgeR stores data in a list-based data object called a DGEList
y = DGEList(counts = counts, group = group)
y = calcNormFactors(y)
design = model.matrix(~ group)
y = estimateDisp(y, design)
```

Perform likelihood ratio tests:
```r
fit = glmFit(y, design)
lrt = glmLRT(fit)
topTags(lrt)
```

***

## limma-voom

```r
library(limma)
design = model.matrix(~ group)
dgel = DGEList(counts)
dgel = calcNormFactors(dgel)
v = voom(dgel, design,plot=FALSE)
fit = lmFit(v, design)
fit = eBayes(fit)
topTable(fit)
```
