# Microarray differential expression analysis

Basic microarray differential expression analysis in R using `limma`.

```r
library(affy)
library(limma)
sample_info = read.AnnotatedDataFrame("samples.csv")
eset = justRMA("/path/to/cel-files", phenoData = sample_info)
design = model.matrix(~ group, pData(eset))
fit = lmFit(eset, design)
efit = eBayes(fit)
topTable(efit, coef = 2)
```
