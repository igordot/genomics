#!/usr/bin/env Rscript

'
Description:
  Remove ambient RNA contamination from 10x Genomics Chromium single-cell RNA-seq data using SoupX.
  Input and output are in the format produced by the Cell Ranger software suite.

Usage:
  scrna-decontaminate-soupx.R <in_dir> <out_dir>

Arguments:
  <in_dir>     input directory
  <out_dir>    output directory (will contain "outs/filtered_feature_bc_matrix")

Options:
  -h, --help   show this screen
' -> doc


# increase output width
options(width = 120)
# print warnings as they occur
options(warn = 1)

# retrieve the command-line arguments
library(docopt)
opts = docopt(doc)

# relevent arguments
in_dir = opts$in_dir
out_dir = opts$out_dir

# check if the input parameters are valid
message("input dir: ", in_dir)
if (!dir.exists(in_dir)) { stop("input dir does not exist") }
message("output dir: ", out_dir)
if (dir.exists(out_dir)) { stop("output dir already exists") }

# load libraries
suppressPackageStartupMessages({
  library(glue)
  library(Matrix)
  library(SoupX)
  library(DropletUtils)
})

# set output directory as working directory
dir.create(out_dir)
if (dir.exists(out_dir)) {
  setwd(out_dir)
} else {
  stop(glue("output dir {out_dir} could not be created"))
}

# log to file
write(glue("analysis: {out_dir}"), file = "create.log", append = TRUE)
write(glue("soupx version: {packageVersion('SoupX')}"), file = "create.log", append = TRUE)

# find "outs" dir (contains "raw_feature_bc_matrix")
cr_outs_dir = list.files(path = in_dir, pattern = "raw_feature_bc_matrix$", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
cr_outs_dir = dirname(cr_outs_dir)
if (length(cr_outs_dir) != 1) stop(glue("no outs directory in {in_dir}"))

# load data and estimate soup profile
sc = SoupX::load10X(cr_outs_dir)

# log the stats
write(glue("counts matrix cells: {ncol(sc$toc)}"), file = "create.log", append = TRUE)
write(glue("counts matrix genes: {nrow(sc$toc)}"), file = "create.log", append = TRUE)
in_umis = Matrix::colSums(sc$toc)
write(glue("unfiltered mean UMIs per cell: {round(mean(in_umis), 3)}"), file = "create.log", append = TRUE)
write(glue("unfiltered median UMIs per cell: {median(in_umis)}"), file = "create.log", append = TRUE)
in_detected_genes = Matrix::colSums(sc$toc > 0)
write(glue("unfiltered min genes per cell: {min(in_detected_genes)}"), file = "create.log", append = TRUE)
write(glue("unfiltered max genes per cell: {max(in_detected_genes)}"), file = "create.log", append = TRUE)
write(glue("unfiltered mean genes per cell: {round(mean(in_detected_genes), 3)}"), file = "create.log", append = TRUE)
write(glue("unfiltered median genes per cell: {median(in_detected_genes)}"), file = "create.log", append = TRUE)

# estimate the level of background contamination (represented as rho)
# creates a plot showing the density of estimates
png("qc.soupx.estimates.png", res = 300, width = 8, height = 5, units = "in")
sc = SoupX::autoEstCont(sc)
dev.off()

write(glue("estimated rho: {sc$fit$rhoEst}"), file = "create.log", append = TRUE)

# clean the data
soupx_out = SoupX::adjustCounts(sc)
dim(soupx_out)

# log the stats
out_umis = Matrix::colSums(soupx_out)
write(glue("decontaminated mean UMIs per cell: {round(mean(out_umis), 3)}"), file = "create.log", append = TRUE)
write(glue("decontaminated median UMIs per cell: {median(out_umis)}"), file = "create.log", append = TRUE)
out_detected_genes = Matrix::colSums(soupx_out > 0)
write(glue("decontaminated min genes per cell: {min(out_detected_genes)}"), file = "create.log", append = TRUE)
write(glue("decontaminated max genes per cell: {max(out_detected_genes)}"), file = "create.log", append = TRUE)
write(glue("decontaminated mean genes per cell: {round(mean(out_detected_genes), 3)}"), file = "create.log", append = TRUE)
write(glue("decontaminated median genes per cell: {median(out_detected_genes)}"), file = "create.log", append = TRUE)

# write count data in the 10x format (path must not exist)
dir.create("./outs")
DropletUtils::write10xCounts(x = soupx_out, path = "./outs/filtered_feature_bc_matrix", version = "3")

# delete Rplots.pdf
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")



# end
