#!/usr/bin/env Rscript


'
Description:
  Generate a color-coded gain/loss plot (heatmap-style) from Control-FREEC output.

Usage:
  cnv-freec-heatmap.R <genome> <cnvs_txt> <out_png>

Arguments:
  <genome>    genome build (UCSC-style such as "hg19" or "mm10")
  <cnvs_txt>  Control-FREEC "_CNVs" file with copy number alterations and p-values added by assess_significance.R
  <out_png>   output png image

Options:
  -h, --help        show this screen
' -> doc


# print warnings as they occur
options(warn = 1)

# retrieve the command-line arguments
suppressPackageStartupMessages(library(docopt))
opts = docopt(doc)

# relevent arguments
genome = opts$genome
cnvs_txt = opts$cnvs_txt
cnvs_png = opts$out_png

# check that input file exists
if (!file.exists(cnvs_txt)) stop("file does not exist: ", cnvs_txt)

# load libraries
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(karyoploteR))

# import CNVs table and filter by Wilcoxon p-value
cnvs = read_tsv(cnvs_txt, guess_max = 999999) %>%
  filter(WilcoxonRankSumTestPvalue < 0.05)

# convert ratio table to GRanges
cnvs_gr = GRanges(cnvs)
seqlevelsStyle(cnvs_gr) = "UCSC"

# separate ratios based on amplifications/deletions
cvns_gain = cnvs_gr[cnvs_gr$status == "gain"]
cnvs_loss = cnvs_gr[cnvs_gr$status == "loss"]

# plot
png(cnvs_png, res = 300, width = 15, height = 2, units = "in")
pp = getDefaultPlotParams(plot.type = 4)
pp$data1inmargin = 0
pp$bottommargin = 80
pp$ideogramheight = 20
kp = plotKaryotype(genome = genome, plot.type = 4, ideogram.plotter = NULL, labels.plotter = NULL, plot.params = pp) %>%
  kpAddCytobandsAsLine() %>%
  kpAddChromosomeNames(srt = 90) %>%
  kpRect(data = cvns_gain, y0 = 0, y1 = 1, col = "firebrick2", border = NA) %>%
  kpRect(data = cnvs_loss, y0 = 0, y1 = 1, col = "royalblue4", border = NA)
dev.off()



# end
