#!/usr/bin/env Rscript


'
Description:
  Plot genome-wide Control-FREEC copy number analysis results with proportional chromosomes in a single line.

Usage:
  cnv-freec-genome-plot.R <genome> <sample_name> <ratio_txt> <out_png>

Arguments:
  <genome>       genome build (UCSC-style such as "hg19" or "mm10")
  <sample_name>  sample name
  <ratio_txt>    Control-FREEC "_ratio.txt" file with ratios and predicted copy number alterations for each window
  <out_png>      output png image

Options:
  -h, --help        show this screen
' -> doc


# output width
options(width = 120)
# print warnings as they occur
options(warn = 1)
# default type for the bitmap devices such as png (should default to "cairo")
options(bitmapType = "cairo")

# retrieve the command-line arguments
suppressPackageStartupMessages(library(docopt))
opts = docopt(doc)

# relevent arguments
genome = opts$genome
sample_name = opts$sample_name
ratio_txt = opts$ratio_txt
ratio_png = opts$out_png

# check that input file exists
if (!file.exists(ratio_txt)) stop("file does not exist: ", ratio_txt)

# load libraries
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(karyoploteR))

# ploidy
ploidy = 2
# maximum copy number level to plot (to avoid very high values)
max_cn = 8

# import ratio table, remove uncertain regions, and cap copy numbers at defined maximum value
ratio =
  read_tsv(ratio_txt, guess_max = 999999, progress = FALSE) %>%
  filter(CopyNumber >= 0) %>%
  mutate(chr = Chromosome, start = Start, end = Start) %>%
  mutate(Ratio = Ratio * ploidy) %>%
  mutate(CopyNumber = ifelse(CopyNumber > max_cn, max_cn, CopyNumber)) %>%
  mutate(Ratio = ifelse(Ratio > max_cn, max_cn, Ratio)) %>%
  select(chr, start, end, Ratio, CopyNumber)

# convert ratio table to GRanges
ratio_gr = GRanges(ratio)
seqlevelsStyle(ratio_gr) = "UCSC"

# separate ratios based on amplifications/deletions
ratio_filtered = ratio_gr[ratio_gr$Ratio > 0]
ratio_norm = ratio_filtered[ratio_filtered$CopyNumber == ploidy]
ratio_amp = ratio_filtered[ratio_filtered$CopyNumber > ploidy]
ratio_del = ratio_filtered[ratio_filtered$CopyNumber < ploidy]

# plot
png(ratio_png, res = 300, width = 15, height = 3, units = "in")
pp = getDefaultPlotParams(plot.type = 4)
pp$data1inmargin = 0
pp$bottommargin = 50
pp$ideogramheight = 20
kp =
  plotKaryotype(genome = genome, plot.type = 4, ideogram.plotter = NULL, labels.plotter = NULL, plot.params = pp) %>%
  kpAxis(ymin = 0, ymax = max_cn, tick.pos = 0:max_cn) %>%
  kpAddCytobandsAsLine() %>%
  kpAddChromosomeNames(srt = 45) %>%
  kpPoints(data = ratio_norm, y = ratio_norm$Ratio,
           cex = 0.3, ymin = 0, ymax = max_cn, col = "darkolivegreen2") %>%
  kpPoints(data = ratio_amp, y = ratio_amp$Ratio,
           cex = 0.3, ymin = 0, ymax = max_cn, col = "firebrick2") %>%
  kpPoints(data = ratio_del, y = ratio_del$Ratio,
           cex = 0.3, ymin = 0, ymax = max_cn, col = "royalblue4") %>%
  kpPoints(data = ratio_gr, y = ratio_gr$CopyNumber,
           cex = 0.5, ymin = 0, ymax = max_cn, col = "black") %>%
  kpAddMainTitle(sample_name)
dev.off()



# end
