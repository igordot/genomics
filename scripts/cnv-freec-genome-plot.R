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
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(karyoploteR)
  library(scales)
})

# ploidy
ploidy = 2
# maximum copy number level to plot (to avoid very high values)
max_cn = 6

# import ratio table, remove uncertain regions, and cap copy numbers at defined maximum value
ratio =
  read_tsv(ratio_txt, guess_max = 999999, show_col_types = FALSE, progress = FALSE) %>%
  dplyr::filter(CopyNumber >= 0) %>%
  dplyr::mutate(chr = Chromosome, start = Start, end = Start) %>%
  dplyr::mutate(Ratio = Ratio * ploidy) %>%
  dplyr::mutate(CopyNumber = ifelse(CopyNumber > max_cn, max_cn, CopyNumber)) %>%
  dplyr::mutate(Ratio = ifelse(Ratio > max_cn, max_cn, Ratio)) %>%
  dplyr::select(chr, start, end, Ratio, CopyNumber)

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
kp = plotKaryotype(genome = genome, plot.type = 4, ideogram.plotter = NULL, labels.plotter = NULL, plot.params = pp, main = sample_name)
kp = kpAxis(kp, ymin = 0, ymax = max_cn, tick.pos = 0:max_cn)
kp = kpAddCytobandsAsLine(kp)
kp = kpAddChromosomeNames(kp, srt = 45)
kp = kpPoints(kp, data = ratio_norm, y = ratio_norm$Ratio,
  cex = 0.3, ymin = 0, ymax = max_cn, col = alpha("darkolivegreen2", 0.3))
if (length(ratio_amp) > 0) {
  kp = kpPoints(kp, data = ratio_amp, y = ratio_amp$Ratio,
    cex = 0.3, ymin = 0, ymax = max_cn, col = alpha("firebrick2", 0.3))
}
if (length(ratio_del) > 0) {
  kp = kpPoints(kp, data = ratio_del, y = ratio_del$Ratio,
    cex = 0.3, ymin = 0, ymax = max_cn, col = alpha("royalblue4", 0.3))
}
kp = kpPoints(kp, data = ratio_gr, y = ratio_gr$CopyNumber,
  cex = 0.5, ymin = 0, ymax = max_cn, col = "gray20")
dev.off()



# end
