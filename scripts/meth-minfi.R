##
## minfi wrapper functions to help streamline the analysis
##


# output width
options(width = 120)
# print warnings as they occur
options(warn = 1)

# dependencies
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readr))

# color scale for plots
plot_colors = c(brewer.pal(9, "Set1"), brewer.pal(8, "Accent"), brewer.pal(8, "Dark2"))

# import data and generate some qc plots
load_data = function(sample_sheet) {

  # check if sample sheet exists
  if (!file.exists(sample_sheet)) stop("sample sheet ", sample_sheet, " does not exist")

  # import sample sheet
  targets = read.csv(sample_sheet, strip.white = TRUE, stringsAsFactors = FALSE)
  targets[targets == ""] = "."

  # sample sheet needs to have "Basename" and "Sentrix_ID" for minfi
  if (!("Basename" %in% colnames(targets))) stop("sample sheet must contain \"Basename\" column")
  if (!("Sentrix_ID" %in% colnames(targets))) stop("sample sheet must contain \"Sentrix_ID\" column")

  # sample sheet needs to have "Sample" and "Condition" for this workflow
  if (!("Sample" %in% colnames(targets))) stop("sample sheet must contain \"Sample\" column")
  if (!("Condition" %in% colnames(targets))) stop("sample sheet must contain \"Condition\" column")

  # add array ID based on basename (to check for batch effects, for example)
  targets$Array = gsub(".*/([0-9]*)_R[0-9][0-9]C[0-9][0-9]", "\\1", targets$Basename)

  message("\n\n ========== minfi::read.metharray.exp() ========== \n\n")

  # red and green channel measurements of the samples (combine() combines two sets of samples)
  raw_set = read.metharray.exp(targets = targets, recursive = TRUE, verbose = TRUE)

  # check that sample names and pData are in the same order (probably not necessary)
  if (!(identical(sampleNames(raw_set), sub(".*/", "", pData(raw_set)$Basename)))) stop("sample names not identical")

  # change sample identifier from "Basename" to "Sample"
  sampleNames(raw_set) = pData(raw_set)$Sample

  # show which array type and corresponding package are being used
  message("array: ", annotation(raw_set)[["array"]])
  message("annotation: ", annotation(raw_set)[["annotation"]])

  # show conditions
  message("samples per condition: ")
  raw_set$Condition %>% table(useNA = "always") %>% print

  message("\n\n ========== minfi::read.qcReport() ========== \n\n")

  # PDF QC report of the most common plots
  qcReport(raw_set, sampGroups=pData(raw_set)$Condition, pdf="plot.qcreport.pdf")

  png("plot.density.raw.condition.png", width = 8, height = 5, units = "in", res = 300)
    densityPlot(raw_set, sampGroups = pData(raw_set)$Condition, pal = plot_colors)
  dev.off()

  png("plot.density.raw.array.png", width = 8, height = 5, units = "in", res = 300)
    densityPlot(raw_set, sampGroups = pData(raw_set)$Array, pal = plot_colors)
  dev.off()

  message("\n\n ========== minfi::detectionP() ========== \n\n")

  # identify failed positions
  det_p = detectionP(raw_set)

  # save detection stats
  det_p_summary = tibble(
    sample = colnames(det_p),
    detected_positions = colSums(det_p < 0.01),
    failed_positions = colSums(det_p >= 0.01),
    failed_positions_pct = round(colMeans(det_p > 0.01), digits = 3)
    )
  det_p_summary = det_p_summary %>% arrange(-failed_positions)
  det_p_summary = det_p_summary %>% mutate(failed_positions_pct = failed_positions_pct * 100)
  write_csv(det_p_summary, "summary.detection.csv")

  return(raw_set)

}

# normalize (functional normalization) raw data using functional normalization and generate some qc plots
normalize_data = function(raw_channel_set) {

  # FunNorm: preprocessFunnorm() -> GenomicRatioSet
  # Noob: preprocessNoob() -> MethylSet -> mapToGenome() -> GenomicMethylSet - > ratioConvert() -> GenomicRatioSet
  # may not be necessary to convert to GenomicRatioSet, getBeta() works with both

  message("\n\n ========== minfi::preprocessRaw() ========== \n\n")

  mset = preprocessRaw(raw_channel_set)

  # plot and save median intensity QC
  qc = getQC(mset)
  mset = addQC(mset, qc = qc)
  png("plot.medianintensity.png", width = 8, height = 8, units = "in", res = 300)
    plotQC(qc)
  dev.off()

  # worst samples (median intensity < 10.5 is failing by default)
  # qc[qc[,"mMed"] < 10.5 | qc[,"uMed"] < 10.5,]

  message("\n\n ========== minfi::preprocessFunnorm() ========== \n\n")

  # functional normalization (FunNorm) - produces GenomicRatioSet
  norm_set = preprocessFunnorm(raw_set)
  # class(norm_set)
  # norm_set
  write(paste0("total probes: ", nrow(norm_set)), file = "norm.log", append = TRUE)

  # identify failed positions
  det_p = detectionP(raw_set)

  # probes detected in at least 90% of the samples
  # normSet = normSet[rowSums(det_p < 0.01) > ncol(det_p)*0.9,]
  # probes detected in all samples
  norm_set = norm_set[rowSums(det_p < 0.01) == ncol(det_p), ]
  write(paste0("detected probes: ", nrow(norm_set)), file = "norm.log", append = TRUE)

  # drop the probes that contain either a SNP at the CpG interrogation or at the single nucleotide extension
  norm_set = addSnpInfo(norm_set)
  # head(granges(norm_set))
  norm_set = dropLociWithSnps(norm_set, snps = c("SBE","CpG"), maf = 0)
  # norm_set
  write(paste0("non-SNP probes: ", nrow(norm_set)), file = "norm.log", append = TRUE)

  # sex prediction plot
  png("plot.sex.png", width = 8, height = 8, units = "in", res = 300)
    plotSex(getSex(norm_set), id = sampleNames(norm_set))
  dev.off()

  # annotation
  annot = getAnnotation(norm_set)

  # remove sex probes
  sex_probes = annot$Name[annot$chr %in% c("chrX", "chrY")]
  norm_set = norm_set[!(rownames(norm_set) %in% sex_probes), ]
  write(paste0("non-sex probes: ", nrow(norm_set)), file = "norm.log", append = TRUE)

  beta = getBeta(norm_set)

  png("plot.density.norm.fnorm.png", width = 8, height = 5, units = "in", res = 300)
    densityPlot(beta, sampGroups = pData(norm_set)$Condition, pal = plot_colors)
  dev.off()

  # mds plots (must be an 'RGChannelSet', a 'MethylSet'  or matrix)

  png("plot.mds.raw.condition.png", width = 8, height = 8, units = "in", res = 300)
    mdsPlot(raw_set, numPositions = 10000, sampNames = sampleNames(raw_set), sampGroups = pData(raw_set)$Condition,
            legendPos = "topright", legendNCol = 1, pal = plot_colors)
  dev.off()

  png("plot.mds.norm.fnorm.condition.png", width = 8, height = 8, units = "in", res = 300)
    mdsPlot(beta, numPositions = 10000, sampNames = sampleNames(norm_set), sampGroups = pData(norm_set)$Condition,
            legendPos = "topright", legendNCol = 1, pal = plot_colors)
  dev.off()

  png("plot.mds.norm.fnorm.array.png", width = 8, height = 8, units = "in", res = 300)
    mdsPlot(beta, numPositions = 10000, sampNames = sampleNames(norm_set), sampGroups = pData(norm_set)$Array,
            legendPos = "topright", legendNCol = 1, pal = plot_colors)
  dev.off()

  return(norm_set)

}



# end
