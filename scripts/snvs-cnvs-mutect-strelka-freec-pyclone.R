#!/usr/bin/env Rscript


'
Description:
Convert mutations from Mutect and Strelka VCFs and CNVs from Control-FREEC to PyClone input table.

Usage:
snvs-cnvs-tsv-freec-pyclone.R <sample_name> <mutect_vcf> <strelka_vcf> <cnvs_txt> <out_txt>

Arguments:
<sample_name>  sample name in the tumor:normal format, must match the VCF sample names
<mutect_vcf>   Mutect (GATK4) VCF
<strelka_vcf>  Strelka VCF
<cnvs_txt>     Control-FREEC "ratio.txt" file with ratios, copy numbers, and genotypes for each window
<out_txt>      output PyClone TSV

Options:
-h, --help        show this screen
' -> doc


# print warnings as they occur
options(warn = 1)

# retrieve the command-line arguments
suppressPackageStartupMessages(library(docopt))
opts = docopt(doc)

# relevent arguments
sample_name = opts$sample_name
mutect_vcf = opts$mutect_vcf
strelka_vcf = opts$strelka_vcf
cnvs_txt = opts$cnvs_txt
out_txt = opts$out_txt

# check that input files exists
if (!file.exists(mutect_vcf)) stop("file does not exist: ", mutect_vcf)
if (!file.exists(strelka_vcf)) stop("file does not exist: ", strelka_vcf)
if (!file.exists(cnvs_txt)) stop("file does not exist: ", cnvs_txt)

# load libraries
suppressPackageStartupMessages({
  library("tidyverse")
  library("glue")
  library("vcfR")
  library("GenomicRanges")
  library("stringr")
})

# Mutect 2.1 (GATK 4)
parse_mutect21 = function(vcfr_obj, sample_T, sample_N) {

  # confirm mutect version
  if (!any(str_detect(vcfr_obj@meta, "##Mutect Version=2.1"))) {
    stop("version mismatch: expecting Mutect 2.1")
  }

  # confirm sample names
  if (!any(str_detect(vcfr_obj@meta, glue("##tumor_sample={sample_T}")))) {
    stop("sample mismatch")
  }
  if (!any(str_detect(vcfr_obj@meta, glue("##normal_sample={sample_N}")))) {
    stop("sample mismatch")
  }

  # if a read is considered uninformative, it is counted towards the DP, but not the AD
  # an uninformative read is not reported in the AD, it is still used in calculations for genotyping
  # if uninformative reads are the only reads, we report the potential variant allele, but keep the AD values 0
  # AD is the number of reads that more likely than not support an allele
  # if you have 10 reads, each with 0.6 probability of having a certain alt allele, you get an AD of 10, whereas you get essentially 0.6 x 10 = 6 reads for the purpose of AF

  ##INFO DP Approximate read depth; some reads may have been filtered
  ##INFO TLOD Tumor LOD score
  ##FORMAT AD Allelic depths for the ref and alt alleles in the order listed
  ##FORMAT AF Allele fractions of alternate alleles in the tumor
  ##FORMAT DP Approximate read depth (reads with MQ=255 or with bad mates are filtered

  # convert vcfR object to a tibble
  muts_tbl = vcfR2tidy(vcfr_obj, single_frame = TRUE, verbose = FALSE,
                       info_fields = c("DP", "TLOD"),
                       format_fields = c("AD", "AF", "DP"))
  muts_tbl = muts_tbl$dat

  # extract relevant metrics
  muts_tbl =
    muts_tbl %>%
    # unique mutation identifier (for joining T and N)
    mutate(mutation_id = as.character(glue("{CHROM}:{POS}:{REF}:{ALT}"))) %>%
    filter(FILTER == "PASS") %>%
    separate(gt_AD, into = c("ref_counts", "alt_counts"), sep = ",", convert = TRUE, extra = "drop") %>%
    # manual AF calculation for comparison
    mutate(myAF = alt_counts / (ref_counts + alt_counts)) %>%
    mutate(
      QUAL = round(as.numeric(TLOD), 1),
      T_DEPTH = ref_counts + alt_counts,
      T_FREQ = round(as.numeric(gt_AF), 3)
    )

  # extract samples T and N to put side by side ("wide" format)
  snvs_n_tbl =
    muts_tbl %>%
    filter(Indiv == sample_N) %>%
    dplyr::rename(N_DEPTH = T_DEPTH, N_FREQ = T_FREQ) %>%
    dplyr::select(mutation_id, N_DEPTH, N_FREQ)
  muts_tbl =
    muts_tbl %>%
    filter(Indiv == sample_T) %>%
    inner_join(snvs_n_tbl, by = "mutation_id")

  # manual filtering
  muts_tbl %>% filter((alt_counts >= 3) & (T_DEPTH >= 10) & (T_FREQ > 0.01) & (T_FREQ > (N_FREQ * 5)))

}

# Strelka 2
parse_strelka2 = function(vcfr_obj, sample_T, sample_N) {

  # confirm strelka version
  if (!any(str_detect(vcfr_obj@meta, "##source_version=2"))) {
    stop("version mismatch: expecting Strelka 2")
  }

  ##INFO QSS Quality score for any somatic snv
  ##INFO SOMATIC Somatic mutation">
  ##INFO QSI Quality score for the ALT haplotype to be present at a significantly different freq in the T and N
  ##FORMAT AU Number of 'A' alleles used in tiers 1,2
  ##FORMAT CU Number of 'C' alleles used in tiers 1,2
  ##FORMAT GU Number of 'G' alleles used in tiers 1,2
  ##FORMAT TU Number of 'T' alleles used in tiers 1,2
  ##FORMAT TAR Reads strongly supporting alternate allele for tiers 1,2
  ##FORMAT TIR Reads strongly supporting indel allele for tiers 1,2

  # convert vcfR object to a tibble
  muts_tbl = vcfR2tidy(vcfr_obj,
                       info_fields = c("SOMATIC", "QSS", "QSI"),
                       format_fields = c("DP", "AU", "CU", "GU", "TU", "TAR", "TIR"),
                       single_frame = TRUE, verbose = FALSE)
  muts_tbl = muts_tbl$dat
  colnames(muts_tbl)

  # extract relevant metrics
  muts_tbl =
    muts_tbl %>%
    mutate(mutation_id = as.character(glue("{CHROM}:{POS}:{REF}:{ALT}"))) %>%
    filter(FILTER == "PASS") %>%
    # extract tier1 counts for each nucleotide
    separate(gt_AU, into = "A_counts", sep = ",", convert = TRUE, extra = "drop") %>%
    separate(gt_CU, into = "C_counts", sep = ",", convert = TRUE, extra = "drop") %>%
    separate(gt_GU, into = "G_counts", sep = ",", convert = TRUE, extra = "drop") %>%
    separate(gt_TU, into = "T_counts", sep = ",", convert = TRUE, extra = "drop") %>%
    separate(gt_TAR, into = "indel_ref_counts", sep = ",", convert = TRUE, extra = "drop") %>%
    separate(gt_TIR, into = "indel_alt_counts", sep = ",", convert = TRUE, extra = "drop") %>%
    # set ref/alt counts
    mutate(
      ref_counts = case_when(
        is.na(QSS) ~ indel_ref_counts,
        REF == "A" ~ A_counts,
        REF == "C" ~ C_counts,
        REF == "G" ~ G_counts,
        REF == "T" ~ T_counts
      )
    ) %>%
    mutate(
      alt_counts = case_when(
        is.na(QSS) ~ indel_alt_counts,
        ALT == "A" ~ A_counts,
        ALT == "C" ~ C_counts,
        ALT == "G" ~ G_counts,
        ALT == "T" ~ T_counts
      )
    ) %>%
    # extract quality
    mutate(
      QUAL = case_when(
        is.na(QSI) ~ QSS,
        is.na(QSS) ~ QSI
      )
    ) %>%
    mutate(T_DEPTH = gt_DP) %>%
    mutate(T_FREQ = round(alt_counts / (ref_counts + alt_counts), 3))

  # extract samples T and N to put side by side ("wide" format)
  snvs_n_tbl =
    muts_tbl %>%
    filter(Indiv == "NORMAL") %>%
    dplyr::rename(N_DEPTH = T_DEPTH, N_FREQ = T_FREQ) %>%
    dplyr::select(mutation_id, N_DEPTH, N_FREQ)
  muts_tbl =
    muts_tbl %>%
    filter(Indiv == "TUMOR") %>%
    inner_join(snvs_n_tbl, by = "mutation_id")

  # manual filtering
  muts_tbl %>% filter((alt_counts >= 3) & (T_DEPTH >= 10) & (T_FREQ > 0.01) & (T_FREQ > (N_FREQ * 5)))

}

# parse either vcf
parse_vcf = function(sample_name, vcf) {

  # split sample name for somatic variants
  if (str_detect(sample_name, ":")) {
    sample_T = str_split_fixed(sample_name, ":", 2)[1]
    sample_N = str_split_fixed(sample_name, ":", 2)[2]
  }

  # import VCF as a vcfR object
  muts_vcfr = read.vcfR(vcf, verbose = FALSE)

  # check if there are any variants
  if (nrow(muts_vcfr@fix) == 0) stop("no variants in imported VCF")

  # determine variant caller based on VCF contents and parse accordingly
  if (any(str_detect(muts_vcfr@meta, "##source=Mutect2"))) {

    # Mutect 2.1 (GATK 4)
    message("parsing Mutect 2.1 VCF")
    caller_type = "somatic"
    vcf_type = "mutect21"
    vcf_tbl = parse_mutect21(vcfr_obj = muts_vcfr, sample_T = sample_T, sample_N = sample_N)

  } else if (any(str_detect(muts_vcfr@meta, "##source=strelka"))) {

    # Strelka 2
    message("parsing Strelka 2 VCF")
    caller_type = "somatic"
    vcf_type = "strelka2"
    vcf_tbl = parse_strelka2(vcfr_obj = muts_vcfr, sample_T = sample_T, sample_N = sample_N)

  } else {

    stop("unknown variant caller")

  }

  # check if table is empty
  if (nrow(vcf_tbl) == 0) stop("output table is empty after parsing")

  # create and sort by mutation_id
  vcf_tbl =
    vcf_tbl %>%
    mutate(mutation_id = as.character(glue("{CHROM}:{POS}:{REF}:{ALT}"))) %>%
    arrange(mutation_id) %>%
    mutate(
      mut_type = case_when(
        str_length(REF) > str_length(ALT) ~ "del",
        str_length(ALT) > str_length(REF) ~ "ins",
        TRUE ~ "pt"
      )
    ) %>%
    filter(mut_type == "pt")

  # output table columns
  vcf_tbl = vcf_tbl %>% mutate(SAMPLE_T = sample_T, SAMPLE_N = sample_N)
  vcf_tbl

}

# parse mutect VCF
muts_mutect_tbl = parse_vcf(sample_name = sample_name, vcf = mutect_vcf)
muts_mutect_tbl$variant_caller = "mutect"

# parse strelka VCF
muts_strelka_tbl = parse_vcf(sample_name = sample_name, vcf = strelka_vcf)
muts_strelka_tbl$variant_caller = "strelka"

# combine both tables
muts_all = bind_rows(muts_mutect_tbl, muts_strelka_tbl)

# determine number of callers (including duplicates)
snvs_tbl =
  muts_all %>%
  add_count(mutation_id, SAMPLE_T, SAMPLE_N) %>%
  mutate(variant_caller = if_else(n == 2, "mutect+strelka", variant_caller)) %>%
  dplyr::select(-QUAL, -n)

# keep one row if called by both callers
# for mutations with entry for each caller, select entry with higher T_FREQ
snvs_tbl =
  snvs_tbl %>%
  group_by(mutation_id, SAMPLE_T, SAMPLE_N) %>%
  arrange(-T_FREQ, -T_DEPTH) %>%
  dplyr::slice(1) %>%
  ungroup()

# chrs to filter out
chr_filter = c("chrX", "chrY", "chrM", "X", "Y", "MT", "M")

# clean up and add PyClone columns
snvs_tbl = snvs_tbl %>%
  filter(!CHROM %in% chr_filter) %>%
  dplyr::rename(chr = CHROM, variant_freq = T_FREQ, var_counts = alt_counts) %>%
  mutate(start = POS, end = POS)

# convert SNVs to GRanges for overlapping
snvs_gr = makeGRangesFromDataFrame(df = snvs_tbl,
                                   ignore.strand = TRUE,
                                   keep.extra.columns = TRUE,
                                   starts.in.df.are.0based = FALSE)
names(snvs_gr) = snvs_gr$mutation_id

# import CNVs file
cnvs_tbl = read_tsv(cnvs_txt, guess_max = 500000, progress = FALSE)

# clean up
cnvs_tbl = cnvs_tbl %>%
  transmute(chr = Chromosome, start = Start - 1, end = Start, genotype = Genotype) %>%
  dplyr::select(chr, start, end, genotype) %>%
  mutate(cnv_id = str_c(chr, ":", end)) %>%
  filter(!chr %in% chr_filter)

# calculate window size (most common start site difference)
win_size = cnvs_tbl %>%
  mutate(win_size = (start - lag(start))) %>%
  filter(win_size > 100) %>%
  na.omit() %>%
  pull(win_size)
win_size = unique(win_size)[which.max(tabulate(match(win_size, unique(win_size))))]

# update coordinates to reflect window size (0-based)
# can be extracted from targeted sequencing results, but not WGS
cnvs_tbl = cnvs_tbl %>% mutate(end = start + win_size)

# add PyClone columns
cnvs_tbl = cnvs_tbl %>%
  filter(genotype != "-", genotype != "") %>%
  mutate(normal_cn = 2,
         minor_cn = str_count(genotype, "B"),
         major_cn = str_count(genotype, "A"))

# convert CNVs to GRanges for overlapping
cnvs_gr = makeGRangesFromDataFrame(df = cnvs_tbl,
                                   ignore.strand = TRUE,
                                   keep.extra.columns = TRUE,
                                   starts.in.df.are.0based = TRUE)
names(cnvs_gr) = cnvs_gr$cnv_id

# set seqlevels in both GRanges to UCSC style
seqlevelsStyle(snvs_gr) = "UCSC"
seqlevelsStyle(cnvs_gr) = "UCSC"

# get overlaps for SNVs and CNVs (consider adjacent windows as well)
overlaps = distanceToNearest(x = snvs_gr, subject = cnvs_gr)
overlaps = overlaps %>% as("data.frame") %>% as_tibble() %>% filter(distance < win_size * 2)

# convert overlaps from identifiers to names
overlaps$mutation_id = snvs_gr[overlaps$queryHits] %>% names()
overlaps$cnv_id = cnvs_gr[overlaps$subjectHits] %>% names()
overlaps = overlaps %>% dplyr::select(-queryHits, -subjectHits)

# merge overlapping SNVs and CNVs
overlaps = overlaps %>%
  inner_join(snvs_tbl, by = "mutation_id") %>%
  inner_join(cnvs_tbl, by = "cnv_id")

# extract columns for PyClone
pyclone_tbl = overlaps %>%
  dplyr::select(mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn, variant_freq, genotype)

# save PyClone table
write_tsv(pyclone_tbl, path = out_txt)



# end
