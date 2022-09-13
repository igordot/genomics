#!/usr/bin/env Rscript

"
Streamlined analysis of high-dimensional flow or mass cytometry data.
This script generates a flowSet object from FCS files.

The inputs are a directory with FCS files and a project directory to store the output.
This is designed to be run twice.
If the project directory does not exist, a metadata sample sheet is generated and should be edited as needed.
If a metadata sample sheet is already present, the specified FCS files are imported.

Usage:
  hdcyto-fcs-flowset.R <fcs_dir> <proj_dir>

" -> doc

# output width
options(width = 120)
# print warnings as they occur
options(warn = 1)

# retrieve the command-line arguments
library(docopt)
opts <- docopt(doc)

# define input/output
fcs_dir <- opts$fcs_dir
if (!dir.exists(fcs_dir)) {
  stop("fcs dir does not exist: ", fcs_dir)
}
proj_dir <- opts$proj_dir
samples_csv <- paste0(proj_dir, "/fcs-annot.csv")
if (dir.exists(proj_dir) && !file.exists(samples_csv)) {
  stop("no sample table in dir: ", samples_csv)
}
data_dir <- paste0(proj_dir, "/r-data")
fs_rds <- paste0(data_dir, "/flowset.rds")
if (file.exists(fs_rds)) {
  stop("flowSet already generated: ", fs_rds)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
})

# check if the analysis directory already exists
if (!dir.exists(proj_dir)) {

  # generate a sample sheet for a new analysis

  message("generating sample sheet in: ", proj_dir)

  # get fcs files
  files <- list.files(path = fcs_dir, pattern = "\\.fcs$", full.names = TRUE, recursive = TRUE)

  # attempt to extract sample names based on file names (should be more generic)
  sample_names <- basename(files)
  sample_names <- str_remove(sample_names, ".fcs")
  # remove well ID
  sample_names <- str_remove(sample_names, "^[A-Z]\\d\\d\\s")
  # remove everything after a space
  sample_names <- str_remove(sample_names, "\\s.*")
  sample_names

  # write the sample table (should be manually edited after creating)
  samples_tbl <-
    tibble(
      filename = basename(files),
      sample_id = sample_names,
      patient_id = sample_names,
      condition = ".",
      full_filename = files
    )
  message("saving sample sheet: ", samples_csv)
  dir.create(proj_dir)
  write_csv(samples_tbl, samples_csv)

  message(glue("{proj_dir} sample sheet generated (should be edited for analysis)"))
} else {

  # run analysis based on the existing sample sheet

  message("importing sample sheet: ", samples_csv)
  samples_tbl <- read_csv(samples_csv, show_col_types = FALSE)
  samples_tbl <- samples_tbl %>% dplyr::arrange(filename)
  samples_tbl
  if (!"filename" %in% names(samples_tbl)) {
    stop("sample table should have 'filename' column")
  }
  if (!"sample_id" %in% names(samples_tbl)) {
    stop("sample table should have 'sample_id' column (for CATALYST)")
  }
  if (!"patient_id" %in% names(samples_tbl)) {
    stop("sample table should have 'patient_id' column (for CATALYST)")
  }
  if (!"condition" %in% names(samples_tbl)) {
    stop("sample table should have 'condition' column (for CATALYST)")
  }
  if (ncol(samples_tbl) < 3) {
    stop("sample table should have sample info")
  }
  if (length(unique(samples_tbl$condition)) < 2) {
    stop("sample table should have multiple conditions")
  }

  suppressPackageStartupMessages(library(flowCore))

  # generate an AnnotatedDataFrame
  samples_adf <- as.data.frame(samples_tbl)
  rownames(samples_adf) <- samples_adf$filename
  samples_adf <- new("AnnotatedDataFrame", data = samples_adf)

  message("importing FCS files: ", fcs_dir)
  # fs <- read.flowSet(files, transformation = FALSE, truncate_max_range = FALSE)
  # fs <- read.flowSet(files, transformation = "scale", alter.names = TRUE)
  # fs <- read.flowSet(files, alter.names = TRUE)
  fs <- suppressWarnings(read.flowSet(path = fcs_dir, alter.names = TRUE, transformation = FALSE, phenoData = samples_adf))

  # str(fs)

  pData(fs)

  # save a table of parameters/antibodies
  # colnames(fs)
  fs[[1]]@parameters@data %>% write_csv(glue("{proj_dir}/fcs-parameters.csv"))

  # save a table of the number of cells in each file
  fsApply(fs, nrow) %>%
    as.data.frame() %>%
    as_tibble(rownames = "filename") %>%
    rename(num_cells = V1) %>%
    arrange(filename) %>%
    write_csv(glue("{proj_dir}/fcs-num-cells.csv"))

  message("subsampling flowSet")

  # randomly subsample to avoid extremely large objects
  set.seed(99)
  # sf <- sampleFilter(size = 100000)
  sf <- sampleFilter(size = 500000)
  fres <- filter(fs, sf)
  # summary(fres)
  fs <- Subset(fs, fres)

  # save a table of the number of cells in each file
  fsApply(fs, nrow) %>%
    as.data.frame() %>%
    as_tibble(rownames = "filename") %>%
    rename(num_cells = V1) %>%
    arrange(filename) %>%
    write_csv(glue("{proj_dir}/flowset-num-cells.csv"))

  # fsApply(fs, nrow)
  # dim(exprs(fs[[1]]))
  # exprs(fs[[1]])[1:6, 1:5]

  # save flowSet object
  dir.create(data_dir, showWarnings = FALSE)
  message("saving flowSet: ", fs_rds)
  saveRDS(fs, fs_rds)

  message("flowSet object generated in: ", proj_dir)
}
