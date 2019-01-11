#!/usr/bin/env Rscript

"
Analysis of 10x Genomics Chromium single cell RNA-seq data using Seurat (version 3.0) starting with Cell Ranger output.

Basic workflow steps:
  1 - create - import counts matrix, perform initial QC, and calculate various variance metrics (slowest step)
  2 - cluster - perform clustering based on number of PCs
  3 - identify - identify clusters based on specified clustering/resolution (higher resolution for more clusters)

Optional steps:
  combine - merge multiple samples/libraries
  integrate - perform integration (batch correction) across multiple sample batches
  align - perform alignment based on number of CCs (after performing CCA)
  de - differential expression between samples/libraries within clusters

Usage:
  scrna-10x-seurat-3.R create <analysis_dir> <sample_name> <sample_dir> [--min_genes=<n> --max_genes=<n> --mt=<n>]
  scrna-10x-seurat-3.R cluster <analysis_dir> <num_pcs>
  scrna-10x-seurat-3.R identify <analysis_dir> <resolution>
  scrna-10x-seurat-3.R combine <analysis_dir> <sample_analysis_dir>...
  scrna-10x-seurat-3.R integrate <analysis_dir> <num_dim> <batch_analysis_dir>...
  scrna-10x-seurat-3.R de <analysis_dir> <resolution>
  scrna-10x-seurat-3.R --help

Options:
  --min_genes=<n>   cutoff for minimum number of genes per cell (2nd percentile if not specified)
  --max_genes=<n>   cutoff for maximum number of genes per cell (98th percentile if not specified)
  --mt=<n>          cutoff for mitochondrial genes percentage per cell [default: 10]
  -h, --help        show this screen
" -> doc


# ========== functions ==========


# load dependencies
load_libraries = function() {

  message("\n\n ========== load libraries ========== \n\n")

  suppressPackageStartupMessages({
    library(magrittr)
    library(glue)
    library(Seurat)
    library(future)
    library(Matrix)
    library(tidyverse)
    library(cowplot)
    library(scales)
    library(RColorBrewer)
    library(ggsci)
    library(eulerr)
    library(UpSetR)
  })

}

# create a single Seurat object from multiple 10x Cell Ranger outputs
# takes vector of one or more sample_names and sample_dirs
# can work with multiple samples, but the appropriate way is to use "combine" with objects that are pre-filtered
load_sample_counts_matrix = function(sample_names, sample_dirs) {

  message("\n\n ========== import cell ranger counts matrix ========== \n\n")

  merged_counts_matrix = NULL

  for (i in 1:length(sample_names)) {

    sample_name = sample_names[i]
    sample_dir = sample_dirs[i]

    message("loading counts matrix for sample: ", sample_name)

    # check if sample dir is valid
    if (!dir.exists(sample_dir)) stop(glue("dir {sample_dir} does not exist"))

    # determine counts matrix directory (HDF5 is not the preferred option)
    # "filtered_gene_bc_matrices" for single library
    # "filtered_gene_bc_matrices_mex" for aggregated
    # Cell Ranger 3.0: "genes" has been replaced by "features" to account for feature barcoding
    # Cell Ranger 3.0: the matrix and barcode files are now gzipped
    data_dir = glue("{sample_dir}/outs")
    if (!dir.exists(data_dir)) stop(glue("dir {sample_dir} does not contain outs directory"))
    data_dir = list.files(path = data_dir, pattern = "matrix.mtx", full.names = TRUE, recursive = TRUE)
    data_dir = str_subset(data_dir, "filtered_.*_bc_matri")[1]
    data_dir = dirname(data_dir)
    if (!dir.exists(data_dir)) stop(glue("dir {sample_dir} does not contain matrix.mtx"))

    message("loading counts matrix dir: ", data_dir)

    counts_matrix = Read10X(data_dir)

    message(glue("library {sample_name} cells: {ncol(counts_matrix)}"))
    message(glue("library {sample_name} genes: {nrow(counts_matrix)}"))
    message(" ")

    # log to file
    write(glue("library {sample_name} cells: {ncol(counts_matrix)}"), file = "create.log", append = TRUE)
    write(glue("library {sample_name} genes: {nrow(counts_matrix)}"), file = "create.log", append = TRUE)

    # clean up counts matrix to make it more readable
    counts_matrix = counts_matrix[sort(rownames(counts_matrix)), ]
    colnames(counts_matrix) = str_c(sample_name, ":", colnames(counts_matrix))

    # combine current matrix with previous
    if (i == 1) {

      # skip if there is no previous matrix
      merged_counts_matrix = counts_matrix

    } else {

      # check if genes are the same for current and previous matrices
      if (!identical(rownames(merged_counts_matrix), rownames(counts_matrix))) {

        # generate a warning, since this is probably a mistake
        warning("counts matrix genes are not the same for different libraries")
        Sys.sleep(1)

        # get common genes
        common_genes = intersect(rownames(merged_counts_matrix), rownames(counts_matrix))
        common_genes = sort(common_genes)
        message("num genes for previous libraries: ", length(rownames(merged_counts_matrix)))
        message("num genes for current library:   ", length(rownames(counts_matrix)))
        message("num genes in common:        ", length(common_genes))

        # exit if the number of overlapping genes is too few
        if (length(common_genes) < (length(rownames(counts_matrix)) * 0.9)) stop("libraries have too few genes in common")

        # subset current and previous matrix to overlapping genes
        merged_counts_matrix = merged_counts_matrix[common_genes, ]
        counts_matrix = counts_matrix[common_genes, ]

      }

      # combine current matrix with previous
      merged_counts_matrix = cbind(merged_counts_matrix, counts_matrix)
      Sys.sleep(1)

    }

  }

  # create a Seurat object
  s_obj = create_seurat_obj(counts_matrix = merged_counts_matrix)

  return(s_obj)

}

# convert a sparse matrix of counts to a Seurat object and generate some QC plots
create_seurat_obj = function(counts_matrix, proj_name = NULL, sample_dir = NULL) {

  message("\n\n ========== save raw counts matrix ========== \n\n")

  # log to file
  write(glue("input cells: {ncol(counts_matrix)}"), file = "create.log", append = TRUE)
  write(glue("input genes: {nrow(counts_matrix)}"), file = "create.log", append = TRUE)

  # remove cells with few genes (there is an additional filter in CreateSeuratObject)
  counts_matrix = counts_matrix[, Matrix::colSums(counts_matrix) >= 250]

  # remove genes without counts (there is an additional filter in CreateSeuratObject)
  counts_matrix = counts_matrix[Matrix::rowSums(counts_matrix) >= 3, ]

  # log to file
  write(glue("detectable cells: {ncol(counts_matrix)}"), file = "create.log", append = TRUE)
  write(glue("detectable genes: {nrow(counts_matrix)}"), file = "create.log", append = TRUE)

  # save counts matrix as a standard text file and gzip
  # counts_matrix_filename = "counts.raw.txt"
  # write.table(as.matrix(counts_matrix), file = counts_matrix_filename, quote = FALSE, sep = "\t", col.names = NA)
  # system(paste0("gzip ", counts_matrix_filename))

  # save counts matrix as a csv file (to be consistent with the rest of the tables)
  raw_data = counts_matrix %>% as.matrix() %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
  write_csv(raw_data, path = "counts.raw.csv.gz")
  rm(raw_data)

  message("\n\n ========== create seurat object ========== \n\n")

  if (is.null(proj_name)) {

    # if name is not set, then it's a manually merged counts matrix
    s_obj = CreateSeuratObject(counts = counts_matrix, min.cells = 5, min.features = 250, project = "proj",
                               names.field = 1, names.delim = ":")
    rm(counts_matrix)

  } else if (proj_name == "aggregated") {

    # multiple libraries combined using Cell Ranger (cellranger aggr)

    # setup taking into consideration aggregated names delimiter
    s_obj = CreateSeuratObject(counts = counts_matrix, min.cells = 5, min.features = 250, project = proj_name,
                               names.field = 2, names.delim = "-")

    # import cellranger aggr sample sheet
    sample_sheet_csv = paste0(sample_dir, "/outs/aggregation_csv.csv")
    sample_sheet = read.csv(sample_sheet_csv, stringsAsFactors = FALSE)
    message("samples: ", paste(sample_sheet[, 1], collapse=", "))

    # change s_obj@meta.data$orig.ident sample identities from numbers to names
    s_obj[["orig.ident"]][, 1] = factor(sample_sheet[s_obj[["orig.ident"]][, 1], 1])
    # set s_obj@ident to the new s_obj@meta.data$orig.ident
    s_obj = set_identity(seurat_obj = s_obj, group_var = "orig.ident")

  } else {

    stop("project name set to unknown value")

  }

  message(glue("imported cells: {ncol(s_obj)}"))
  message(glue("imported genes: {nrow(s_obj)}"))
  message(" ")

  # log to file
  write(glue("imported cells: {ncol(s_obj)}"), file = "create.log", append = TRUE)
  write(glue("imported genes: {nrow(s_obj)}"), file = "create.log", append = TRUE)

  # rename nCount_RNA and nFeature_RNA slots to make them more clear
  s_obj$num_UMIs = s_obj$nCount_RNA
  s_obj$num_genes = s_obj$nFeature_RNA

  # nFeature_RNA and nCount_RNA are automatically calculated for every object by Seurat
  # calculate the percentage of mitochondrial genes here and store it in percent.mito using the AddMetaData
  mt_genes = grep("^MT-", rownames(GetAssayData(s_obj)), ignore.case = TRUE, value = TRUE)
  percent_mt = Matrix::colSums(GetAssayData(s_obj)[mt_genes, ]) / Matrix::colSums(GetAssayData(s_obj))
  percent_mt = round(percent_mt * 100, digits = 3)

  # add columns to object@meta.data, and is a great place to stash QC stats
  s_obj = AddMetaData(s_obj, metadata = percent_mt, col.name = "pct_mito")

  message("\n\n ========== nFeature_RNA/nCount_RNA/percent_mito plots ========== \n\n")

  # create a named color scheme to ensure names and colors are in the proper order
  sample_names = s_obj$orig.ident %>% as.character() %>% sort() %>% unique()
  colors_samples_named = colors_samples[1:length(sample_names)]
  names(colors_samples_named) = sample_names

  vln_theme =
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "none"
    )
  suppressMessages({
    dist_unfilt_nft_plot =
      VlnPlot(
        s_obj, features = "num_genes", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_unfilt_nct_plot =
      VlnPlot(
        s_obj, features = "num_UMIs", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_unfilt_pmt_plot =
      VlnPlot(
        s_obj, features = "pct_mito", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_unfilt_plot = plot_grid(dist_unfilt_nft_plot, dist_unfilt_nct_plot, dist_unfilt_pmt_plot, ncol = 3)
    ggsave("qc.distribution.unfiltered.png", plot = dist_unfilt_plot, width = 10, height = 6, units = "in")
  })
  Sys.sleep(1)

  cor_ncr_nfr_plot =
    FeatureScatter(
      s_obj, feature1 = "num_UMIs", feature2 = "num_genes", group.by = "orig.ident", cols = colors_samples
    ) +
    theme(aspect.ratio = 1)
  cor_ncr_pmt_plot =
    FeatureScatter(
      s_obj, feature1 = "num_UMIs", feature2 = "pct_mito", group.by = "orig.ident", cols = colors_samples
    ) +
    theme(aspect.ratio = 1)
  cor_nfr_pmt_plot =
    FeatureScatter(
      s_obj, feature1 = "num_genes", feature2 = "pct_mito", group.by = "orig.ident", cols = colors_samples
    ) +
    theme(aspect.ratio = 1)
  cor_unfilt_plot = plot_grid(cor_ncr_nfr_plot, cor_ncr_pmt_plot, cor_nfr_pmt_plot, ncol = 3)
  ggsave("qc.correlations.unfiltered.png", plot = cor_unfilt_plot, width = 18, height = 5, units = "in")
  Sys.sleep(1)

  # check distribution of gene counts and mitochondrial percentage
  low_quantiles = c(0.05, 0.02, 0.01, 0.001)
  high_quantiles = c(0.95, 0.98, 0.99, 0.999)
  message("num genes low percentiles:")
  s_obj$num_genes %>% quantile(low_quantiles) %>% round(1) %>% print()
  message(" ")
  message("num genes high percentiles:")
  s_obj$num_genes %>% quantile(high_quantiles) %>% round(1) %>% print()
  message(" ")
  message("pct mito high percentiles:")
  s_obj$pct_mito %>% quantile(high_quantiles) %>% round(1) %>% print()
  message(" ")

  # save unfiltered cell metadata
  s_obj@meta.data %>%
    rownames_to_column("cell") %>% as_tibble() %>%
    mutate(sample_name = orig.ident) %>%
    write_excel_csv(path = "metadata.unfiltered.csv")

  return(s_obj)

}

# filter data by number of genes and mitochondrial percentage
filter_data = function(seurat_obj, min_genes = NULL, max_genes = NULL, max_mt = 10) {

  s_obj = seurat_obj

  message("\n\n ========== filter data matrix ========== \n\n")

  # log the unfiltered gene numbers to file
  write(glue("unfiltered min genes: {min(s_obj$num_genes)}"), file = "create.log", append = TRUE)
  write(glue("unfiltered max genes: {max(s_obj$num_genes)}"), file = "create.log", append = TRUE)
  write(glue("unfiltered mean num genes: {round(mean(s_obj$num_genes), 3)}"), file = "create.log", append = TRUE)
  write(glue("unfiltered median num genes: {median(s_obj$num_genes)}"), file = "create.log", append = TRUE)

  # convert arguments to integers (command line arguments end up as characters)
  min_genes = as.numeric(min_genes)
  max_genes = as.numeric(max_genes)
  max_mt = as.numeric(max_mt)

  # default cutoffs (gene numbers rounded to nearest 10)
  # as.numeric() converts NULLs to 0 length numerics, so can't use is.null()
  if (!length(min_genes)) min_genes = s_obj$num_genes %>% quantile(0.02, names = FALSE) %>% round(-1)
  if (!length(max_genes)) max_genes = s_obj$num_genes %>% quantile(0.98, names = FALSE) %>% round(-1)
  if (!length(max_mt)) max_mt = 10

  message(glue("min genes cutoff: {min_genes}"))
  message(glue("max genes cutoff: {max_genes}"))
  message(glue("max mitochondrial percentage cutoff: {max_mt}"))
  message(" ")

  # log the cutoffs to file
  write(glue("min genes cutoff: {min_genes}"), file = "create.log", append = TRUE)
  write(glue("max genes cutoff: {max_genes}"), file = "create.log", append = TRUE)
  write(glue("max mitochondrial percentage cutoff: {max_mt}"), file = "create.log", append = TRUE)

  message(glue("imported cells: {ncol(s_obj)}"))
  message(glue("imported genes: {nrow(s_obj)}"))

  # filter
  cells_subset =
    seurat_obj@meta.data %>%
    rownames_to_column("cell") %>%
    filter(nFeature_RNA > min_genes & nFeature_RNA < max_genes & pct_mito < max_mt) %>%
    pull(cell)
  s_obj = subset(s_obj, cells = cells_subset)

  message("filtered cells: ", ncol(s_obj))
  message("filtered genes: ", nrow(s_obj))

  # create a named color scheme to ensure names and colors are in the proper order
  sample_names = s_obj$orig.ident %>% as.character() %>% sort() %>% unique()
  colors_samples_named = colors_samples[1:length(sample_names)]
  names(colors_samples_named) = sample_names

  vln_theme =
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "none"
    )
  suppressMessages({
    dist_filt_nft_plot =
      VlnPlot(
        s_obj, features = "num_genes", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_filt_nct_plot =
      VlnPlot(
        s_obj, features = "num_UMIs", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_filt_pmt_plot =
      VlnPlot(
        s_obj, features = "pct_mito", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_filt_plot = plot_grid(dist_filt_nft_plot, dist_filt_nct_plot, dist_filt_pmt_plot, ncol = 3)
    ggsave("qc.distribution.filtered.png", plot = dist_filt_plot, width = 10, height = 6, units = "in")
  })
  Sys.sleep(1)

  # after removing unwanted cells from the dataset, normalize the data
  # LogNormalize:
  # - normalizes the gene expression measurements for each cell by the total expression
  # - multiplies this by a scale factor (10,000 by default)
  # - log-transforms the result
  s_obj = NormalizeData(s_obj, normalization.method = "LogNormalize", scale.factor = 100000, verbose = FALSE)

  # save counts matrix as a basic gzipped text file
  # object@data stores normalized and log-transformed single cell expression
  # used for visualizations, such as violin and feature plots, most diff exp tests, finding high-variance genes
  counts_norm = GetAssayData(s_obj) %>% as.matrix() %>% round(3)
  counts_norm = counts_norm %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
  write_csv(counts_norm, path = "counts.normalized.csv.gz")

  # log to file
  write(glue("filtered cells: {ncol(s_obj)}"), file = "create.log", append = TRUE)
  write(glue("filtered genes: {nrow(s_obj)}"), file = "create.log", append = TRUE)
  write(glue("filtered mean num genes: {round(mean(s_obj$num_genes), 3)}"), file = "create.log", append = TRUE)
  write(glue("filtered median num genes: {median(s_obj$num_genes)}"), file = "create.log", append = TRUE)

  return(s_obj)

}

# merge multiple Seurat objects
combine_seurat_obj = function(original_wd, sample_analysis_dirs) {

  if (length(sample_analysis_dirs) < 2) stop("must have at least 2 samples to merge")

  message("\n\n ========== combine samples ========== \n\n")

  seurat_obj_list = list()
  for (i in 1:length(sample_analysis_dirs)) {

    sample_analysis_dir = sample_analysis_dirs[i]
    sample_analysis_dir = glue("{original_wd}/{sample_analysis_dir}")
    sample_seurat_rds = glue("{sample_analysis_dir}/seurat_obj.rds")

    # check if analysis dir is valid
    if (!dir.exists(sample_analysis_dir)) stop(glue("dir {sample_analysis_dir} does not exist"))
    # check if seurat object exists
    if (!file.exists(sample_seurat_rds)) stop(glue("seurat object rds {sample_seurat_rds} does not exist"))

    # load seurat object
    seurat_obj_list[[i]] = readRDS(sample_seurat_rds)

    # clean up object
    seurat_obj_list[[i]]@assays$RNA@var.features = vector()
    seurat_obj_list[[i]]@assays$RNA@scale.data = matrix()
    seurat_obj_list[[i]]@reductions = list()
    seurat_obj_list[[i]]@meta.data = seurat_obj_list[[i]]@meta.data %>% select(-starts_with("snn_res"))

    # print single sample sample stats
    # sample_name = seurat_obj_list[[i]]@meta.data[1, "orig.ident"] %>% as.character()
    sample_name = seurat_obj_list[[i]]$orig.ident[1] %>% as.character()
    message(glue("sample {sample_name} dir: {basename(sample_analysis_dir)}"))
    write(glue("sample {sample_name} dir: {basename(sample_analysis_dir)}"), file = "create.log", append = TRUE)
    message(glue("sample {sample_name} cells: {ncol(seurat_obj_list[[i]])}"))
    write(glue("sample {sample_name} cells: {ncol(seurat_obj_list[[i]])}"), file = "create.log", append = TRUE)
    message(glue("sample {sample_name} genes: {nrow(seurat_obj_list[[i]])}"))
    write(glue("sample {sample_name} genes: {nrow(seurat_obj_list[[i]])}"), file = "create.log", append = TRUE)
    message(" ")

  }

  # merge
  merged_obj = merge(seurat_obj_list[[1]], seurat_obj_list[2:length(seurat_obj_list)])
  rm(seurat_obj_list)

  # print combined sample stats
  message(glue("combined unfiltered cells: {ncol(merged_obj)}"))
  write(glue("combined unfiltered cells: {ncol(merged_obj)}"), file = "create.log", append = TRUE)
  message(glue("combined unfiltered genes: {nrow(merged_obj)}"))
  write(glue("combined unfiltered genes: {nrow(merged_obj)}"), file = "create.log", append = TRUE)

  # filter poorly expressed genes (detected in less than 10 cells)
  filtered_genes = Matrix::rowSums(GetAssayData(merged_obj, slot = "counts") > 0)
  filtered_genes = filtered_genes[filtered_genes >= 10] %>% names() %>% sort()
  merged_obj = subset(merged_obj, features = filtered_genes)

  # print combined sample stats
  message(glue("combined cells: {ncol(merged_obj)}"))
  write(glue("combined cells: {ncol(merged_obj)}"), file = "create.log", append = TRUE)
  message(glue("combined genes: {nrow(merged_obj)}"))
  write(glue("combined genes: {nrow(merged_obj)}"), file = "create.log", append = TRUE)

  # save raw counts matrix
  counts_raw = GetAssayData(merged_obj, slot = "counts") %>% as.matrix()
  counts_raw = counts_raw %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
  write_csv(counts_raw, path = "counts.raw.csv.gz")
  rm(counts_raw)

  # save counts matrix as a basic gzipped text file
  counts_norm = GetAssayData(merged_obj) %>% as.matrix() %>% round(3)
  counts_norm = counts_norm %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
  write_csv(counts_norm, path = "counts.normalized.csv.gz")
  rm(counts_norm)

  # create a named color scheme to ensure names and colors are in the proper order
  sample_names = merged_obj$orig.ident %>% as.character() %>% sort() %>% unique()
  colors_samples_named = colors_samples[1:length(sample_names)]
  names(colors_samples_named) = sample_names

  vln_theme =
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "none"
    )
  suppressMessages({
    dist_nft_plot =
      VlnPlot(
        merged_obj, features = "num_genes", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_nct_plot =
      VlnPlot(
        merged_obj, features = "num_UMIs", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_pmt_plot =
      VlnPlot(
        merged_obj, features = "pct_mito", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_plot = plot_grid(dist_nft_plot, dist_nct_plot, dist_pmt_plot, ncol = 3)
    ggsave("qc.distribution.png", plot = dist_plot, width = 15, height = 6, units = "in")
  })
  Sys.sleep(1)

  return(merged_obj)

}

# integrate multiple Seurat objects
integrate_seurat_obj = function(original_wd, sample_analysis_dirs, num_dim) {

  # check if the inputs seems reasonable
  if (length(sample_analysis_dirs) < 2) stop("must have at least 2 samples to merge")
  num_dim = as.integer(num_dim)
  if (num_dim < 5) stop("too few dims: ", num_dim)
  if (num_dim > 50) stop("too many dims: ", num_dim)

  message("\n\n ========== integrate samples ========== \n\n")

  seurat_obj_list = list()
  var_genes_list = list()
  exp_genes = c()
  for (i in 1:length(sample_analysis_dirs)) {

    sample_analysis_dir = sample_analysis_dirs[i]
    sample_analysis_dir = glue("{original_wd}/{sample_analysis_dir}")
    sample_seurat_rds = glue("{sample_analysis_dir}/seurat_obj.rds")

    # check if analysis dir is valid
    if (!dir.exists(sample_analysis_dir)) stop(glue("dir {sample_analysis_dir} does not exist"))
    # check if seurat object exists
    if (!file.exists(sample_seurat_rds)) stop(glue("seurat object rds {sample_seurat_rds} does not exist"))

    # load seurat object
    seurat_obj_list[[i]] = readRDS(sample_seurat_rds)
    # sample_name = seurat_obj_list[[i]]@meta.data[1, "orig.ident"] %>% as.character()
    sample_name = seurat_obj_list[[i]]$orig.ident[1] %>% as.character()

    # clean up object
    seurat_obj_list[[i]]@assays$RNA@scale.data = matrix()
    seurat_obj_list[[i]]@reductions = list()
    seurat_obj_list[[i]]@meta.data = seurat_obj_list[[i]]@meta.data %>% select(-starts_with("snn_res"))

    # save expressed genes keeping only genes present in all the datasets (for genes to integrate in IntegrateData)
    if (length(exp_genes) > 0) {
      exp_genes = intersect(exp_genes, rownames(seurat_obj_list[[i]])) %>% sort()
    } else {
      exp_genes = rownames(seurat_obj_list[[i]])
    }

    # save variable genes
    var_genes_list[[sample_name]] = VariableFeatures(seurat_obj_list[[i]])

    # print single sample sample stats
    message(glue("sample {sample_name} dir: {basename(sample_analysis_dir)}"))
    write(glue("sample {sample_name} dir: {basename(sample_analysis_dir)}"), file = "create.log", append = TRUE)
    message(glue("sample {sample_name} cells: {ncol(seurat_obj_list[[i]])}"))
    write(glue("sample {sample_name} cells: {ncol(seurat_obj_list[[i]])}"), file = "create.log", append = TRUE)
    message(glue("sample {sample_name} genes: {nrow(seurat_obj_list[[i]])}"))
    write(glue("sample {sample_name} genes: {nrow(seurat_obj_list[[i]])}"), file = "create.log", append = TRUE)
    message(" ")

  }

  # euler plot of variable gene overlaps (becomes unreadable and can take days for many overlaps)
  if (length(var_genes_list) < 8) {
    colors_euler = colors_samples[1:length(var_genes_list)]
    euler_fit = euler(var_genes_list, shape = "ellipse")
    euler_plot = plot(euler_fit,
                      fills = list(fill = colors_euler, alpha = 0.7),
                      edges = list(col = colors_euler))
    png("variance.vargenes.euler.png", res = 200, width = 5, height = 5, units = "in")
      print(euler_plot)
    dev.off()
  }

  # upset plot of variable gene overlaps
  png("variance.vargenes.upset.png", res = 200, width = 8, height = 5, units = "in")
    upset(fromList(var_genes_list), nsets = 50, nintersects = 15, order.by = "freq", mb.ratio = c(0.5, 0.5))
  dev.off()

  message("\n\n ========== Seurat::FindIntegrationAnchors() ========== \n\n")

  # find the integration anchors
  anchors = FindIntegrationAnchors(object.list = seurat_obj_list, anchor.features = 2000, dims = 1:num_dim)
  rm(seurat_obj_list)

  message("\n\n ========== Seurat::IntegrateData() ========== \n\n")

  # integrating all genes may cause issues and may not add any relevant information
  # integrated_obj = IntegrateData(anchorset = anchors, dims = 1:num_dim, features.to.integrate = exp_genes)
  integrated_obj = IntegrateData(anchorset = anchors, dims = 1:num_dim)
  rm(anchors)

  # after running IntegrateData, the Seurat object will contain a new Assay with the integrated expression matrix
  # the original (uncorrected values) are still stored in the object in the “RNA” assay

  # switch to integrated assay
  DefaultAssay(integrated_obj) = "integrated"

  # print integrated sample stats
  message(glue("integrated unfiltered cells: {ncol(integrated_obj)}"))
  write(glue("integrated unfiltered cells: {ncol(integrated_obj)}"), file = "create.log", append = TRUE)
  message(glue("integrated unfiltered genes: {nrow(integrated_obj)}"))
  write(glue("integrated unfiltered genes: {nrow(integrated_obj)}"), file = "create.log", append = TRUE)

  # filter poorly expressed genes (detected in less than 10 cells)
  filtered_genes = Matrix::rowSums(GetAssayData(integrated_obj, assay = "RNA", slot = "counts") > 0)
  filtered_genes = filtered_genes[filtered_genes >= 10] %>% names() %>% sort()
  integrated_obj = subset(integrated_obj, features = filtered_genes)

  # print integrated sample stats
  message(glue("integrated cells: {ncol(integrated_obj)}"))
  write(glue("integrated cells: {ncol(integrated_obj)}"), file = "create.log", append = TRUE)
  message(glue("integrated genes: {nrow(integrated_obj)}"))
  write(glue("integrated genes: {nrow(integrated_obj)}"), file = "create.log", append = TRUE)

  # save raw counts matrix
  counts_raw = GetAssayData(integrated_obj, assay = "RNA", slot = "counts") %>% as.matrix()
  counts_raw = counts_raw %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
  write_csv(counts_raw, path = "counts.raw.csv.gz")
  rm(counts_raw)

  # save normalized counts matrix
  counts_norm = GetAssayData(integrated_obj, assay = "RNA") %>% as.matrix() %>% round(3)
  counts_norm = counts_norm %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
  write_csv(counts_norm, path = "counts.normalized.csv.gz")
  rm(counts_norm)

  # save integrated counts matrix
  # counts_int = GetAssayData(integrated_obj, assay = "integrated") %>% as.matrix() %>% round(3)
  # counts_int = counts_int %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
  # write_csv(counts_int, path = "counts.integrated.csv.gz")
  # rm(counts_int)

  vln_theme =
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "none"
    )
  suppressMessages({
    dist_nft_plot =
      VlnPlot(
        integrated_obj, features = "num_genes", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_nct_plot =
      VlnPlot(
        integrated_obj, features = "num_UMIs", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_pmt_plot =
      VlnPlot(
        integrated_obj, features = "pct_mito", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_plot = plot_grid(dist_nft_plot, dist_nct_plot, dist_pmt_plot, ncol = 3)
    ggsave("qc.distribution.png", plot = dist_plot, width = 15, height = 6, units = "in")
  })
  Sys.sleep(1)

  return(integrated_obj)

}

# calculate various variance metrics and perform basic analysis
# PC selection approaches:
# - PCHeatmap - more supervised, exploring PCs to determine relevant sources of heterogeneity
# - PCElbowPlot - heuristic that is commonly used and can be calculated instantly
# - JackStrawPlot - implements a statistical test based on a random null model, but is time-consuming
# jackStraw procedure is very slow, so skip for large projects (>10,000 cells)
calculate_variance = function(seurat_obj, jackstraw_max_cells = 10000) {

  s_obj = seurat_obj

  message("\n\n ========== Seurat::FindVariableGenes() ========== \n\n")

  # identify features that are outliers on a 'mean variability plot'
  # Seurat v3 implements an improved method based on a variance stabilizing transformation ("vst")
  s_obj = FindVariableFeatures(s_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

  # export highly variable feature information (mean, variance, variance standardized)
  hvf_tbl = HVFInfo(s_obj) %>% round(3) %>% rownames_to_column("gene") %>% arrange(-variance.standardized)
  write_excel_csv(hvf_tbl, path = "variance.csv")

  # plot variance
  var_plot = VariableFeaturePlot(s_obj, pt.size = 0.5)
  var_plot = LabelPoints(var_plot, points = head(hvf_tbl$gene, 30), repel = TRUE, xnudge = 0, ynudge = 0)
  ggsave("variance.features.png", plot = var_plot, width = 12, height = 5, units = "in")

  message("\n\n ========== Seurat::ScaleData() ========== \n\n")

  # regress out unwanted sources of variation
  # regressing uninteresting sources of variation can improve dimensionality reduction and clustering
  # could include technical noise, batch effects, biological sources of variation (cell cycle stage)
  # scaled z-scored residuals of these models are stored in scale.data slot
  # used for dimensionality reduction and clustering
  # RegressOut function has been deprecated, and replaced with the vars.to.regress argument in ScaleData
  # s_obj = ScaleData(s_obj, features = rownames(s_obj), vars.to.regress = c("num_UMIs", "pct_mito"), verbose = FALSE)
  s_obj = ScaleData(s_obj, vars.to.regress = c("num_UMIs", "pct_mito"), verbose = FALSE)

  message("\n\n ========== Seurat::PCA() ========== \n\n")

  # use fewer PCs for small datasets
  num_pcs = 50
  if (ncol(s_obj) < 100) num_pcs = 20
  if (ncol(s_obj) < 25) num_pcs = 5

  # PCA on the scaled data
  # PCA calculation stored in object[["pca"]]
  s_obj = RunPCA(s_obj, assay = "RNA", features = VariableFeatures(s_obj), npcs = num_pcs, verbose = FALSE)

  # plot the output of PCA analysis (shuffle cells so any one group does not appear overrepresented due to ordering)
  pca_plot =
    DimPlot(
      s_obj, cells = sample(colnames(s_obj)), group.by = "orig.ident", reduction = "pca",
      pt.size = 0.5, cols = colors_samples
    ) +
    theme(aspect.ratio = 1)
  ggsave("variance.pca.png", plot = pca_plot, width = 8, height = 6, units = "in")

  message("\n\n ========== Seurat::DimHeatmap() ========== \n\n")

  # PCHeatmap (former) allows for easy exploration of the primary sources of heterogeneity in a dataset
  png("variance.pca.heatmap.png", res = 300, width = 10, height = 16, units = "in")
    DimHeatmap(s_obj, reduction = "pca", dims = 1:15, nfeatures = 20, cells = 250, fast = TRUE)
  dev.off()

  message("\n\n ========== Seurat::PCElbowPlot() ========== \n\n")

  # a more ad hoc method for determining PCs to use, draw cutoff where there is a clear elbow in the graph
  elbow_plot = ElbowPlot(s_obj, reduction = "pca", ndims = num_pcs)
  ggsave("variance.pca.elbow.png", plot = elbow_plot, width = 8, height = 5, units = "in")

  # resampling test inspired by the jackStraw procedure - very slow, so skip for large projects (>10,000 cells)
  if (ncol(s_obj) < jackstraw_max_cells) {

    message("\n\n ========== Seurat::JackStraw() ========== \n\n")

    # determine statistical significance of PCA scores
    s_obj = JackStraw(s_obj, assay = "RNA", reduction = "pca", dims = num_pcs, verbose = FALSE)

    # compute Jackstraw scores significance
    s_obj = ScoreJackStraw(s_obj, reduction = "pca", dims = 1:num_pcs, do.plot = FALSE)

    # plot the results of the JackStraw analysis for PCA significance
    # significant PCs will show a strong enrichment of genes with low p-values (solid curve above the dashed line)
    jackstraw_plot =
      JackStrawPlot(s_obj, reduction = "pca", dims = 1:num_pcs) +
      guides(col = guide_legend(ncol = 2))
    ggsave("variance.pca.jackstraw.png", plot = jackstraw_plot, width = 12, height = 6, units = "in")

  }

  return(s_obj)

}

# calculate various variance metrics and perform basic analysis (integrated analysis workflow)
calculate_variance_integrated = function(seurat_obj, num_dim) {

  s_obj = seurat_obj

  num_dim = as.integer(num_dim)
  if (num_dim < 5) stop("too few dims: ", num_dim)
  if (num_dim > 50) stop("too many dims: ", num_dim)

  message("\n\n ========== Seurat::ScaleData() ========== \n\n")

  # s_obj = ScaleData(s_obj, features = rownames(s_obj), verbose = FALSE)
  s_obj = ScaleData(s_obj, verbose = FALSE)

  message("\n\n ========== Seurat::PCA() ========== \n\n")

  # PCA on the scaled data
  s_obj = RunPCA(s_obj, npcs = num_dim, verbose = FALSE)

  # plot the output of PCA analysis (shuffle cells so any one group does not appear overrepresented due to ordering)
  pca_plot =
    DimPlot(
      s_obj, reduction = "pca", cells = sample(colnames(s_obj)), group.by = "orig.ident",
      pt.size = 0.5, cols = colors_samples
    ) +
    theme(aspect.ratio = 1)
  ggsave("variance.pca.png", plot = pca_plot, width = 8, height = 6, units = "in")

  message("\n\n ========== Seurat::RunTSNE() ========== \n\n")

  # use tSNE as a tool to visualize, not for clustering directly on tSNE components
  # cells within the graph-based clusters determined above should co-localize on the tSNE plot
  s_obj = RunTSNE(s_obj, reduction = "pca", dims.use = 1:num_dim)

  # reduce point size for larger datasets
  dr_pt_size = get_dr_point_size(s_obj)

  # tSNE using original sample names (shuffle cells so any one group does not appear overrepresented due to ordering)
  s_obj = set_identity(seurat_obj = s_obj, group_var = "orig.ident")
  plot_tsne =
    DimPlot(s_obj, reduction = "tsne", cells = sample(colnames(s_obj)), pt.size = dr_pt_size, cols = colors_samples) +
    theme(aspect.ratio = 1)
  ggsave(glue("dr.tsne.{num_dim}.sample.png"), plot = plot_tsne, width = 8, height = 6, units = "in")
  Sys.sleep(1)
  ggsave(glue("dr.tsne.{num_dim}.sample.pdf"), plot = plot_tsne, width = 8, height = 6, units = "in")
  Sys.sleep(1)

  message("\n\n ========== Seurat::RunUMAP() ========== \n\n")

  # runs the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique
  s_obj = RunUMAP(s_obj, reduction = "pca", dims = 1:num_dim, verbose = FALSE)

  # tSNE using original sample names (shuffle cells so any one group does not appear overrepresented due to ordering)
  s_obj = set_identity(seurat_obj = s_obj, group_var = "orig.ident")
  plot_umap =
    DimPlot(s_obj, reduction = "umap", cells = sample(colnames(s_obj)), pt.size = dr_pt_size, cols = colors_samples) +
    theme(aspect.ratio = 1)
  ggsave(glue("dr.umap.{num_dim}.sample.png"), plot = plot_umap, width = 8, height = 6, units = "in")
  Sys.sleep(1)
  ggsave(glue("dr.umap.{num_dim}.sample.pdf"), plot = plot_umap, width = 8, height = 6, units = "in")
  Sys.sleep(1)

  # compile all cell metadata into a single table
  metadata_tbl = s_obj@meta.data %>% rownames_to_column("cell") %>% as_tibble() %>% mutate(sample_name = orig.ident)
  tsne_tbl = s_obj[["tsne"]]@cell.embeddings %>% round(3) %>% as.data.frame() %>% rownames_to_column("cell")
  umap_tbl = s_obj[["umap"]]@cell.embeddings %>% round(3) %>% as.data.frame() %>% rownames_to_column("cell")
  cells_metadata = metadata_tbl %>% full_join(tsne_tbl, by = "cell") %>% full_join(umap_tbl, by = "cell")
  cells_metadata = cells_metadata %>% arrange(cell)
  write_excel_csv(cells_metadata, path = "metadata.csv")

  return(s_obj)

}

# determine point size for tSNE/UMAP plots (smaller for larger datasets)
get_dr_point_size = function(seurat_obj) {

  pt_size = 1.8
  if (ncol(seurat_obj) > 1000) pt_size = 1.2
  if (ncol(seurat_obj) > 5000) pt_size = 1.0
  if (ncol(seurat_obj) > 10000) pt_size = 0.8
  if (ncol(seurat_obj) > 25000) pt_size = 0.6

  return(pt_size)

}

# perform graph-based clustering and tSNE
calculate_clusters = function(seurat_obj, num_dim, reduction_type = "pca") {

  # check if number of dimensions seems reasonable
  if (num_dim < 5) stop("too few dims: ", num_dim)
  if (num_dim > 50) stop("too many dims: ", num_dim)

  s_obj = seurat_obj

  message("\n\n ========== Seurat::RunTSNE() ========== \n\n")

  # use tSNE as a tool to visualize, not for clustering directly on tSNE components
  # cells within the graph-based clusters determined above should co-localize on the tSNE plot
  s_obj = RunTSNE(s_obj, reduction = reduction_type, dims.use = 1:num_dim)

  # reduce point size for larger datasets
  dr_pt_size = get_dr_point_size(s_obj)

  # tSNE using original sample names (shuffle cells so any one group does not appear overrepresented due to ordering)
  s_obj = set_identity(seurat_obj = s_obj, group_var = "orig.ident")
  plot_tsne =
    DimPlot(s_obj, reduction = "tsne", cells = sample(colnames(s_obj)), pt.size = dr_pt_size, cols = colors_samples) +
    theme(aspect.ratio = 1)
  ggsave(glue("dr.tsne.{num_dim}.sample.png"), plot = plot_tsne, width = 8, height = 6, units = "in")
  Sys.sleep(1)
  ggsave(glue("dr.tsne.{num_dim}.sample.pdf"), plot = plot_tsne, width = 8, height = 6, units = "in")
  Sys.sleep(1)

  message("\n\n ========== Seurat::RunUMAP() ========== \n\n")

  # runs the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique
  s_obj = RunUMAP(s_obj, reduction = "pca", dims = 1:num_dim, verbose = FALSE)

  # tSNE using original sample names (shuffle cells so any one group does not appear overrepresented due to ordering)
  s_obj = set_identity(seurat_obj = s_obj, group_var = "orig.ident")
  plot_umap =
    DimPlot(s_obj, reduction = "umap", cells = sample(colnames(s_obj)), pt.size = dr_pt_size, cols = colors_samples) +
    theme(aspect.ratio = 1)
  ggsave(glue("dr.umap.{num_dim}.sample.png"), plot = plot_umap, width = 8, height = 6, units = "in")
  Sys.sleep(1)
  ggsave(glue("dr.umap.{num_dim}.sample.pdf"), plot = plot_umap, width = 8, height = 6, units = "in")
  Sys.sleep(1)

  message("\n\n ========== Seurat::FindNeighbors() ========== \n\n")

  message("reduction type: ", reduction_type)
  message("num dims: ", num_dim)

  # construct a Shared Nearest Neighbor (SNN) Graph for a given dataset
  s_obj = FindNeighbors(s_obj, reduction = reduction_type, dims = 1:num_dim, graph.name = "snn", compute.SNN = TRUE, force.recalc = TRUE)

  message("\n\n ========== Seurat::FindClusters() ========== \n\n")

  message("initial metadata fields: ", str_c(colnames(s_obj@meta.data), collapse = ", "))

  # resolutions for graph-based clustering
  # increased resolution values lead to more clusters (recommendation: 0.6-1.2 for 3K cells, 2-4 for 33K cells)
  res_range = seq(0.1, 2.5, 0.1)
  if (ncol(s_obj) > 1000) res_range = c(res_range, 3, 4, 5, 6, 7, 8)

  # algorithm: 1 = original Louvain; 2 = Louvain with multilevel refinement; 3 = SLM
  # identify clusters of cells by SNN modularity optimization based clustering algorithm
  s_obj = FindClusters(s_obj, algorithm = 3, resolution = res_range, graph.name = "snn", verbose = FALSE)

  message("new metadata fields: ", str_c(colnames(s_obj@meta.data), collapse = ", "))

  # PrintFindClustersParams to print a nicely formatted formatted summary of the parameters that were chosen

  # create a separate sub-directory for cluster resolution plots
  clusters_dir = "clusters-resolutions"
  if (!dir.exists(clusters_dir)) dir.create(clusters_dir)

  # for calculated cluster resolutions: remove redundant (same number of clusters), rename, and plot
  res_cols = str_subset(colnames(s_obj@meta.data), "snn_res")
  res_cols = sort(res_cols)
  res_num_clusters_prev = 1
  for (res in res_cols) {

    # proceed if current resolution has more clusters than previous and less than 80 (limited by the color scheme)
    res_vector = s_obj@meta.data[, res] %>% as.character()
    res_num_clusters_cur = res_vector %>% n_distinct()
    if (res_num_clusters_cur > res_num_clusters_prev && res_num_clusters_cur < 80) {

      # check if the resolution still has original labels (characters starting with 0)
      if (min(res_vector) == "0") {

        # convert to character vector
        s_obj@meta.data[, res] = as.character(s_obj@meta.data[, res])
        # relabel identities so they start with 1 and not 0
        s_obj@meta.data[, res] = as.numeric(s_obj@meta.data[, res]) + 1
        # pad with 0s to avoid sorting issues
        s_obj@meta.data[, res] = str_pad(s_obj@meta.data[, res], width = 2, side = "left", pad = "0")
        # pad with "C" to avoid downstream numeric conversions
        s_obj@meta.data[, res] = str_c("C", s_obj@meta.data[, res])
        # encode as a factor
        s_obj@meta.data[, res] = factor(s_obj@meta.data[, res])

      }

      # resolution value based on resolution column name
      res_val = sub("snn_res\\.", "", res)

      # plot file name
      res_str = gsub("\\.", "", res)
      dr_filename = glue("{clusters_dir}/dr.{reduction_type}.{num_dim}.{res_str}.clust{res_num_clusters_cur}")

      s_obj = plot_clusters(seurat_obj = s_obj, resolution = res_val, filename_base = dr_filename)

      # add blank line to make output easier to read
      message(" ")

    } else {

      # remove resolution if the number of clusters is same as previous
      s_obj@meta.data = s_obj@meta.data %>% select(-one_of(res))

    }

    # update resolution cluster count for next iteration
    res_num_clusters_prev = res_num_clusters_cur

  }

  message("updated metadata fields: ", str_c(colnames(s_obj@meta.data), collapse = ", "))

  # compile all cell metadata into a single table
  metadata_tbl = s_obj@meta.data %>% rownames_to_column("cell") %>% as_tibble() %>% mutate(sample_name = orig.ident)
  tsne_tbl = s_obj[["tsne"]]@cell.embeddings %>% round(3) %>% as.data.frame() %>% rownames_to_column("cell")
  umap_tbl = s_obj[["umap"]]@cell.embeddings %>% round(3) %>% as.data.frame() %>% rownames_to_column("cell")
  cells_metadata = metadata_tbl %>% full_join(tsne_tbl, by = "cell") %>% full_join(umap_tbl, by = "cell")
  cells_metadata = cells_metadata %>% arrange(cell)
  write_excel_csv(cells_metadata, path = "metadata.csv")

  return(s_obj)

}

# plot tSNE with color-coded clusters at specified resolution
plot_clusters = function(seurat_obj, resolution, filename_base) {

  s_obj = seurat_obj

  # set identities based on specified resolution
  s_obj = set_identity(seurat_obj = s_obj, group_var = resolution)

  # print stats
  num_clusters = Idents(s_obj) %>% as.character() %>% n_distinct()
  message("resolution: ", resolution)
  message("num clusters: ", num_clusters)

  # generate plot if there is a reasonable number of clusters
  if (num_clusters > 1 && num_clusters < 80) {

    # shuffle cells so they appear randomly and one group does not show up on top
    plot_tsne =
      DimPlot(
        s_obj, reduction = "tsne", cells = sample(colnames(s_obj)),
        pt.size = get_dr_point_size(s_obj), cols = colors_clusters
      ) +
      theme(aspect.ratio = 1)
    ggsave(glue("{filename_base}.tsne.png"), plot = plot_tsne, width = 8, height = 6, units = "in")
    Sys.sleep(1)
    ggsave(glue("{filename_base}.tsne.pdf"), plot = plot_tsne, width = 8, height = 6, units = "in")
    Sys.sleep(1)

    plot_umap =
      DimPlot(
        s_obj, reduction = "umap", cells = sample(colnames(s_obj)),
        pt.size = get_dr_point_size(s_obj), cols = colors_clusters
      ) +
      theme(aspect.ratio = 1)
    ggsave(glue("{filename_base}.umap.png"), plot = plot_umap, width = 8, height = 6, units = "in")
    Sys.sleep(1)
    ggsave(glue("{filename_base}.umap.pdf"), plot = plot_umap, width = 8, height = 6, units = "in")
    Sys.sleep(1)

    if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")

  }

  return(s_obj)

}

# check grouping variable/resolution against existing meta data columns
check_group_var = function(seurat_obj, group_var) {

  s_obj = seurat_obj

  # check if the grouping variable is one of meta data columns
  if (!(group_var %in% colnames(s_obj@meta.data))) {

    # check if grouping variable is the resolution value (X.X instead of res.X.X)
    res_column = str_c("snn_res.", group_var)
    if (res_column %in% colnames(s_obj@meta.data)) {
      group_var = res_column
    } else {
      stop("unknown grouping variable: ", group_var)
    }

  }

  return(group_var)

}

# set identity based on a specified variable/resolution
set_identity = function(seurat_obj, group_var) {

  s_obj = seurat_obj

  group_var = check_group_var(seurat_obj = s_obj, group_var = group_var)

  # set identities based on selected grouping variable
  message("setting grouping variable: ", group_var)
  Idents(s_obj) = group_var

  return(s_obj)

}

# plot a set of genes
plot_genes = function(seurat_obj, genes, name) {

  # color gradient for FeaturePlot-based plots
  gradient_colors = c("gray85", "red2")

  # tSNE plots color-coded by expression level (should be square to match the original tSNE plots)
  feat_plot =
    FeaturePlot(
      seurat_obj, features = genes, reduction = "tsne", cells = sample(colnames(seurat_obj)),
      pt.size = 0.5, cols = gradient_colors, ncol = 4
    )
  ggsave(glue("{name}.tsne.png"), plot = feat_plot, width = 16, height = 10, units = "in")
  ggsave(glue("{name}.tsne.pdf"), plot = feat_plot, width = 16, height = 10, units = "in")

  # dot plot visualization
  dot_plot = DotPlot(seurat_obj, features = genes, dot.scale = 12, cols = gradient_colors)
  ggsave(glue("{name}.dotplot.png"), plot = dot_plot, width = 20, height = 8, units = "in")
  ggsave(glue("{name}.dotplot.pdf"), plot = dot_plot, width = 20, height = 8, units = "in")

  # gene violin plots (size.use below 0.2 doesn't seem to make a difference)
  # skip PDF since every cell has to be plotted and they become too big
  vln_plot = VlnPlot(seurat_obj, features = genes, pt.size = 0.1, combine = TRUE, cols = colors_clusters, ncol = 4)
  ggsave(glue("{name}.violin.png"), plot = vln_plot, width = 16, height = 10, units = "in")

  # expression levels per cluster for bar plots (averaging and output are in non-log space)
  cluster_avg_exp = AverageExpression(seurat_obj, assay = "RNA", features = genes, verbose = FALSE)[["RNA"]]
  cluster_avg_exp_long = cluster_avg_exp %>% rownames_to_column("gene") %>% gather(cluster, avg_exp, -gene)

  # bar plots
  # create a named color scheme to ensure names and colors are in the proper order
  clust_names = levels(Idents(seurat_obj))
  color_scheme_named = colors_clusters[1:length(clust_names)]
  names(color_scheme_named) = clust_names
  barplot_plot = ggplot(cluster_avg_exp_long, aes(x = cluster, y = avg_exp, fill = cluster)) +
    geom_col(color = "black") +
    theme(legend.position = "none") +
    scale_fill_manual(values = color_scheme_named) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_cowplot() +
    facet_wrap(~ gene, ncol = 4, scales = "free")
  ggsave(glue("{name}.barplot.png"), plot = barplot_plot, width = 16, height = 10, units = "in")
  ggsave(glue("{name}.barplot.pdf"), plot = barplot_plot, width = 16, height = 10, units = "in")

}

# calculate cluster stats (number of cells, average expression, cell-gene matrix)
calculate_cluster_stats = function(seurat_obj, label) {

  message("\n\n ========== calculate cluster stats ========== \n\n")

  message("cluster names: ", str_c(levels(Idents(seurat_obj)), collapse = ", "))

  # compile relevant cell metadata into a single table
  seurat_obj$cluster = Idents(seurat_obj)
  metadata_tbl = seurat_obj@meta.data %>% rownames_to_column("cell") %>% as_tibble() %>%
    select(cell, num_UMIs, num_genes, pct_mito, sample_name = orig.ident, cluster)
  tsne_tbl = seurat_obj[["tsne"]]@cell.embeddings %>% round(3) %>% as.data.frame() %>% rownames_to_column("cell")
  umap_tbl = seurat_obj[["umap"]]@cell.embeddings %>% round(3) %>% as.data.frame() %>% rownames_to_column("cell")
  cells_metadata = metadata_tbl %>% full_join(tsne_tbl, by = "cell") %>% full_join(umap_tbl, by = "cell")
  cells_metadata = cells_metadata %>% arrange(cell)
  write_excel_csv(cells_metadata, path = glue("metadata.{label}.csv"))

  # get number of cells split by cluster and by sample
  summary_cluster_sample =
    cells_metadata %>%
    select(cluster, sample_name) %>%
    mutate(num_cells_total = n()) %>%
    group_by(sample_name) %>%
    mutate(num_cells_sample = n()) %>%
    group_by(cluster) %>%
    mutate(num_cells_cluster = n()) %>%
    group_by(cluster, sample_name) %>%
    mutate(num_cells_cluster_sample = n()) %>%
    ungroup() %>%
    distinct() %>%
    mutate(
      pct_cells_cluster = num_cells_cluster / num_cells_total,
      pct_cells_cluster_sample = num_cells_cluster_sample / num_cells_sample
    ) %>%
    mutate(
      pct_cells_cluster = round(pct_cells_cluster * 100, 1),
      pct_cells_cluster_sample = round(pct_cells_cluster_sample * 100, 1)
    ) %>%
    arrange(cluster, sample_name)

  # get number of cells split by cluster (ignore samples)
  summary_cluster = summary_cluster_sample %>% select(-contains("sample")) %>% distinct()
  write_excel_csv(summary_cluster, path = glue("summary.{label}.csv"))

  # gene expression for an "average" cell in each identity class (averaging and output are in non-log space)
  cluster_avg_exp = AverageExpression(seurat_obj, assay = "RNA", verbose = FALSE)[["RNA"]]
  cluster_avg_exp = cluster_avg_exp %>% round(3) %>% rownames_to_column("gene") %>% arrange(gene)
  write_excel_csv(cluster_avg_exp, path = glue("expression.mean.{label}.csv"))

  Sys.sleep(1)

  # export results split by sample if multiple samples are present
  num_samples = cells_metadata %>% pull(sample_name) %>% n_distinct()
  if (num_samples > 1) {

    # number of cells split by cluster and by sample
    write_excel_csv(summary_cluster_sample, path = glue("summary.{label}.per-sample.csv"))

    # cluster averages split by sample
    sample_avg_exp = AverageExpression(seurat_obj, assay = "RNA", add.ident = "orig.ident", verbose = FALSE)[["RNA"]]
    sample_avg_exp = sample_avg_exp %>% round(3) %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
    write_excel_csv(sample_avg_exp, path = glue("expression.mean.{label}.per-sample.csv"))

  }

}

# calculate cluster markers (compared to all other cells) and plot top ones
# tests:
# - roc: ROC test returns the classification power (ranging from 0 - random, to 1 - perfect)
# - wilcox: Wilcoxon rank sum test (default in Seurat 2)
# - bimod: Likelihood-ratio test for single cell gene expression (McDavid, Bioinformatics, 2013) (default in Seurat 1)
# - tobit: Tobit-test for differential gene expression (Trapnell, Nature Biotech, 2014)
# - MAST: GLM-framework that treates cellular detection rate as a covariate (Finak, Genome Biology, 2015)
# pairwise option compares each cluster to each of the other clusters to yield markers that are both local and global
calculate_cluster_markers = function(seurat_obj, label, test, pairwise = FALSE) {

  message("\n\n ========== calculate cluster markers ========== \n\n")

  message("cluster set: ", label)
  message("marker test: ", test)

  # get cluster names
  clusters = Idents(seurat_obj) %>% as.character() %>% unique() %>% sort()

  # use only clusters with more than 10 cells
  clusters = clusters[table(Idents(seurat_obj)) > 10]

  if (!pairwise) {

    # standard cluster markers calculation

    markers_dir = "markers-global"

    # capture output to avoid excessive warnings
    markers_log =
      capture.output({
        all_markers =
          FindAllMarkers(
            seurat_obj, assay = "RNA", test.use = test, logfc.threshold = log(1.2), min.pct = 0.2,
            only.pos = FALSE, min.diff.pct = -Inf, verbose = FALSE
          )
      }, type = "message")

    # do some light filtering and clean up (ROC test returns slightly different output)
    if (test == "roc") {

      all_markers =
        all_markers %>%
        select(cluster, gene, avg_logFC, myAUC, power) %>%
        filter(power > 0.4) %>%
        mutate(avg_logFC = round(avg_logFC, 3), myAUC = round(myAUC, 3), power = round(power, 3)) %>%
        arrange(cluster, -power)
      top_markers = all_markers %>% filter(avg_logFC > 0) %>% group_by(cluster) %>% top_n(50, power)

    } else {

      all_markers =
        all_markers %>%
        select(cluster, gene, avg_logFC, p_val, p_val_adj) %>%
        filter(p_val_adj < 0.001) %>%
        mutate(avg_logFC = round(avg_logFC, 3)) %>%
        arrange(cluster, p_val_adj, p_val)
      top_markers = all_markers %>% filter(avg_logFC > 0) %>% group_by(cluster) %>% top_n(50, avg_logFC)

    }

  } else {

    # pairwise (each cluster versus each other cluster) cluster markers calculation

    markers_dir = "markers-pairwise"

    # initialize empty results tibble
    unfiltered_markers = tibble(
      cluster = character(),
      cluster2 = character(),
      gene = character(),
      avg_logFC = numeric(),
      p_val = numeric(),
      p_val_adj = numeric()
    )

    # check each cluster combination
    for (cluster1 in clusters) {
      for (cluster2 in setdiff(clusters, cluster1)) {

        # find differentially expressed genes between two specific clusters
        # low fold change cutoff to maximize chance of appearing in all comparisons
        # capture output to avoid excessive warnings
        markers_log =
          capture.output({
            cur_markers =
              FindMarkers(
                seurat_obj, assay = "RNA", ident.1 = cluster1, ident.2 = cluster2, test.use = test,
                logfc.threshold = log(1.1), min.pct = 0.1,
                only.pos = TRUE, min.diff.pct = -Inf, verbose = FALSE
              )
          }, type = "message")

        # clean up markers table (would need to be modified for "roc" test)
        cur_markers =
          cur_markers %>%
          rownames_to_column("gene") %>%
          mutate(cluster = cluster1) %>%
          mutate(cluster2 = cluster2) %>%
          filter(p_val_adj < 0.01) %>%
          mutate(avg_logFC = round(avg_logFC, 3)) %>%
          select(one_of(colnames(unfiltered_markers)))

        # add current cluster combination genes to the table of all markers
        unfiltered_markers = bind_rows(unfiltered_markers, cur_markers)

      }
    }

    # adjust test name for output
    test = glue("pairwise.{test}")

    # sort the markers to make the table more readable
    unfiltered_markers =
      unfiltered_markers %>%
      distinct() %>%
      add_count(cluster, gene) %>%
      rename(cluster_gene_n = n) %>%
      arrange(cluster, gene, cluster2)

    # filter for genes that are significant compared to all other clusters
    all_markers =
      unfiltered_markers %>%
      filter(cluster_gene_n == (length(clusters) - 1)) %>%
      select(-cluster_gene_n)

    # extract the lowest and highest fold changes and p-values
    all_markers =
      all_markers %>%
      group_by(cluster, gene) %>%
      summarize_at(
        c("avg_logFC", "p_val", "p_val_adj"),
        funs(min, max)
      ) %>%
      ungroup() %>%
      arrange(cluster, -avg_logFC_min)
    all_markers

    top_markers = all_markers %>% group_by(cluster) %>% top_n(50, avg_logFC_min)

  }

  # create a separate sub-directory for all markers
  if (!dir.exists(markers_dir)) dir.create(markers_dir)

  # save unfiltered markers for pairwise comparisons
  if (pairwise) {
    unfiltered_markers_csv = glue("{markers_dir}/markers.{label}.{test}.unfiltered.csv")
    message("unfiltered markers: ", unfiltered_markers_csv)
    write_excel_csv(unfiltered_markers, path = unfiltered_markers_csv)
    Sys.sleep(1)
  }

  all_markers_csv = glue("{markers_dir}/markers.{label}.{test}.all.csv")
  message("all markers: ", all_markers_csv)
  write_excel_csv(all_markers, path = all_markers_csv)
  Sys.sleep(1)

  top_markers_csv = glue("{markers_dir}/markers.{label}.{test}.top.csv")
  message("top markers: ", top_markers_csv)
  write_excel_csv(top_markers, path = top_markers_csv)
  Sys.sleep(1)

  # get marker genes for each cluster
  for (cluster_name in clusters) {

    # plot top genes if enough were identified
    filename_label = glue("{markers_dir}/markers.{label}-{cluster_name}.{test}")
    cluster_markers = top_markers %>% filter(cluster == cluster_name)
    if (nrow(cluster_markers) > 9) {
      Sys.sleep(1)
      top_cluster_markers = cluster_markers %>% head(12) %>% pull(gene)
      plot_genes(seurat_obj, genes = top_cluster_markers, name = filename_label)
    }

  }

}

# calculate differentially expressed genes within each cluster
calculate_cluster_de_genes = function(seurat_obj, label, test) {

  message("\n\n ========== calculate cluster DE genes ========== \n\n")

  # create a separate sub-directory for differential expression results
  de_dir = "diff-expression"
  if (!dir.exists(de_dir)) dir.create(de_dir)

  # common settings
  num_de_genes = 50

  # cluster names
  clusters = seurat_obj@ident %>% as.character() %>% unique() %>% sort()

  # get DE genes for each cluster
  for (clust_name in clusters) {

    message(glue("calculating DE genes for cluster {clust_name}"))

    # subset to the specific cluster
    clust_obj = SubsetData(seurat_obj, ident.use = clust_name)

    # revert back to original sample/library labels
    clust_obj = SetAllIdent(clust_obj, id = "orig.ident")

    message("cluster cells: ", ncol(clust_obj@data))
    message("cluster groups: ", paste(levels(clust_obj@ident), collapse = ", "))

    # continue if cluster has multiple groups and more than 100 cells and more than 10 cells per group
    if (length(unique(clust_obj@ident)) > 1 && ncol(clust_obj@data) > 50 && min(table(clust_obj@ident)) > 10) {

      # iterate through sample/library combinations (relevant if more than two)
      group_combinations = combn(levels(clust_obj@ident), m = 2, simplify = TRUE)
      for (combination_num in 1:ncol(group_combinations)) {

        # determine combination
        s1 = group_combinations[1, combination_num]
        s2 = group_combinations[2, combination_num]
        comparison_label = glue("{s1}-vs-{s2}")
        message(glue("comparison: {clust_name} {comparison_label}"))

        filename_label = glue("{de_dir}/de.{label}-{clust_name}.{comparison_label}.{test}")

        # find differentially expressed genes (default Wilcoxon rank sum test)
        de_genes = FindMarkers(clust_obj, ident.1 = s1, ident.2 = s2, test.use = test,
                               logfc.threshold = log(1.1), min.pct = 0.1,
                               only.pos = FALSE, print.bar = FALSE)

        # do some light filtering and clean up
        de_genes = de_genes %>%
          rownames_to_column("gene") %>%
          select(gene, avg_logFC, p_val, p_val_adj) %>%
          filter(p_val < 0.01) %>%
          mutate(avg_logFC = round(avg_logFC, 3)) %>%
          arrange(p_val_adj, p_val)

        message(glue("num DE genes (p<0.01): {nrow(de_genes)}"))

        # save stats table
        write_excel_csv(de_genes, path = glue("{filename_label}.stats.csv"))

        # heatmap of top genes if any significant genes are present
        if (nrow(de_genes) > 5) {
          top_de_genes = de_genes %>% top_n(num_de_genes, -p_val_adj) %>% arrange(avg_logFC) %>% pull(gene)
          plot_hm = DoHeatmap(clust_obj, genes.use = top_de_genes,
                              use.scaled = TRUE, remove.key = FALSE, slim.col.label = TRUE,
                              col.low = "lemonchiffon", col.mid = "gold2", col.high = "red3",
                              cex.row = 10, cex.col = 0.5, group.cex = 15, disp.min = -2, disp.max = 2)
          heatmap_prefix = glue("{filename_label}.heatmap.top{num_de_genes}")
          ggsave(glue("{heatmap_prefix}.png"), plot = plot_hm, width = 15, height = 10, units = "in")
          Sys.sleep(1)
          ggsave(glue("{heatmap_prefix}.pdf"), plot = plot_hm, width = 15, height = 10, units = "in")
          Sys.sleep(1)
        }

      }

    } else {

      message("skip cluster: ", clust_name)

    }

    message(" ")

  }

}


# ========== main ==========


# output width
options(width = 120)
# print warnings as they occur
options(warn = 1)
# default type for the bitmap devices such as png (should default to "cairo")
options(bitmapType = "cairo")

# retrieve the command-line arguments
suppressPackageStartupMessages(library(docopt))
opts = docopt(doc)

# show docopt options
# print(opts)

# dependencies
load_libraries()

# evaluate R expressions asynchronously when possible (such as ScaleData)
plan("multiprocess", workers = 4)
# increase the limit of the data to be shuttled between the processes from default 500MB to 50GB
options(future.globals.maxSize = 50e9)

# global settings
colors_samples = c(brewer.pal(5, "Set1"), brewer.pal(8, "Dark2"), pal_igv("default")(51))
colors_clusters = c(pal_d3("category10")(10), pal_d3("category20b")(20), pal_igv("default")(51))

# analysis info
analysis_step = "unknown"
out_dir = opts$analysis_dir

# create analysis directory if starting new analysis or exit if analysis already exists
if (opts$create || opts$combine || opts$integrate) {

  if (opts$create) analysis_step = "create"
  if (opts$combine) analysis_step = "combine"
  if (opts$integrate) analysis_step = "integrate"
  message(glue("\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"))

  if (dir.exists(out_dir)) {
    stop(glue("output analysis dir {out_dir} already exists"))
  } else {
    dir.create(out_dir)
  }

  # original working dir (before moving to analysis dir)
  original_wd = getwd()

}

# set analysis directory as working directory
if (dir.exists(out_dir)) {
  setwd(out_dir)
} else {
  stop(glue("output analysis dir {out_dir} does not exist"))
}

# check which command was used
if (opts$create) {

  # log to file
  write(glue("analysis: {out_dir}"), file = "create.log", append = TRUE)
  write(glue("seurat version: {packageVersion('Seurat')}"), file = "create.log", append = TRUE)

  # create new seurat object based on input sample names and sample directories
  # can work with multiple samples, but the appropriate way is to use "combine" with objects that are pre-filtered
  seurat_obj = load_sample_counts_matrix(opts$sample_name, opts$sample_dir)

  # filter by number of genes and mitochondrial genes percentage (optional parameters)
  seurat_obj = filter_data(seurat_obj, min_genes = opts$min_genes, max_genes = opts$max_genes, max_mt = opts$mt)

  # calculate various variance metrics
  seurat_obj = calculate_variance(seurat_obj)

  saveRDS(seurat_obj, file = "seurat_obj.rds")

} else if (opts$combine) {

  # merge multiple samples/libraries based on previous analysis directories
  seurat_obj = combine_seurat_obj(original_wd = original_wd, sample_analysis_dirs = opts$sample_analysis_dir)

  saveRDS(seurat_obj, file = "seurat_obj.rds")

  # calculate various variance metrics
  seurat_obj = calculate_variance(seurat_obj)

  saveRDS(seurat_obj, file = "seurat_obj.rds")

} else if (opts$integrate) {

  # run integration
  seurat_obj = integrate_seurat_obj(original_wd, sample_analysis_dirs = opts$batch_analysis_dir, num_dim = opts$num_dim)
  seurat_obj = calculate_variance_integrated(seurat_obj, num_dim = opts$num_dim)

  saveRDS(seurat_obj, file = "seurat_obj.rds")

} else {

  # all commands besides "create" and "cca" start with an existing seurat object
  if (file.exists("seurat_obj.rds")) {

    message("loading seurat_obj")
    seurat_obj = readRDS("seurat_obj.rds")

  } else {

    stop("seurat obj does not already exist (run 'create' step first)")

  }

  if (opts$cluster) {

    analysis_step = "cluster"
    message(glue("\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"))

    # determine clusters
    seurat_obj = calculate_clusters(seurat_obj, num_dim = as.integer(opts$num_pcs), reduction_type = "pca")
    saveRDS(seurat_obj, file = "seurat_obj.rds")

  }

  if (opts$identify || opts$de) {

    # set resolution in the seurat object
    group_var = check_group_var(seurat_obj = seurat_obj, group_var = opts$resolution)
    seurat_obj = set_identity(seurat_obj = seurat_obj, group_var = group_var)

    # use a grouping-specific sub-directory for all output
    grouping_label = gsub("\\.", "", group_var)
    num_clusters = Idents(seurat_obj) %>% as.character() %>% n_distinct()
    clust_label = glue("clust{num_clusters}")
    res_dir = glue("clusters-{grouping_label}-{clust_label}")
    if (!dir.exists(res_dir)) dir.create(res_dir)
    setwd(res_dir)

    if (opts$identify) {

      analysis_step = "identify"
      message(glue("\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"))

      # create tSNE plot (should already exist in the main directory)
      dr_filename = glue("dr.{grouping_label}.{clust_label}")
      seurat_obj = plot_clusters(seurat_obj, resolution = opts$resolution, filename_base = dr_filename)

      # cluster stat tables (number of cells and average expression)
      calculate_cluster_stats(seurat_obj, label = clust_label)

      # calculate and plot standard cluster markers
      calculate_cluster_markers(seurat_obj, label = clust_label, test = "roc")
      calculate_cluster_markers(seurat_obj, label = clust_label, test = "wilcox")
      calculate_cluster_markers(seurat_obj, label = clust_label, test = "MAST")

      # calculate and plot pairwise cluster markers (very slow, so skip for high number of clusters)
      num_clusters = Idents(seurat_obj) %>% as.character() %>% n_distinct()
      if (num_clusters < 20) {
        calculate_cluster_markers(seurat_obj, label = clust_label, test = "wilcox", pairwise = TRUE)
        calculate_cluster_markers(seurat_obj, label = clust_label, test = "MAST", pairwise = TRUE)
      }

    }

    # check if there are multiple samples for differential expression analysis
    num_samples = seurat_obj$orig.ident %>% n_distinct()
    if (num_samples > 1) {

      # differential expression
      if (opts$de) {

          analysis_step = "diff"
          message(glue("\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"))

          calculate_cluster_de_genes(seurat_obj, label = clust_label, test = "wilcox")
          calculate_cluster_de_genes(seurat_obj, label = clust_label, test = "bimod")
          calculate_cluster_de_genes(seurat_obj, label = clust_label, test = "MAST")

      }

    }

  }

}

message(glue("\n\n ========== finished analysis step {analysis_step} for {out_dir} ========== \n\n"))

# delete Rplots.pdf
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")



# end
