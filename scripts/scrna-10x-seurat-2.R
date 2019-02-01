#!/usr/bin/env Rscript


"
Analysis of 10x Genomics Chromium single cell RNA-seq data using Seurat (version 2.3) starting with Cell Ranger output.

Basic workflow steps:
  1 - create - import counts matrix, perform initial QC, and calculate various variance metrics (slowest step)
  2 - cluster - perform clustering based on number of PCs
  3 - identify - identify clusters based on specified clustering/resolution (higher resolution for more clusters)

Optional steps:
  combine - merge multiple samples/libraries
  cca - perform CCA for alignment (batch correction) across multiple sample batches
  align - perform alignment based on number of CCs (after performing CCA)
  de - differential expression between samples/libraries within clusters

Usage:
  scrna-10x-seurat-2.R create <analysis_dir> <sample_name> <sample_dir> [--min_genes=<n> --max_genes=<n> --mt=<n>]
  scrna-10x-seurat-2.R cluster <analysis_dir> <num_pcs>
  scrna-10x-seurat-2.R identify <analysis_dir> <resolution>
  scrna-10x-seurat-2.R combine <analysis_dir> <sample_analysis_dir>...
  scrna-10x-seurat-2.R cca <analysis_dir> <batch_analysis_dir>...
  scrna-10x-seurat-2.R align <analysis_dir> <num_ccs>
  scrna-10x-seurat-2.R de <analysis_dir> <resolution>
  scrna-10x-seurat-2.R --help

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
    library(Matrix)
    library(tidyverse)
    library(cowplot)
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

    # determine counts matrix directory
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
    colnames(counts_matrix) = paste0(sample_name, ":", colnames(counts_matrix))

    # combine current matrix with previous
    if (i == 1) {

      # skip if there is no previous matrix
      merged_counts_matrix = counts_matrix

    } else {

      # check if genes are the same for current and previous matrices
      if (!identical(rownames(merged_counts_matrix), rownames(counts_matrix))) {

        # generate a warning, since this is probably a mistake
        warning("counts matrix genes are not the same for different libraries")
        Sys.sleep(5)

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
  counts_matrix_filename = "counts.raw.txt"
  write.table(as.matrix(counts_matrix), file = counts_matrix_filename, quote = FALSE, sep = "\t", col.names = NA)
  system(paste0("gzip ", counts_matrix_filename))

  # save counts matrix as a csv file (to be consistent with the rest of the tables)
  raw_data = counts_matrix %>% as.matrix() %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
  write_csv(raw_data, path = "counts.raw.csv.gz")
  rm(raw_data)

  message("\n\n ========== create seurat object ========== \n\n")

  if (is.null(proj_name)) {

    # if name is not set, then it's a manually merged counts matrix
    # save.raw parameter causes errors in v2
    s_obj = CreateSeuratObject(raw.data = counts_matrix, min.cells = 5, min.genes = 250, project = "proj",
                               names.field = 1, names.delim = ":")
    rm(counts_matrix)

  } else if (proj_name == "aggregated") {

    # multiple libraries combined using Cell Ranger (cellranger aggr)

    # setup taking into consideration aggregated names delimiter
    s_obj = CreateSeuratObject(raw.data = counts_matrix, min.cells = 5, min.genes = 250, project = proj_name,
                               names.field = 2, names.delim = "-")

    # import cellranger aggr sample sheet
    sample_sheet_csv = paste0(sample_dir, "/outs/aggregation_csv.csv")
    sample_sheet = read.csv(sample_sheet_csv, stringsAsFactors = FALSE)
    message("samples: ", paste(sample_sheet[, 1], collapse=", "))

    # change s_obj@meta.data$orig.ident sample identities from numbers to names
    s_obj@meta.data$orig.ident = factor(sample_sheet[s_obj@meta.data$orig.ident, 1])
    # set s_obj@ident to the new s_obj@meta.data$orig.ident
    s_obj = SetAllIdent(s_obj, id = "orig.ident")

  } else {

    stop("project name set to unknown value")

  }

  message(glue("imported cells: {ncol(s_obj@data)}"))
  message(glue("imported genes: {nrow(s_obj@data)}"))
  message(" ")

  # log to file
  write(glue("imported cells: {ncol(s_obj@data)}"), file = "create.log", append = TRUE)
  write(glue("imported genes: {nrow(s_obj@data)}"), file = "create.log", append = TRUE)

  # nGene and nUMI are automatically calculated for every object by Seurat
  # calculate the percentage of mitochondrial genes here and store it in percent.mito using the AddMetaData
  mt_genes = grep("^MT-", rownames(s_obj@data), ignore.case = TRUE, value = TRUE)
  percent_mt = Matrix::colSums(s_obj@raw.data[mt_genes, ]) / Matrix::colSums(s_obj@raw.data)
  percent_mt = round(percent_mt * 100, digits = 3)

  # add columns to object@meta.data, and is a great place to stash QC stats
  s_obj = AddMetaData(s_obj, metadata = percent_mt, col.name = "percent.mito")

  message("\n\n ========== nGene/nUMI/percent_mito plots ========== \n\n")

  dist_unfilt_plot = VlnPlot(s_obj, c("nGene", "nUMI", "percent.mito"), nCol = 3, group.by = "orig.ident",
                             point.size.use = 0.1, x.lab.rot = TRUE, size.title.use = 12, cols.use = colors_samples)
  ggsave("qc.distribution.unfiltered.png", plot = dist_unfilt_plot, width = 10, height = 6, units = "in")
  Sys.sleep(1)

  # check for high mitochondrial percentage or low UMI content
  png("qc.correlations.unfiltered.png", res = 200, width = 10, height = 5, units = "in")
    par(mfrow = c(1, 2))
    GenePlot(s_obj, gene1 = "nUMI", gene2 = "percent.mito", cex.use = 0.5, col.use = colors_samples)
    GenePlot(s_obj, gene1 = "nUMI", gene2 = "nGene", cex.use = 0.5, col.use = colors_samples)
  dev.off()
  Sys.sleep(1)

  # check distribution of gene counts and mitochondrial percentage
  low_quantiles = c(0.05, 0.02, 0.01, 0.001)
  high_quantiles = c(0.95, 0.98, 0.99, 0.999)
  message("nGene low percentiles:")
  s_obj@meta.data$nGene %>% quantile(low_quantiles) %>% round(digits = 1) %>% print()
  message(" ")
  message("nGene high percentiles:")
  s_obj@meta.data$nGene %>% quantile(high_quantiles) %>% round(digits = 1) %>% print()
  message(" ")
  message("percent.mito high percentiles:")
  s_obj@meta.data$percent.mito %>% quantile(high_quantiles) %>% round(digits = 1) %>% print()
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
  write(glue("unfiltered min genes: {min(s_obj@meta.data$nGene)}"), file = "create.log", append = TRUE)
  write(glue("unfiltered max genes: {max(s_obj@meta.data$nGene)}"), file = "create.log", append = TRUE)
  write(glue("unfiltered mean num genes: {mean(s_obj@meta.data$nGene)}"), file = "create.log", append = TRUE)
  write(glue("unfiltered median num genes: {median(s_obj@meta.data$nGene)}"), file = "create.log", append = TRUE)

  # convert arguments to integers (command line arguments end up as characters)
  min_genes = as.numeric(min_genes)
  max_genes = as.numeric(max_genes)
  max_mt = as.numeric(max_mt)

  # default cutoffs (gene numbers rounded to nearest 10)
  # as.numeric() converts NULLs to 0 length numerics, so can't use is.null()
  if (!length(min_genes)) min_genes = s_obj@meta.data$nGene %>% quantile(0.02, names = FALSE) %>% round(-1)
  if (!length(max_genes)) max_genes = s_obj@meta.data$nGene %>% quantile(0.98, names = FALSE) %>% round(-1)
  if (!length(max_mt)) max_mt = 10

  message(glue("min genes cutoff: {min_genes}"))
  message(glue("max genes cutoff: {max_genes}"))
  message(glue("max mitochondrial percentage cutoff: {max_mt}"))
  message(" ")

  # log the cutoffs to file
  write(glue("min genes cutoff: {min_genes}"), file = "create.log", append = TRUE)
  write(glue("max genes cutoff: {max_genes}"), file = "create.log", append = TRUE)
  write(glue("max mitochondrial percentage cutoff: {max_mt}"), file = "create.log", append = TRUE)

  message("imported cells: ", ncol(s_obj@data))
  message("imported genes: ", nrow(s_obj@data))

  s_obj = FilterCells(s_obj,
                      subset.names = c("nGene", "percent.mito"),
                      low.thresholds = c(min_genes, -Inf),
                      high.thresholds = c(max_genes, max_mt))

  message("filtered cells: ", ncol(s_obj@data))
  message("filtered genes: ", nrow(s_obj@data))

  dist_filt_plot = VlnPlot(s_obj, c("nGene", "nUMI", "percent.mito"), nCol = 3, group.by = "orig.ident",
                           point.size.use = 0.1, x.lab.rot = TRUE, size.title.use = 12, cols.use = colors_samples)
  ggsave("qc.distribution.filtered.png", plot = dist_filt_plot, width = 10, height = 6, units = "in")

  # after removing unwanted cells from the dataset, normalize the data
  # LogNormalize:
  # - normalizes the gene expression measurements for each cell by the total expression
  # - multiplies this by a scale factor (10,000 by default)
  # - log-transforms the result
  s_obj = NormalizeData(s_obj, normalization.method = "LogNormalize", scale.factor = 100000, display.progress = FALSE)

  # save counts matrix as a basic gzipped text file
  # object@data stores normalized and log-transformed single cell expression
  # used for visualizations, such as violin and feature plots, most diff exp tests, finding high-variance genes
  counts_norm = s_obj@data %>% as.matrix() %>% round(3)
  counts_norm = counts_norm %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
  write_csv(counts_norm, path = "counts.normalized.csv.gz")

  # log to file
  write(glue("filtered cells: {ncol(s_obj@data)}"), file = "create.log", append = TRUE)
  write(glue("filtered genes: {nrow(s_obj@data)}"), file = "create.log", append = TRUE)
  write(glue("filtered mean num genes: {mean(s_obj@meta.data$nGene)}"), file = "create.log", append = TRUE)
  write(glue("filtered median num genes: {median(s_obj@meta.data$nGene)}"), file = "create.log", append = TRUE)

  return(s_obj)

}

# merge multiple Seurat objects
combine_seurat_obj = function(original_wd, sample_analysis_dirs) {

  if (length(sample_analysis_dirs) < 2) stop("must have at least 2 samples to merge")

  message("\n\n ========== combine samples ========== \n\n")

  for (i in 1:length(sample_analysis_dirs)) {

    sample_analysis_dir = sample_analysis_dirs[i]
    sample_analysis_dir = glue("{original_wd}/{sample_analysis_dir}")
    sample_seurat_rds = glue("{sample_analysis_dir}/seurat_obj.rds")

    # check if analysis dir is valid
    if (!dir.exists(sample_analysis_dir)) stop(glue("dir {sample_analysis_dir} does not exist"))
    # check if seurat object exists
    if (!file.exists(sample_seurat_rds)) stop(glue("seurat object rds {sample_seurat_rds} does not exist"))

    # load seurat object
    single_obj = readRDS(sample_seurat_rds)

    # clean up object
    single_obj@var.genes = vector()
    single_obj@hvg.info = data.frame()
    single_obj@scale.data = matrix()
    single_obj@dr = list()
    single_obj@meta.data = single_obj@meta.data %>% select(-starts_with("res"))

    # print single sample sample stats
    sample_name = single_obj@meta.data[1, "orig.ident"]
    message(glue("sample {sample_name} dir: {basename(sample_analysis_dir)}"))
    write(glue("sample {sample_name} dir: {basename(sample_analysis_dir)}"), file = "create.log", append = TRUE)
    message(glue("sample {sample_name} cells: {ncol(single_obj@data)}"))
    write(glue("sample {sample_name} cells: {ncol(single_obj@data)}"), file = "create.log", append = TRUE)
    message(glue("sample {sample_name} genes: {nrow(single_obj@data)}"))
    write(glue("sample {sample_name} genes: {nrow(single_obj@data)}"), file = "create.log", append = TRUE)
    message(" ")

    # merge samples
    if (i == 1) {

      # copy if first one
      merged_obj = single_obj
      rm(single_obj)

    } else {

      # require more cells with detected expression of each gene for larger datasets
      min_cells = 5
      if (ncol(merged_obj@data) > 10000) min_cells = 10

      # save original normalized data to preserve it (in case custom normalization was used)
      merged_norm_mat = Seurat:::RowMergeSparseMatrices(merged_obj@data, single_obj@data)
      # merge objects
      merged_obj = MergeSeurat(object1 = merged_obj, object2 = single_obj, min.cells = min_cells,
                               names.field = 1, names.delim = ":", do.normalize = FALSE)
      # restore original normalized data
      merged_obj@data = merged_norm_mat[rownames(merged_obj@data), colnames(merged_obj@data)]
      # remove large objects
      rm(single_obj)
      rm(merged_norm_mat)

    }

  }

  # print combined sample stats
  message(glue("combined cells: {ncol(merged_obj@data)}"))
  write(glue("combined cells: {ncol(merged_obj@data)}"), file = "create.log", append = TRUE)
  message(glue("combined genes: {nrow(merged_obj@data)}"))
  write(glue("combined genes: {nrow(merged_obj@data)}"), file = "create.log", append = TRUE)

  # save raw counts matrix
  counts_raw = merged_obj@raw.data %>% as.matrix() %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
  write_csv(counts_raw, path = "counts.raw.csv.gz")
  rm(counts_raw)

  # save counts matrix as a basic gzipped text file
  # object@data stores normalized and log-transformed single cell expression
  counts_norm = merged_obj@data %>% as.matrix() %>% round(3)
  counts_norm = counts_norm %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
  write_csv(counts_norm, path = "counts.normalized.csv.gz")
  rm(counts_norm)

  # create a named color scheme to ensure names and colors are in the proper order
  sample_names = merged_obj@meta.data$orig.ident %>% as.character() %>% sort() %>% unique()
  colors_samples_named = colors_samples[1:length(sample_names)]
  names(colors_samples_named) = sample_names

  dist_plot = VlnPlot(merged_obj, c("nGene", "nUMI", "percent.mito"), nCol = 3, group.by = "orig.ident",
                      point.size.use = 0.1, x.lab.rot = TRUE, size.title.use = 12, size.x.use = 12,
                      do.sort = TRUE, remove.legend = TRUE, cols.use = colors_samples_named)
  ggsave("qc.distribution.png", plot = dist_plot, width = 15, height = 6, units = "in")

  return(merged_obj)

}

# calculate various variance metrics
# PC selection approaches:
# - PCHeatmap - more supervised, exploring PCs to determine relevant sources of heterogeneity
# - PCElbowPlot - heuristic that is commonly used and can be calculated instantly
# - JackStrawPlot - implements a statistical test based on a random null model, but is time-consuming
# jackStraw procedure is very slow, so skip for large projects (>10,000 cells)
calculate_variance = function(seurat_obj, jackstraw_max_cells = 10000) {

  s_obj = seurat_obj

  message("\n\n ========== Seurat::FindVariableGenes() ========== \n\n")

  # detection of variable genes across the single cells
  # FindVariableGenes:
  # - calculates the average expression and dispersion (a normalized measure of cell-to-cell variation) for each gene
  # - places these genes into bins
  # - calculates a z-score for dispersion within each bin
  # default settings: x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 1, y.high.cutoff = Inf
  # pbmc3k tutorial: x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5
  png("variance.meanvar.png", res = 200, width = 8, height = 5, units = "in")
    s_obj = FindVariableGenes(s_obj, mean.function = ExpMean, dispersion.function = LogVMR,
                              x.low.cutoff = 0.05, x.high.cutoff = 5, y.cutoff = 0.5, display.progress = FALSE)
  dev.off()
  message(glue("variable genes: {length(s_obj@var.genes)}"))

  # log to file
  write(glue("variable genes: {length(s_obj@var.genes)}"), file = "create.log", append = TRUE)

  message("\n\n ========== Seurat::ScaleData() ========== \n\n")

  # regress out unwanted sources of variation
  # regressing uninteresting sources of variation can improve dimensionality reduction and clustering
  # could include technical noise, batch effects, biological sources of variation (cell cycle stage)
  # scaled z-scored residuals of these models are stored in scale.data slot
  # used for dimensionality reduction and clustering
  # RegressOut function has been deprecated, and replaced with the vars.to.regress argument in ScaleData
  s_obj = ScaleData(s_obj, vars.to.regress = c("nUMI", "percent.mito"),
                    do.par = TRUE, num.cores = 4, check.for.norm = FALSE, display.progress = FALSE)

  message("\n\n ========== Seurat::PCA() ========== \n\n")

  # use fewer PCs for small datasets
  num_pcs = 50
  if (ncol(s_obj@data) < 100) num_pcs = 20
  if (ncol(s_obj@data) < 25) num_pcs = 5

  # PCA on the scaled data
  # PCA calculation stored in object@dr$pca
  s_obj = RunPCA(s_obj, pc.genes = s_obj@var.genes, pcs.compute = num_pcs,
                 do.print = FALSE, pcs.print = 3, genes.print = 5)

  # ProjectPCA scores each gene in the dataset based on their correlation with the calculated components
  # it can be used to identify markers that are strongly correlated with cellular heterogeneity
  s_obj = ProjectPCA(s_obj, pcs.store = num_pcs, do.print = FALSE)

  # Examine and visualize both cells and genes that define the PCA:
  # - PrintPCA(s_obj, pcs.print = 1:5, genes.print = 5, use.full = TRUE)
  # - VizPCA(s_obj, 1:2)
  # - PCAPlot
  # - PCHeatmap

  # plot the output of PCA analysis (shuffle cells so any one group does not appear overrepresented due to ordering)
  png("variance.pca.png", res = 200, width = 8, height = 6, units = "in")
    PCAPlot(s_obj, 1, 2, pt.size = 1.5, cells.use = sample(colnames(s_obj@data)), cols.use = colors_samples)
  dev.off()

  message("\n\n ========== Seurat::PCHeatmap() ========== \n\n")

  # PCHeatmap() allows for easy exploration of the primary sources of heterogeneity in a dataset
  png("variance.pc.heatmap.png", res = 200, width = 10, height = 14, units = "in")
    PCHeatmap(s_obj, pc.use = 1:15, cells.use = min(ncol(s_obj@data), 500), num.genes = 20,
              do.balanced = TRUE, label.columns = FALSE, use.full = FALSE, do.return = FALSE)
  dev.off()

  message("\n\n ========== Seurat::PCElbowPlot() ========== \n\n")

  # a more ad hoc method for determining PCs to use, draw cutoff where there is a clear elbow in the graph
  plot_elbow = PCElbowPlot(s_obj, num.pc = num_pcs)
  ggsave("variance.pc.elbow.png", plot = plot_elbow, width = 8, height = 5, units = "in")

  # resampling test inspired by the jackStraw procedure - very slow, so skip for large projects (>10,000 cells)
  if (ncol(s_obj@data) < jackstraw_max_cells) {

    message("\n\n ========== Seurat::JackStraw() ========== \n\n")

    # identify significant PCs as those who have a strong enrichment of low p-value genes
    s_obj = JackStraw(s_obj, num.pc = num_pcs, do.par = TRUE, num.cores = 4)

    # compare the distribution of p-values for each PC with a uniform distribution (dashed line)
    # significant PCs will show a strong enrichment of genes with low p-values (solid curve above the dashed line)
    plot_jackstraw = JackStrawPlot(s_obj, PCs = 1:num_pcs, nCol = 10)
    ggsave("variance.pc.jackstraw.png", plot = plot_jackstraw, width = 15, height = 5, units = "in")

  }

  return(s_obj)

}

# determine point size for tSNE plots (smaller for larger datasets)
get_tsne_point_size = function(seurat_obj) {

  pt_size = 1.8
  if (ncol(seurat_obj@data) > 1000) pt_size = 1.2
  if (ncol(seurat_obj@data) > 5000) pt_size = 1.0
  if (ncol(seurat_obj@data) > 10000) pt_size = 0.8
  if (ncol(seurat_obj@data) > 25000) pt_size = 0.6

  return(pt_size)

}

# perform graph-based clustering and tSNE
calculate_clusters = function(seurat_obj, num_dim, reduction_type = "pca") {

  # check if number of dimensions seems reasonable
  if (num_dim < 5) stop("too few dims: ", num_dim)
  if (num_dim > 50) stop("too many dims: ", num_dim)

  s_obj = seurat_obj

  message("\n\n ========== Seurat::FindClusters() ========== \n\n")

  message("reduction type: ", reduction_type)
  message("num dims: ", num_dim)

  message("initial meta.data fields: ", paste(colnames(s_obj@meta.data), collapse = ", "))

  # resolutions for graph-based clustering
  # increased resolution values lead to more clusters (recommendation: 0.6-1.2 for 3K cells, 2-4 for 33K cells)
  res_range = seq(0.1, 2.5, 0.1)
  if (ncol(s_obj@data) > 1000) res_range = c(res_range, 3, 4, 5, 6, 7, 8)

  # algorithm: 1 = original Louvain; 2 = Louvain with multilevel refinement; 3 = SLM
  # save the SNN so that the algorithm can be rerun using the same graph, but with different resolutions
  # FindClusters() uses "reduction.type" and RunTSNE() uses "reduction.use", so that may change in the future
  s_obj = FindClusters(s_obj, reduction.type = reduction_type, dims.use = 1:num_dim,
                       algorithm = 3, resolution = res_range,
                       save.SNN = TRUE, force.recalc = TRUE, print.output = FALSE)

  message("new meta.data fields: ", paste(colnames(s_obj@meta.data), collapse = ", "))

  # PrintFindClustersParams to print a nicely formatted formatted summary of the parameters that were chosen

  message("\n\n ========== Seurat::RunTSNE() ========== \n\n")

  # use tSNE as a tool to visualize, not for clustering directly on tSNE components
  # cells within the graph-based clusters determined above should co-localize on the tSNE plot
  s_obj = RunTSNE(s_obj, reduction.use = reduction_type, dims.use = 1:num_dim, do.fast = TRUE)

  # reduce point size for larger datasets
  tsne_pt_size = get_tsne_point_size(s_obj)

  # tSNE using original sample names (shuffle cells so any one group does not appear overrepresented due to ordering)
  s_obj = SetAllIdent(s_obj, id = "orig.ident")
  plot_tsne = TSNEPlot(s_obj, cells.use = sample(colnames(s_obj@data)), pt.size = tsne_pt_size,
                       colors.use = colors_samples, do.return = TRUE)
  plot_tsne = plot_tsne + theme(aspect.ratio = 1)
  ggsave(glue("tsne.{reduction_type}.{num_dim}.sample.png"), plot = plot_tsne, width = 8, height = 6, units = "in")
  Sys.sleep(1)
  ggsave(glue("tsne.{reduction_type}.{num_dim}.sample.pdf"), plot = plot_tsne, width = 8, height = 6, units = "in")
  Sys.sleep(1)

  # create a separate sub-directory for cluster resolution plots
  clusters_dir = "clusters-resolutions"
  if (!dir.exists(clusters_dir)) dir.create(clusters_dir)

  # for calculated cluster resolutions: remove redundant (same number of clusters), rename, and plot
  res_cols = grep("^res.", colnames(s_obj@meta.data), value = TRUE)
  res_cols = sort(res_cols)
  res_num_clusters_prev = 1
  for (res in res_cols) {

    # proceed if current resolution has more clusters than previous and less than 80 (limited by the color scheme)
    res_num_clusters_cur = s_obj@meta.data[, res] %>% unique() %>% length()
    if (res_num_clusters_cur > res_num_clusters_prev && res_num_clusters_cur < 80) {

      # check if the resolution still has original labels (characters starting with 0)
      if (min(s_obj@meta.data[, res]) == "0") {

        # relabel identities so they start with 1 and not 0
        s_obj@meta.data[, res] = as.numeric(s_obj@meta.data[, res]) + 1
        # pad with 0s to avoid sorting issues
        s_obj@meta.data[, res] = str_pad(s_obj@meta.data[, res], width = 2, side = "left", pad = "0")
        # pad with "C" to avoid downstream numeric conversions
        s_obj@meta.data[, res] = str_c("C", s_obj@meta.data[, res])

      }

      # resolution value based on resolution column name
      res_val = sub("res\\.", "", res)

      # plot file name
      res_str = gsub("\\.", "", res)
      num_clusters = s_obj@meta.data %>% pull(res) %>% unique() %>% length()
      filename = glue("{clusters_dir}/tsne.{reduction_type}.{num_dim}.{res_str}.clust{num_clusters}")

      s_obj = plot_clusters(seurat_obj = s_obj, resolution = res_val, filename_base = filename)

      # add blank line to make output easier to read
      message(" ")

    } else {

      # remove resolution if the number of clusters is same as previous
      s_obj@meta.data = s_obj@meta.data %>% select(-one_of(res))

    }

    # update resolution cluster count for next iteration
    res_num_clusters_prev = res_num_clusters_cur

  }

  message("updated meta.data fields: ", paste(colnames(s_obj@meta.data), collapse = ", "))

  # compile all cell metadata into a single table
  metadata_tbl = s_obj@meta.data %>%
    rownames_to_column("cell") %>% as_tibble() %>%
    mutate(sample_name = orig.ident)
  tsne_tbl = s_obj@dr$tsne@cell.embeddings %>%
    round(3) %>% as.data.frame() %>% rownames_to_column("cell") %>% as_tibble()
  cells_metadata = full_join(metadata_tbl, tsne_tbl, by = "cell") %>% arrange(cell)
  write_excel_csv(cells_metadata, path = "metadata.csv")

  return(s_obj)

}

# plot tSNE with color-coded clusters at specified resolution
plot_clusters = function(seurat_obj, resolution, filename_base) {

  s_obj = seurat_obj

  # set identities based on specified resolution
  s_obj = set_identity(seurat_obj = s_obj, grouping_column = resolution)

  # print stats
  num_clusters = s_obj@ident %>% unique() %>% as.character() %>% length()
  message("resolution: ", resolution)
  message("num clusters: ", num_clusters)

  # generate plot if there is a reasonable number of clusters
  if (num_clusters > 1 && num_clusters < 80) {

    # shuffle cells so they appear randomly and one group does not show up on top
    plot_tsne = TSNEPlot(s_obj, cells.use = sample(colnames(s_obj@data)), pt.size = get_tsne_point_size(s_obj),
                         colors.use = colors_clusters, do.return = TRUE)
    plot_tsne = plot_tsne + theme(aspect.ratio = 1)
    ggsave(glue("{filename_base}.png"), plot = plot_tsne, width = 8, height = 6, units = "in")
    Sys.sleep(1)
    ggsave(glue("{filename_base}.pdf"), plot = plot_tsne, width = 8, height = 6, units = "in")
    Sys.sleep(1)

    if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")

  }

  message(" ")

  return(s_obj)

}

# check grouping variable/resolution against existing meta data columns
check_grouping_column = function(seurat_obj, grouping_column) {

  s_obj = seurat_obj

  # check if the grouping variable is one of meta data columns
  if (!(grouping_column %in% colnames(s_obj@meta.data))) {

    # check if grouping variable is the resolution value (X.X instead of res.X.X)
    res_column = paste0("res.", grouping_column)
    if (res_column %in% colnames(s_obj@meta.data)) {
      grouping_column = res_column
    } else {
      stop("unknown grouping variable: ", grouping_column)
    }

  }

  return(grouping_column)

}

# set identity based on a specified variable/resolution
set_identity = function(seurat_obj, grouping_column) {

  s_obj = seurat_obj

  grouping_column = check_grouping_column(seurat_obj = s_obj, grouping_column = grouping_column)

  # set identities based on selected grouping variable
  message("setting grouping variable: ", grouping_column)
  s_obj = SetAllIdent(s_obj, id = grouping_column)

  return(s_obj)

}

# plot a set of genes
plot_genes = function(seurat_obj, genes, name) {

  # color gradient for FeaturePlot-based plots
  featplot_colors = c("gray85", "red2")

  # tSNE plots color-coded by expression level (should be square to match the original tSNE plots)
  # do.return returns a list, so can't use ggsave
  png(glue("{name}.tsne.png"), res = 200, width = 10, height = 8, units = "in")
    FeaturePlot(seurat_obj, genes, pt.size = 0.5, cols.use = featplot_colors, no.axes = TRUE, no.legend = TRUE)
  dev.off()
  pdf(file = glue("{name}.tsne.pdf"), width = 10, height = 8)
    FeaturePlot(seurat_obj, genes, pt.size = 0.5, cols.use = featplot_colors, no.axes = TRUE, no.legend = TRUE)
  dev.off()

  # dot plot visualization
  dotplot_plot = DotPlot(seurat_obj, genes.plot = genes,
                         dot.scale = 12, plot.legend = TRUE,
                         cols.use = featplot_colors, do.return = TRUE)
  ggsave(glue("{name}.dotplot.png"), plot = dotplot_plot, width = 20, height = 5, units = "in")
  ggsave(glue("{name}.dotplot.pdf"), plot = dotplot_plot, width = 20, height = 5, units = "in")

  # gene violin plots (size.use below 0.2 doesn't seem to make a difference)
  # skip PDF since every cell has to be plotted and they become too big
  vln_plot = VlnPlot(seurat_obj, features.plot = genes,
                     point.size.use = 0.05, size.title.use = 12, size.x.use = 12, size.y.use = 10,
                     single.legend = FALSE, cols.use = colors_clusters, do.return = TRUE)
  ggsave(glue("{name}.violin.png"), plot = vln_plot, width = 15, height = 8, units = "in")

  # expression levels per cluster for bar plots (averaging and output are in non-log space)
  cluster_avg_exp = AverageExpression(seurat_obj, genes.use = genes, show.progress = FALSE)
  cluster_avg_exp_long = cluster_avg_exp %>% rownames_to_column("gene") %>% gather(cluster, avg_exp, -gene)

  # bar plots
  # create a named color scheme to ensure names and colors are in the proper order
  clust_names = levels(seurat_obj@ident)
  color_scheme_named = colors_clusters[1:length(clust_names)]
  names(color_scheme_named) = clust_names
  barplot_plot = ggplot(cluster_avg_exp_long, aes(x = cluster, y = avg_exp, fill = cluster)) +
    geom_col(color = "black") +
    theme(legend.position = "none") +
    scale_fill_manual(values = color_scheme_named) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_cowplot() +
    facet_wrap(~ gene, ncol = 3, scales = "free")
  ggsave(glue("{name}.barplot.png"), plot = barplot_plot, width = 15, height = 8, units = "in")
  ggsave(glue("{name}.barplot.pdf"), plot = barplot_plot, width = 15, height = 8, units = "in")

}

# calculate cluster stats (number of cells, average expression, cell-gene matrix)
calculate_cluster_stats = function(seurat_obj, label) {

  message("\n\n ========== calculate cluster stats ========== \n\n")

  message("cluster names: ", paste(levels(seurat_obj@ident), collapse = ", "))

  # compile relevant cell metadata into a single table
  seurat_obj = StashIdent(object = seurat_obj, save.name = "cluster")
  metadata_tbl = seurat_obj@meta.data %>%
    rownames_to_column("cell") %>% as_tibble() %>%
    select(cell, nGene, nUMI, percent.mito, orig.ident, cluster) %>%
    rename(sample_name = orig.ident)
  tsne_tbl = seurat_obj@dr$tsne@cell.embeddings %>%
    round(3) %>% as.data.frame() %>% rownames_to_column("cell") %>% as_tibble()
  cells_metadata = full_join(metadata_tbl, tsne_tbl, by = "cell") %>% arrange(cell)
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
  cluster_avg_exp = AverageExpression(seurat_obj, show.progress = FALSE)
  cluster_avg_exp = cluster_avg_exp %>% round(3) %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
  write_excel_csv(cluster_avg_exp, path = glue("expression.mean.{label}.csv"))

  Sys.sleep(1)

  # export results split by sample if multiple samples are present
  num_samples = cells_metadata %>% pull(sample_name) %>% n_distinct()
  if (num_samples > 1) {

    # number of cells split by cluster and by sample
    write_excel_csv(summary_cluster_sample, path = glue("summary.{label}.per-sample.csv"))

    # cluster averages split by sample
    sample_avg_exp = AverageExpression(seurat_obj, add.ident = "orig.ident", show.progress = FALSE)
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
  clusters = seurat_obj@ident %>% as.character() %>% unique() %>% sort()

  # use only clusters with more than 10 cells
  clusters = clusters[table(seurat_obj@ident) > 10]

  if (!pairwise) {

    # standard cluster markers calculation

    markers_dir = "markers-global"

    # capture output to avoid excessive warnings
    markers_log = capture.output({
      all_markers = FindAllMarkers(seurat_obj, logfc.threshold = log(1.2), min.pct = 0.20, min.diff.pct = -Inf,
                                   test.use = test, only.pos = FALSE, print.bar = FALSE)
    }, type = "message")

    # do some light filtering and clean up (ROC test returns slightly different output)
    if (test == "roc") {

      all_markers =
        all_markers %>%
        select(cluster, gene, avg_logFC, myAUC, power) %>%
        filter(power > 0.4) %>%
        mutate(avg_logFC = round(avg_logFC, 3), myAUC = round(myAUC, 3), power = round(power, 3)) %>%
        arrange(cluster, -power)
      top_markers = all_markers %>% filter(avg_logFC > 0) %>% group_by(cluster) %>% top_n(20, power)

    } else {

      all_markers =
        all_markers %>%
        select(cluster, gene, avg_logFC, p_val, p_val_adj) %>%
        filter(p_val_adj < 0.001) %>%
        mutate(avg_logFC = round(avg_logFC, 3)) %>%
        arrange(cluster, p_val_adj, p_val)
      top_markers = all_markers %>% filter(avg_logFC > 0) %>% group_by(cluster) %>% top_n(20, avg_logFC)

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
        markers_log = capture.output({
          cur_markers = FindMarkers(seurat_obj, ident.1 = cluster1, ident.2 = cluster2, test.use = test,
                                    logfc.threshold = log(1.1), min.pct = 0.1,
                                    only.pos = TRUE, print.bar = FALSE)
        }, type = "message")

        # clean up markers table (would need to be modified for "roc" test)
        cur_markers =
          cur_markers %>%
          rownames_to_column("gene") %>%
          mutate(cluster = cluster1) %>%
          mutate(cluster2 = cluster2) %>%
          filter(p_val_adj < 0.05) %>%
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

    top_markers = all_markers %>% group_by(cluster) %>% top_n(20, avg_logFC_min)

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

# perform CCA for Seurat "alignment"
# takes original working dir (before moving to analysis dir) and a vector of batch names and batch dirs
# need to run "create" first so that input objects can be manually modified after creation and before alignment
calculate_cca = function(original_wd, batch_analysis_dirs) {

  num_batches = length(unique(batch_analysis_dirs))
  if (num_batches < 2) stop("minimum 2 batches needed for CCA")

  message("\n\n ========== load batch seurat objects ========== \n\n")

  seurat_obj_list = list()
  var_genes_list = list()
  for (i in 1:length(batch_analysis_dirs)) {

    batch_name = glue("B{i}")
    batch_analysis_dir = batch_analysis_dirs[i]
    batch_analysis_dir = glue("{original_wd}/{batch_analysis_dir}")
    batch_seurat_rds = glue("{batch_analysis_dir}/seurat_obj.rds")

    message("loading analysis for batch: ", batch_analysis_dir)

    # check if analysis dir is valid
    if (!dir.exists(batch_analysis_dir)) stop(glue("dir {batch_analysis_dir} does not exist"))
    # check if seurat object exists
    if (!file.exists(batch_seurat_rds)) stop(glue("seurat object rds {batch_seurat_rds} does not exist"))

    s_obj = readRDS(batch_seurat_rds)

    batch_samples = unique(s_obj@meta.data$orig.ident)
    if (length(batch_samples) == 1) { batch_name = glue("{batch_name}:{batch_samples}") }

    # save orig.ident in meta.data (RunMultiCCA sets orig.ident to cell name)
    s_obj@meta.data$batch.orig.ident = s_obj@meta.data$orig.ident

    # save batch name in meta.data
    s_obj@meta.data$batch = batch_name
    seurat_obj_list[[batch_name]] = s_obj

    write(glue("batch: {batch_name}"), file = "create.log", append = TRUE)
    write(glue("sample dir: {basename(batch_analysis_dir)}"), file = "create.log", append = TRUE)
    write(glue("sample cells: {ncol(s_obj@data)}"), file = "create.log", append = TRUE)

    # top 1,000 genes with the highest dispersion (var/mean) from both datasets
    # "rownames in hvg.info are sorted by the variance/mean ratio, which usually works quite well for UMI data"
    # "we chose to take the top 2k genes as it was easy for each dataset to contribute the same number of genes"
    # https://github.com/satijalab/seurat/issues/227
    batch_var_genes = s_obj@hvg.info %>% head(1000) %>% rownames()
    var_genes_list[[batch_name]] = batch_var_genes

  }

  # euler plot of variable gene overlaps (becomes unreadable can take days for many overlaps)
  if (length(var_genes_list) < 10) {
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

  # convert list of variable genes per batch into a single vector
  var_genes = var_genes_list %>% unlist(use.names = FALSE) %>% as.character()

  # check how well the variable genes overlap between datasets
  message(glue("variable genes (1+ batches): {length(which(table(var_genes) > 0))}"))
  write(glue("variable genes (1+ batches): {length(which(table(var_genes) > 0))}"), file = "create.log", append = TRUE)
  message(glue("variable genes (2+ batches): {length(which(table(var_genes) > 1))}"))
  write(glue("variable genes (2+ batches): {length(which(table(var_genes) > 1))}"), file = "create.log", append = TRUE)
  message(glue("variable genes (3+ batches): {length(which(table(var_genes) > 2))}"))
  write(glue("variable genes (3+ batches): {length(which(table(var_genes) > 2))}"), file = "create.log", append = TRUE)

  # select genes to use for CCA: highly variable in more than one datasets
  if (num_batches > 2) {
    var_genes = names(which(table(var_genes) > 1))
  } else {
    var_genes = union(var_genes_list[[1]], var_genes_list[[2]])
  }

  # select genes to use for CCA: present in all datasets
  for (i in 1:length(seurat_obj_list)) {
    var_genes = var_genes[var_genes %in% rownames(seurat_obj_list[[i]]@scale.data)]
  }

  message(glue("variable genes for CCA: {length(var_genes)}"))
  write(glue("variable genes for CCA: {length(var_genes)}"), file = "create.log", append = TRUE)

  message("seurat objects: ", length(seurat_obj_list))

  message("\n\n ========== Seurat::RunMultiCCA() ========== \n\n")

  # run multi-set CCA (must give at least 3 objects/matrices)
  if (num_batches > 2) {
    s_obj = RunMultiCCA(seurat_obj_list, genes.use = var_genes, num.ccs = 40)
  } else {
    s_obj = RunCCA(object = seurat_obj_list[[1]], object2 = seurat_obj_list[[2]],
                   group.by = "batch", genes.use = var_genes, num.cc = 40)
  }

  write(glue("combined cells: {ncol(s_obj@data)}"), file = "create.log", append = TRUE)

  # restore orig.ident in meta.data (RunMultiCCA sets orig.ident to cell name)
  s_obj@meta.data$orig.ident = s_obj@meta.data$batch.orig.ident
  s_obj = SetAllIdent(s_obj, id = "orig.ident")

  # get rid of old clusters
  s_obj@meta.data = s_obj@meta.data %>% select(-starts_with("res"))

  # save raw counts matrix
  counts_raw = s_obj@raw.data %>% as.matrix() %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
  write_csv(counts_raw, path = "counts.raw.csv.gz")

  dist_plot = VlnPlot(s_obj, c("nGene", "nUMI", "percent.mito"), nCol = 3, group.by = "batch",
                      point.size.use = 0.1, x.lab.rot = TRUE, size.title.use = 12, size.x.use = 12,
                      do.sort = TRUE, remove.legend = TRUE, cols.use = colors_samples)
  ggsave("qc.distribution.batch.png", plot = dist_plot, width = 15, height = 6, units = "in")

  message("\n\n ========== Seurat::MetageneBicorPlot() ========== \n\n")

  # MetageneBicorPlot plots correlation strength for each CC and should saturates after a reasonable number of CCs
  bicor_plot = MetageneBicorPlot(s_obj, grouping.var = "batch", dims.eval = 1:40)
  ggsave("variance.cca.bicor.png", plot = bicor_plot, width = 8, height = 5, units = "in")

  # explore the CC dimensions as heatmap
  png("variance.cca.heatmap.png", res = 200, width = 8, height = 16, units = "in")
    DimHeatmap(object = s_obj, reduction.type = "cca", cells.use = 500, dim.use = 1:30, do.balanced = TRUE)
  dev.off()

  # plot pre-alignment CCA-based tSNE
  dim_ident = DimPlot(object = s_obj, cells.use = sample(colnames(s_obj@data)), reduction.use = "cca",
                      group.by = "batch", pt.size = 0.8, cols.use = colors_samples, do.return = TRUE)
  ggsave("variance.cca.tsne.ident.png", plot = dim_ident, width = 7, height = 6, units = "in")

  return(s_obj)

}

# align the CCA subspaces, which creates a new dimensional reduction "cca.aligned"
align_cca = function(seurat_obj, num_ccs) {

  # check if number of CCs seems reasonable
  if (num_ccs < 5) stop("too few CCs: ", num_ccs)
  if (num_ccs > 50) stop("too many CCs: ", num_ccs)

  s_obj = seurat_obj

  message("\n\n ========== Seurat::CalcVarExpRatio() ========== \n\n")

  # search for cells whose expression profile cannot be well-explained by CCA compared to PCA
  # CalcVarExpRatio: calculate the ratio of variance explained by ICA or PCA to CCA
  s_obj = CalcVarExpRatio(object = s_obj, reduction.type = "pca", grouping.var = "batch", dims.use = 1:num_ccs)

  # discard cells where the variance explained by CCA is <2-fold (ratio < 0.5) compared to PCA
  s_obj_unfiltered = s_obj
  s_obj = SubsetData(object = s_obj_unfiltered, subset.name = "var.ratio.pca", accept.low = 0.5)
  s_obj_discarded = SubsetData(object = s_obj_unfiltered, subset.name = "var.ratio.pca", accept.high = 0.5)

  message(glue("unfiltered cells: {nrow(s_obj_unfiltered@meta.data)}"))

  message(glue("discarded cells: {nrow(s_obj_discarded@meta.data)}"))
  write(glue("discarded cells: {nrow(s_obj_discarded@meta.data)}"), file = "create.log", append = TRUE)
  message(glue("discarded median num genes: {median(s_obj_discarded@meta.data$nGene)}"))
  write(glue("discarded median num genes: {median(s_obj_discarded@meta.data$nGene)}"), file = "create.log", append = TRUE)

  message(glue("passing cells: {nrow(s_obj@meta.data)}"))
  write(glue("passing cells: {nrow(s_obj@meta.data)}"), file = "create.log", append = TRUE)
  message(glue("passing mean num genes: {mean(s_obj@meta.data$nGene)}"))
  write(glue("passing mean num genes: {mean(s_obj@meta.data$nGene)}"), file = "create.log", append = TRUE)
  message(glue("passing median num genes: {median(s_obj@meta.data$nGene)}"))
  write(glue("passing median num genes: {median(s_obj@meta.data$nGene)}"), file = "create.log", append = TRUE)

  message("\n\n ========== Seurat::AlignSubspace() ========== \n\n")

  # align the CCA subspaces, which returns a new dimensional reduction called cca.aligned (slow part)
  s_obj = AlignSubspace(object = s_obj, reduction.type = "cca", grouping.var = "batch", dims.align = 1:num_ccs)
  s_obj = SetAllIdent(s_obj, id = "orig.ident")

  # check dimensional reductions
  # str(s_cca@dr)

  # save normalized and log-transformed counts matrix as a basic gzipped text file
  counts_norm = s_obj@data %>% as.matrix() %>% round(3)
  counts_norm = counts_norm %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
  write_csv(counts_norm, path = "counts.normalized.csv.gz")

  message("\n\n ========== Seurat::RunPCA() ========== \n\n")

  # compute PCA to compare with CCA alignment
  # find variable genes using the usual method, not the hvg.info genes from each sample used for CCA
  s_fvg_obj = FindVariableGenes(s_obj, mean.function = ExpMean, dispersion.function = LogVMR,
                                x.low.cutoff = 0.05, x.high.cutoff = 3, y.cutoff = 0.5, display.progress = FALSE)
  s_obj = RunPCA(object = s_obj, pc.genes = s_fvg_obj@var.genes, pcs.compute = 50,
                 do.print = FALSE, pcs.print = 3, genes.print = 5)

  # tSNE based on PCA with different number of dimensions
  for (num_pcs in c(10, 20, 30)) {

    s_obj = RunTSNE(object = s_obj, reduction.use = "pca", dims.use = 1:num_pcs, do.fast = TRUE)
    plot_pca_tsne = TSNEPlot(s_obj, cells.use = sample(colnames(s_obj@data)), group.by = "batch",
                             pt.size = get_tsne_point_size(s_obj), colors.use = colors_samples, do.return = TRUE)
    plot_pca_tsne = plot_pca_tsne + theme(aspect.ratio = 1)
    ggsave(glue("tsne.pca.{num_pcs}.batch.png"), plot = plot_pca_tsne, width = 8, height = 6, units = "in")
    Sys.sleep(1)
    ggsave(glue("tsne.pca.{num_pcs}.batch.pdf"), plot = plot_pca_tsne, width = 8, height = 6, units = "in")
    Sys.sleep(1)

    # plot original samples/libraries if each batch does not correspond to a single sample/library
    if (length(unique(s_obj@meta.data$orig.ident)) > length(unique(s_obj@meta.data$batch))) {
      plot_pca_tsne = TSNEPlot(s_obj, cells.use = sample(colnames(s_obj@data)), group.by = "orig.ident",
                               pt.size = get_tsne_point_size(s_obj), colors.use = colors_samples, do.return = TRUE)
      plot_pca_tsne = plot_pca_tsne + theme(aspect.ratio = 1)
      ggsave(glue("tsne.pca.{num_pcs}.sample.png"), plot = plot_pca_tsne, width = 8, height = 6, units = "in")
      Sys.sleep(1)
      ggsave(glue("tsne.pca.{num_pcs}.sample.pdf"), plot = plot_pca_tsne, width = 8, height = 6, units = "in")
      Sys.sleep(1)
    }

  }

  # tSNE based on cca.aligned
  s_obj = RunTSNE(object = s_obj, reduction.use = "cca.aligned", dims.use = 1:num_ccs, do.fast = TRUE)
  plot_cca_tsne = TSNEPlot(s_obj, cells.use = sample(colnames(s_obj@data)), group.by = "batch",
                           pt.size = get_tsne_point_size(s_obj), colors.use = colors_samples, do.return = TRUE)
  plot_cca_tsne = plot_cca_tsne + theme(aspect.ratio = 1)
  ggsave(glue("tsne.cca.aligned.{num_ccs}.batch.png"), plot = plot_cca_tsne, width = 8, height = 6, units = "in")
  Sys.sleep(1)
  ggsave(glue("tsne.cca.aligned.{num_ccs}.batch.pdf"), plot = plot_cca_tsne, width = 8, height = 6, units = "in")
  Sys.sleep(1)

  # plot original samples/libraries if each batch does not correspond to a single sample/library
  if (length(unique(s_obj@meta.data$orig.ident)) > length(unique(s_obj@meta.data$batch))) {
    plot_cca_tsne = TSNEPlot(s_obj, cells.use = sample(colnames(s_obj@data)), group.by = "orig.ident",
                             pt.size = get_tsne_point_size(s_obj), colors.use = colors_samples, do.return = TRUE)
    plot_cca_tsne = plot_cca_tsne + theme(aspect.ratio = 1)
    ggsave(glue("tsne.cca.aligned.{num_ccs}.sample.png"), plot = plot_cca_tsne, width = 8, height = 6, units = "in")
    Sys.sleep(1)
    ggsave(glue("tsne.cca.aligned.{num_ccs}.sample.pdf"), plot = plot_cca_tsne, width = 8, height = 6, units = "in")
    Sys.sleep(1)
  }

  return(s_obj)

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

# global settings
colors_samples = c(brewer.pal(5, "Set1"), brewer.pal(8, "Dark2"), pal_igv("default")(51))
colors_clusters = c(pal_d3("category10")(10), pal_d3("category20b")(20), pal_igv("default")(51))

# analysis info
analysis_step = "unknown"
out_dir = opts$analysis_dir

# create analysis directory if starting new analysis or exit if analysis already exists
if (opts$create || opts$combine || opts$cca) {

  if (opts$create) analysis_step = "create"
  if (opts$combine) analysis_step = "combine"
  if (opts$cca) analysis_step = "cca"
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

  # create new seurat object based on input sample names and sample directories
  # can work with multiple samples, but the appropriate way is to use "combine" with objects that are pre-filtered
  seurat_obj = load_sample_counts_matrix(opts$sample_name, opts$sample_dir)

  # save right after loading data in case analysis needs to be customized and later steps are disrupted
  saveRDS(seurat_obj, file = "seurat_obj.rds")

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

} else if (opts$cca) {

  # run CCA based on input batch names and previous analysis directories
  # manually specifying batch names allows multiple samples/libraries per batch
  seurat_obj = calculate_cca(original_wd, opts$batch_analysis_dir)
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

    # consider using BuildClusterTree() to merge clusters (see pbmc33k examples)

  }

  if (opts$identify || opts$de) {

    # set resolution in the seurat object
    grouping_column = check_grouping_column(seurat_obj = seurat_obj, grouping_column = opts$resolution)
    seurat_obj = set_identity(seurat_obj = seurat_obj, grouping_column = grouping_column)

    # use a grouping-specific sub-directory for all output
    grouping_label = gsub("\\.", "", grouping_column)
    num_clusters = seurat_obj@ident %>% as.character() %>% unique() %>% length()
    clust_label = glue("clust{num_clusters}")
    res_dir = glue("clusters-{grouping_label}-{clust_label}")
    if (!dir.exists(res_dir)) dir.create(res_dir)
    setwd(res_dir)

    if (opts$identify) {

      analysis_step = "identify"
      message(glue("\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"))

      # create tSNE plot (should already exist in the main directory)
      tsne_filename = glue("tsne.{grouping_label}.{clust_label}")
      seurat_obj = plot_clusters(seurat_obj, resolution = opts$resolution, filename_base = tsne_filename)

      # cluster stat tables (number of cells and average expression)
      calculate_cluster_stats(seurat_obj, label = clust_label)

      # calculate and plot standard cluster markers
      calculate_cluster_markers(seurat_obj, label = clust_label, test = "roc")
      calculate_cluster_markers(seurat_obj, label = clust_label, test = "wilcox")
      calculate_cluster_markers(seurat_obj, label = clust_label, test = "MAST")

      # calculate and plot pairwise cluster markers (very slow, so skip for high number of clusters)
      num_clusters = seurat_obj@ident %>% n_distinct()
      if (num_clusters < 20) {
        calculate_cluster_markers(seurat_obj, label = clust_label, test = "wilcox", pairwise = TRUE)
        calculate_cluster_markers(seurat_obj, label = clust_label, test = "MAST", pairwise = TRUE)
      }

    }

    # cluster stat tables with sample/library info included
    num_samples = seurat_obj@meta.data %>% pull(orig.ident) %>% n_distinct()
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

  if (opts$align) {

    analysis_step = "align"
    message(glue("\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"))

    # save pre-aligned seurat object file since alignment will remove some cells (where CCA and PCA disagree)
    file.rename(from = "seurat_obj.rds", to = "seurat_obj.unfiltered.rds")

    # perform alignment and determine clusters
    seurat_obj = align_cca(seurat_obj, num_ccs = as.integer(opts$num_ccs))
    seurat_obj = calculate_clusters(seurat_obj, num_dim = as.integer(opts$num_ccs), reduction_type = "cca.aligned")
    saveRDS(seurat_obj, file = "seurat_obj.rds")

  }

}

message(glue("\n\n ========== finished analysis step {analysis_step} for {out_dir} ========== \n\n"))

# delete Rplots.pdf
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")



# end
