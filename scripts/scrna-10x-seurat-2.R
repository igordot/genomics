#!/usr/bin/env Rscript


"
Analysis of 10x Genomics Chromium single cell RNA-seq data using Seurat (version 2.1) starting with Cell Ranger output.

Main steps:
  1 - create - import counts matrix, perform initial QC, and calculate various variance metrics (slowest step)
  2 - cluster - perform clustering based on number of PCs
  3 - identify - identify clusters based on resolution (increased resolution results in more clusters)

Extra steps:
  diff - differential expression between samples/libraries within clusters

Usage:
  scrna-10x-seurat.R create <analysis_dir> (<sample_name> <sample_dir>)... [--min_genes=<n> --max_genes=<n> --mt=<n>]
  scrna-10x-seurat.R cluster <analysis_dir> <num_pcs>
  scrna-10x-seurat.R identify <analysis_dir> <resolution>
  scrna-10x-seurat.R diff <analysis_dir> <resolution>
  scrna-10x-seurat.R --help

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

  suppressPackageStartupMessages(library(magrittr))
  suppressPackageStartupMessages(library(glue))
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(Matrix))
  suppressPackageStartupMessages(library(RColorBrewer))
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(cowplot))

}

# convert a sparse matrix of counts to a Seurat object and generate some QC plots
create_seurat_obj = function(counts_matrix, proj_name = NULL, sample_dir = NULL) {

  message("\n\n ========== save raw counts matrix ========== \n\n")

  # remove genes with very few counts
  counts_matrix = counts_matrix[Matrix::rowSums(counts_matrix) >= 5, ]

  message("input genes: ", nrow(counts_matrix))
  message("input cells: ", ncol(counts_matrix))
  message(" ")

  # log to file
  write(glue("input genes: {nrow(counts_matrix)}"), file = "create.log", append = TRUE)
  write(glue("input cells: {ncol(counts_matrix)}"), file = "create.log", append = TRUE)

  # save counts matrix as a standard text file and gzip
  counts_matrix_filename = "counts.raw.txt"
  write.table(as.matrix(counts_matrix), file = counts_matrix_filename, quote = FALSE, sep = "\t", col.names = NA)
  system(paste0("gzip ", counts_matrix_filename))

  message("\n\n ========== create seurat object ========== \n\n")

  # keep all genes expressed in >= 3 cells, keep all cells with >= 200 genes
  # CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, project = "10X_PBMC")

  if (is.null(proj_name)) {

    # if name is not set, then it's a manually merged counts matrix
    # save.raw parameter causes errors in v2
    s_obj = CreateSeuratObject(raw.data = counts_matrix, min.cells = 10, min.genes = 300, project = "proj",
                               names.field = 1, names.delim = ":")

  } else if (proj_name == "aggregated") {

    # multiple libraries combined using Cell Ranger (cellranger aggr)

    # setup taking into consideration aggregated names delimiter
    s_obj = CreateSeuratObject(raw.data = counts_matrix, min.cells = 10, min.genes = 300, project = proj_name,
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

  message(glue("imported genes: {nrow(s_obj@data)}"))
  message(glue("imported cells: {ncol(s_obj@data)}"))
  message(" ")

  # log to file
  write(glue("imported genes: {nrow(s_obj@data)}"), file = "create.log", append = TRUE)
  write(glue("imported cells: {ncol(s_obj@data)}"), file = "create.log", append = TRUE)

  # nGene and nUMI are automatically calculated for every object by Seurat
  # calculate the percentage of mitochondrial genes here and store it in percent.mito using the AddMetaData
  mt_genes = grep("^MT-", rownames(s_obj@data), ignore.case = TRUE, value = TRUE)
  percent_mt = Matrix::colSums(s_obj@raw.data[mt_genes, ]) / Matrix::colSums(s_obj@raw.data)
  percent_mt = round(percent_mt * 100, digits = 3)

  # add columns to object@meta.data, and is a great place to stash QC stats
  s_obj = AddMetaData(s_obj, metadata = percent_mt, col.name = "percent.mito")

  message("\n\n ========== nGene/nUMI/percent_mito plots ========== \n\n")

  dist_unfilt_plot = VlnPlot(s_obj, c("nGene", "nUMI", "percent.mito"), nCol = 3, group.by = "orig.ident",
                             point.size.use = 0.2, x.lab.rot = TRUE, size.title.use = 12, cols.use = color_scheme)
  cowplot::ggsave("qc.distribution.unfiltered.png", plot = dist_unfilt_plot, width = 10, height = 5, units = "in")

  # check for high mitochondrial percentage or low UMI content
  png("qc.correlations.unfiltered.png", res = 200, width = 10, height = 5, units = "in")
    par(mfrow = c(1, 2))
    GenePlot(s_obj, gene1 = "nUMI", gene2 = "percent.mito", cex.use = 1.2, col.use = color_scheme)
    GenePlot(s_obj, gene1 = "nUMI", gene2 = "nGene", cex.use = 1.2, col.use = color_scheme)
  dev.off()

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

  return(s_obj)

}

# create a single Seurat object from multiple 10x Cell Ranger outputs
# takes a vector of sample_names and sample_dirs (can be just one)
load_sample_counts_matrix = function(sample_names, sample_dirs) {

  message("\n\n ========== import cell ranger counts matrix ========== \n\n")

  merged_counts_matrix = NULL

  for (i in 1:length(sample_names)) {

    sample_name = sample_names[i]
    sample_dir = sample_dirs[i]

    message("loading counts matrix for sample: ", sample_name)

    # check if sample dir is valid
    if (!dir.exists(sample_dir)) stop("dir ", sample_dir, " does not exist")

    # determine counts matrix directory
    # "filtered_gene_bc_matrices" for single library
    # "filtered_gene_bc_matrices_mex" for aggregated
    data_dir = paste0(sample_dir, "/outs")
    data_dir = list.files(path = data_dir, pattern = "matrix.mtx", full.names = TRUE, recursive = TRUE)
    data_dir = grep("filtered_gene_bc_matrices", data_dir, value = TRUE)[1]
    data_dir = dirname(data_dir)
    if (!dir.exists(data_dir)) stop("dir ", data_dir, " does not exist")

    message("loading counts matrix dir: ", data_dir)

    counts_matrix = Read10X(data_dir)

    message(glue("library {sample_name} genes: {nrow(counts_matrix)}"))
    message(glue("library {sample_name} cells: {ncol(counts_matrix)}"))
    message(" ")

    # log to file
    write(glue("library {sample_name} genes: {nrow(counts_matrix)}"), file = "create.log", append = TRUE)
    write(glue("library {sample_name} cells: {ncol(counts_matrix)}"), file = "create.log", append = TRUE)

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

# filter data by number of genes and mitochondrial percentage
filter_data = function(seurat_obj, min_genes = NULL, max_genes = NULL, max_mt = 10) {

  s_obj = seurat_obj

  message("\n\n ========== filter data matrix ========== \n\n")

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

  # log to file
  write(glue("min genes cutoff: {min_genes}"), file = "create.log", append = TRUE)
  write(glue("max genes cutoff: {max_genes}"), file = "create.log", append = TRUE)
  write(glue("max mitochondrial percentage cutoff: {max_mt}"), file = "create.log", append = TRUE)

  message("imported genes: ", nrow(s_obj@data))
  message("imported cells: ", ncol(s_obj@data))

  s_obj = FilterCells(s_obj,
                      subset.names = c("nGene", "percent.mito"),
                      low.thresholds = c(min_genes, -Inf),
                      high.thresholds = c(max_genes, max_mt))

  message("filtered genes: ", nrow(s_obj@data))
  message("filtered cells: ", ncol(s_obj@data))

  dist_filt_plot = VlnPlot(s_obj, c("nGene", "nUMI", "percent.mito"), nCol = 3, group.by = "orig.ident",
                           point.size.use = 0.2, x.lab.rot = TRUE, size.title.use = 12, cols.use = color_scheme)
  cowplot::ggsave("qc.distribution.filtered.png", plot = dist_filt_plot, width = 10, height = 5, units = "in")

  # after removing unwanted cells from the dataset, normalize the data
  # LogNormalize:
  # - normalizes the gene expression measurements for each cell by the total expression
  # - multiplies this by a scale factor (10,000 by default)
  # - log-transforms the result
  s_obj = NormalizeData(s_obj, normalization.method = "LogNormalize", scale.factor = 10000, display.progress = FALSE)

  # save counts matrix as a standard text file and gzip
  counts_norm = s_obj@data %>% as.matrix() %>% round(digits = 3)
  norm_matrix_filename = "counts.normalized.txt"
  write.table(counts_norm, file = norm_matrix_filename, quote = FALSE, sep = "\t", col.names = NA)
  system(paste0("gzip ", norm_matrix_filename))

  # log to file
  write(glue("filtered genes: {nrow(s_obj@data)}"), file = "create.log", append = TRUE)
  write(glue("filtered cells: {ncol(s_obj@data)}"), file = "create.log", append = TRUE)

  return(s_obj)

}

# calculate various variance metrics
# PC selection approaches:
# - PCHeatmap - more supervised, exploring PCs to determine relevant sources of heterogeneity
# - PCElbowPlot - heuristic that is commonly used and can be calculated instantly
# - JackStrawPlot - implements a statistical test based on a random null model, but is time-consuming
calculate_variance = function(seurat_obj) {

  s_obj = seurat_obj

  message("\n\n ========== Seurat::FindVariableGenes() ========== \n\n")

  # detection of variable genes across the single cells
  # FindVariableGenes:
  # - calculates the average expression and dispersion for each gene
  # - places these genes into bins
  # - calculates a z-score for dispersion within each bin
  # FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR,
  #                   x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  png("variance.meanvar.png", res = 200, width = 8, height = 5, units = "in")
    # more stringent bottom cutoff (still lower than default 0.1)
    s_obj = FindVariableGenes(s_obj, mean.function = ExpMean, dispersion.function = LogVMR,
                              x.low.cutoff = 0.05, x.high.cutoff = 3, y.cutoff = 0.5, display.progress = FALSE)
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
  s_obj = ScaleData(s_obj, vars.to.regress = c("nUMI", "percent.mito"), display.progress = FALSE)

  message("\n\n ========== Seurat::PCA() ========== \n\n")

  # PCA on the scaled data (can use PCAFast for bigger dataset)
  # PCA calculation stored in object@dr$pca
  s_obj = RunPCA(s_obj, pc.genes = s_obj@var.genes, pcs.compute = 50,
                 do.print = FALSE, pcs.print = 3, genes.print = 5)

  # ProjectPCA scores each gene in the dataset based on their correlation with the calculated components
  # it can be used to identify markers that are strongly correlated with cellular heterogeneity
  s_obj = ProjectPCA(s_obj, pcs.store = 50, do.print = FALSE)

  # Examine and visualize both cells and genes that define the PCA:
  # - PrintPCA(s_obj, pcs.print = 1:5, genes.print = 5, use.full = TRUE)
  # - VizPCA(s_obj, 1:2)
  # - PCAPlot
  # - PCHeatmap

  # Graphs the output of a PCA analysis Cells are colored by their identity class
  png("variance.pca.png", res = 200, width = 8, height = 6, units = "in")
    PCAPlot(s_obj, 1, 2, pt.size = 1.5, cols.use = color_scheme)
  dev.off()

  message("\n\n ========== Seurat::PCHeatmap() ========== \n\n")

  # PCHeatmap() allows for easy exploration of the primary sources of heterogeneity in a dataset
  png("variance.pc.heatmap.png", res = 200, width = 10, height = 14, units = "in")
    PCHeatmap(s_obj, pc.use = 1:15, cells.use = 500, num.genes = 20,
              do.balanced = TRUE, label.columns = FALSE, use.full = FALSE, do.return = FALSE)
  dev.off()

  message("\n\n ========== Seurat::PCElbowPlot() ========== \n\n")

  # A more ad hoc method for determining which PCs to use
  # draw your cutoff where there is a clear elbow in the graph
  plot_elbow = PCElbowPlot(s_obj, num.pc = 30)
  cowplot::ggsave("variance.pc.elbow.png", plot = plot_elbow, width = 8, height = 5, units = "in")

  # resampling test inspired by the jackStraw procedure - very slow, so skip for large projects (>10,000 cells)
  if (ncol(s_obj@data) < 10000) {

    message("\n\n ========== Seurat::JackStraw() ========== \n\n")

    # identify significant PCs as those who have a strong enrichment of low p-value genes
    s_obj = JackStraw(s_obj, num.pc = 40, do.print = FALSE)

    # compare the distribution of p-values for each PC with a uniform distribution (dashed line)
    # significant PCs will show a strong enrichment of genes with low p-values (solid curve above the dashed line)
    plot_jackstraw = JackStrawPlot(s_obj, PCs = 1:40, nCol = 8)
    cowplot::ggsave("variance.pc.jackstraw.png", plot = plot_jackstraw, width = 12, height = 8, units = "in")

  }

  return(s_obj)

}

# perform graph-based clustering and tSNE
calculate_clusters = function(seurat_obj, num_pcs) {

  # check if number of PCs seems reasonable
  if (num_pcs < 5) stop("too few PCs: ", num_pcs)
  if (num_pcs > 50) stop("too many PCs: ", num_pcs)

  s_obj = seurat_obj

  message("\n\n ========== Seurat::FindClusters() ========== \n\n")

  message("initial meta.data fields: ", paste(colnames(s_obj@meta.data), collapse = ", "))
  message("initial identities: ", paste(levels(s_obj@ident), collapse = ", "))

  # graph-based clustering approach
  # clusters are saved in the object@ident slot
  # save the SNN so that the SLM algorithm can be rerun using the same graph, but with different resolutions
  # increased resolution values lead to more clusters (0.6-1.2 for 3K cells, 2-4 for 33K cells)
  # try multiple resolutions here (will be saved as s_obj@data.info columns)
  # s_obj = FindClusters(seurat_obj, pc.use = 1:num_pcs, print.output = FALSE, save.SNN = TRUE)
  s_obj = FindClusters(s_obj, reduction.type = "pca", dims.use = 1:num_pcs, resolution = seq(0.1, 2.1, 0.2),
                       save.SNN = TRUE, print.output = FALSE)

  message("new meta.data fields: ", paste(colnames(s_obj@meta.data), collapse = ", "))
  message("new identities: ", paste(levels(s_obj@ident), collapse = ", "))

  # PrintFindClustersParams to print a nicely formatted formatted summary of the parameters that were chosen

  message("\n\n ========== Seurat::RunTSNE() ========== \n\n")

  # use tSNE as a tool to visualize, not for clustering directly on tSNE components
  # cells within the graph-based clusters determined above should co-localize on the tSNE plot
  s_obj = RunTSNE(s_obj, dims.use = 1:num_pcs, do.fast = TRUE)

  # reduce point size for larger datasets
  tsne_pt_size = 1
  if (ncol(s_obj@data) > 10000) tsne_pt_size = 0.6
  if (ncol(s_obj@data) > 15000) tsne_pt_size = 0.4

  # plot tSNE based on original sample names
  s_obj = SetAllIdent(s_obj, id = "orig.ident")
  plot_tsne = TSNEPlot(s_obj, do.return = TRUE, pt.size = tsne_pt_size, colors.use = color_scheme)
  cowplot::ggsave(glue("tsne.pcs{num_pcs}.original.png"), plot = plot_tsne, width = 7, height = 6, units = "in")
  cowplot::ggsave(glue("tsne.pcs{num_pcs}.original.pdf"), plot = plot_tsne, width = 7, height = 6, units = "in")

  # plot clusters for calculated resolutions
  res_cols = grep("^res.", colnames(s_obj@meta.data), value = TRUE)
  for (res in res_cols) {

    # resolution value based on resolution column name
    res_val = sub("res\\.", "", res)

    # plot file name
    res_str = gsub("\\.", "", res)
    num_clusters = s_obj@meta.data %>% pull(res) %>% unique() %>% length()
    filename = glue("tsne.pcs{num_pcs}.{res_str}.clust{num_clusters}")

    s_obj = plot_clusters(seurat_obj = s_obj, resolution = res_val, plot_filename = filename)

    # perform clustering and add blank line to make output easier to read
    message(" ")

  }

  return(s_obj)

}

# plot tSNE with color-coded clusters at specified resolution
plot_clusters = function(seurat_obj, resolution, plot_filename) {

  s_obj = seurat_obj

  # set identities based on specified resolution
  s_obj = set_resolution(seurat_obj = s_obj, resolution = resolution)

  # print stats
  num_clusters = s_obj@ident %>% unique() %>% as.character() %>% length()
  message("resolution: ", resolution)
  message("clusters: ", num_clusters)

  # generate plot if there is a reasonable number of clusters
  if (num_clusters > 2 && num_clusters < 20) {

    # reduce point size for larger datasets
    tsne_pt_size = 1
    if (ncol(s_obj@data) > 10000) tsne_pt_size = 0.6
    if (ncol(s_obj@data) > 15000) tsne_pt_size = 0.4

    plot_tsne = TSNEPlot(s_obj, do.return = TRUE, pt.size = tsne_pt_size, colors.use = color_scheme)
    plot_filename_png = glue("{plot_filename}.png")
    message("save tSNE plot: ", plot_filename_png)
    ggsave(plot_filename_png, plot = plot_tsne, width = 7, height = 6, units = "in")
    plot_filename_pdf = glue("{plot_filename}.pdf")
    message("save tSNE plot: ", plot_filename_pdf)
    ggsave(plot_filename_pdf, plot = plot_tsne, width = 7, height = 6, units = "in")

  }

  message(" ")

  return(s_obj)

}

# set identity based on a specified resolution
set_resolution = function(seurat_obj, resolution) {

  s_obj = seurat_obj
  res_col = paste0("res.", resolution)

  # stop if resolution was not pre-computed
  if (res_col %in% colnames(s_obj@meta.data)) {
    message("setting resolution: ", resolution)
  } else {
    stop("unknown resolution: ", resolution)
  }

  # set identities based on selected resolution
  s_obj = SetAllIdent(s_obj, id = res_col)

  # relabel identities so they start with 1 and not 0
  levels(s_obj@ident) = as.numeric(levels(s_obj@ident)) + 1

  return(s_obj)

}

# plot a set of genes
plot_genes = function(seurat_obj, genes, name) {

  # gradient color scale for FeaturePlot-based plots (default: c("yellow", "red"))
  featplot_colors = c("gray95", "red2")

  # tSNE plots color-coded by expression level (should be square to match the original tSNE plots)
  # do.return returns a list, so can't use ggsave
  png(glue("{name}.tsne.png"), res = 200, width = 10, height = 8, units = "in")
    FeaturePlot(seurat_obj, genes, pt.size = 0.5, cols.use = featplot_colors, no.axes = TRUE, no.legend = TRUE)
  dev.off()
  pdf(glue("{name}.tsne.pdf"), width = 10, height = 8)
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
                     single.legend = FALSE, cols.use = color_scheme, do.return = TRUE)
  ggsave(glue("{name}.violin.png"), plot = vln_plot, width = 15, height = 8, units = "in")

  # expression levels per cluster for bar plots
  cluster_avg_exp = AverageExpression(seurat_obj, genes.use = genes, show.progress = FALSE)
  cluster_avg_exp_long = cluster_avg_exp %>% rownames_to_column("gene") %>% gather(cluster, avg_exp, -gene)

  # bar plots
  # create a named color scheme to ensure names and colors are in the proper order
  clust_names = levels(seurat_obj@ident)
  color_scheme_named = color_scheme[1:length(clust_names)]
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

  # get number of cells per cluster
  num_cells_per_cluster = tibble(cluster = seurat_obj@ident)
  num_cells_per_cluster = num_cells_per_cluster %>% group_by(cluster) %>% summarise(num_cells = n())
  num_cells_per_cluster_filename = glue("summary.{label}.csv")
  write_excel_csv(num_cells_per_cluster, path = num_cells_per_cluster_filename)

  # create a separate sub-directory for expression values
  exp_dir = "expression"
  if (!dir.exists(exp_dir)) dir.create(exp_dir)

  # gene expression for an "average" single cell in each identity class
  # add add.ident if you want to observe cluster averages separated by replicate
  # output is in log-space, but averaging is done in non-log space
  cluster_avg_exp = AverageExpression(seurat_obj, show.progress = FALSE)
  colnames(cluster_avg_exp) = paste0("clust_", colnames(cluster_avg_exp))
  cluster_avg_exp = cluster_avg_exp %>% round(digits = 3) %>% as.data.frame() %>% rownames_to_column("gene")
  cluster_avg_exp_filename = glue("{exp_dir}/exp.{label}.mean.csv")
  write_excel_csv(cluster_avg_exp, path = cluster_avg_exp_filename)

  # create cell-gene matrix for each cluster
  clusters = seurat_obj@ident %>% unique() %>% sort()
  for (clust_name in clusters) {

    message("cluster name: ", clust_name)

    # subset seurat_obj to only one group
    seurat_subset = SubsetData(seurat_obj, ident.use = clust_name)

    # create matrix and remove genes with no counts
    cluster_exp = seurat_subset@data %>% as.matrix() %>% round(digits = 3)
    cluster_exp = cluster_exp[rowSums(cluster_exp) > 0, ] %>% as.data.frame() %>% rownames_to_column("gene")

    message("cluster genes: ", nrow(cluster_exp))
    message("cluster cells: ", ncol(cluster_exp) - 1)

    clust_exp_csv = glue("{exp_dir}/exp.{label}-{clust_name}.cells.csv.gz")
    write_excel_csv(cluster_exp, path = clust_exp_csv)
  }

}

# calculate cluster markers and plot top ones
calculate_cluster_markers = function(seurat_obj, label, test) {

  message("\n\n ========== calculate cluster markers ========== \n\n")

  # create a separate sub-directory for all markers
  markers_dir = "markers"
  if (!dir.exists(markers_dir)) dir.create(markers_dir)

  message("marker test: ", test)
  all_markers = FindAllMarkers(seurat_obj, logfc.threshold = 0.25, min.pct = 0.25, min.diff.pct = -Inf,
                               test.use = test, only.pos = FALSE, print.bar = FALSE)

  # do some light filtering and clean up (different tests return slighly different output)
  if (test == "roc") {
    # ROC test returns the classification power (ranging from 0 - random, to 1 - perfect)
    all_markers = all_markers %>% select(cluster, gene, avg_logFC, myAUC, power) %>%
      filter(myAUC > 0.3 & power > 0.3) %>%
      mutate(avg_logFC = round(avg_logFC, 3), myAUC = round(myAUC, 3), power = round(power, 3)) %>%
      arrange(cluster, -myAUC)
    top_markers = all_markers %>% group_by(cluster) %>% top_n(20, myAUC)
  } else {
    # wilcox: Wilcoxon rank sum test (default in Seurat 2)
    # bimod: likelihood-ratio test for single cell gene expression based on zero-inflated data (default in Seurat 1)
    # tobit: Tobit-test for differential gene expression as in Trapnell et al., Nature Biotech, 2014
    all_markers = all_markers %>% select(cluster, gene, avg_logFC, p_val, p_val_adj) %>%
      filter(p_val_adj < 0.01 & abs(avg_logFC) > 1) %>%
      mutate(avg_logFC = round(avg_logFC, 3)) %>%
      arrange(cluster, p_val_adj, p_val)
    top_markers = all_markers %>% group_by(cluster) %>% top_n(20, abs(avg_logFC))
  }

  all_markers_csv = glue("{markers_dir}/markers.{label}.{test}.all.csv")
  message("all markers: ", all_markers_csv)
  write_excel_csv(all_markers, path = all_markers_csv)

  top_markers_csv = glue("{markers_dir}/markers.{label}.{test}.top.csv")
  message("top markers: ", all_markers_csv)
  write_excel_csv(top_markers, path = top_markers_csv)

  # cluster names
  clusters = seurat_obj@ident %>% unique() %>% as.character() %>% sort()

  # get marker genes for each cluster
  for (clust_name in clusters) {

    Sys.sleep(1)

    # plot top genes if enough were identified
    filename_label = glue("{markers_dir}/markers.{label}-{clust_name}.{test}")
    cluster_markers = all_markers %>% filter(cluster == clust_name)
    if (nrow(cluster_markers) > 9) {

      if (test == "roc") {
        cluster_markers = cluster_markers %>% arrange(-myAUC)
      } else {
        cluster_markers = cluster_markers %>% arrange(p_val)
      }

      top_cluster_markers = cluster_markers %>% head(12) %$% gene
      plot_genes(seurat_obj, genes = top_cluster_markers, name = filename_label)

    }

  }

}

# calculate differentially expressed genes within each cluster
calculate_cluster_de_genes = function(seurat_obj, label) {

  message("\n\n ========== calculate cluster DE genes ========== \n\n")

  # common settings
  num_de_genes = 100
  heatmap_colors = colorRampPalette(brewer.pal(9, "YlOrRd"))(50)

  # cluster names
  clusters = seurat_obj@ident %>% unique() %>% as.character() %>% sort()

  # get DE genes for each cluster
  for (clust_name in clusters) {

    message(glue("calculating DE genes for cluster {clust_name}"))

    # subset to the specific cluster
    clust_obj = SubsetData(seurat_obj, ident.use = clust_name)

    # revert back to original sample/library labels
    clust_obj = SetAllIdent(clust_obj, id = "orig.ident")

    message("cluster cells: ", ncol(clust_obj@data))
    message("cluster groups: ", paste(levels(clust_obj@ident), collapse = ", "))


    # do differential expression if cluster contains multiple groups and multiple cells per group
    if (length(unique(clust_obj@ident)) > 1 && min(table(clust_obj@ident)) > 5) {

      # iterate through sample/library combinations (relevant if more than two)
      group_combinations = combn(levels(clust_obj@ident), m = 2, simplify = TRUE)
      for (combination_num in 1:ncol(group_combinations)) {

        # determine combination
        s1 = group_combinations[1, combination_num]
        s2 = group_combinations[2, combination_num]
        comparison_name = paste0(s1, "-vs-", s2)
        message("comparison: ", comparison_name)

        filename_label = paste0("de.", label, "-", clust_name, ".", comparison_name, ".top", num_de_genes)

        # find markers (differentially expressed genes) for identity classes (likelihood-ratio test)
        de_genes = FindMarkers(clust_obj, ident.1 = s1, ident.2 = s2, test.use = "bimod",
                               logfc.threshold = 0.25, min.pct = 0.25,
                               only.pos = FALSE, print.bar = FALSE)

        # most significant genes
        top_genes = de_genes %>% rownames_to_column("gene") %>% top_n(num_de_genes, -p_val) %>% dplyr::arrange(-avg_diff)

        # png heatmap
        heatmap_png = paste0("heatmap.", filename_label, ".png")
        message("heatmap png: ", heatmap_png)
        png(heatmap_png, res = 200, width = 10, height = 10, units = "in")
        DoHeatmap(clust_obj, genes.use = top_genes$gene, order.by.ident = TRUE, slim.col.label = TRUE,
                  do.scale = TRUE, remove.key = TRUE,
                  col.use = heatmap_colors, cex.col = 0.5, disp.min = -1, disp.max = 2)
        dev.off()
        Sys.sleep(1)

        # pdf heatmap
        heatmap_pdf = paste0("heatmap.", filename_label, ".pdf")
        message("heatmap pdf: ", heatmap_pdf)
        pdf(heatmap_pdf, width = 10, height = 10)
        DoHeatmap(clust_obj, genes.use = top_genes$gene, order.by.ident = TRUE, slim.col.label = TRUE,
                  do.scale = TRUE, remove.key = TRUE,
                  col.use = heatmap_colors, cex.col = 0.5, disp.min = -1, disp.max = 2)
        dev.off()
        Sys.sleep(1)

        # save stats table (sorted by p-value)
        top_genes = top_genes %>% dplyr::arrange(p_val)
        genes_csv = paste0(filename_label, ".genes.csv")
        message("genes: ", genes_csv)
        write_excel_csv(top_genes, path = genes_csv)

      }

    } else {

      message("skip cluster: ", clust_name)

    }

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

# global settings
color_scheme = c(brewer.pal(9, "Set1"), brewer.pal(8, "Accent"), brewer.pal(8, "Dark2"))

# analysis info
analysis_step = "unknown"
out_dir = opts$analysis_dir

# create analysis directory if starting new analysis or exit if analysis already exists
if (opts$create) {

  analysis_step = "create"
  message(glue("\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"))

  if (dir.exists(out_dir)) {
    stop("output dir ", out_dir, " already exists")
  } else {
    dir.create(out_dir)
  }

}

# set analysis directory as working directory
setwd(out_dir)

# check which command was used
if (opts$create) {

  # log to file
  write(glue("analysis: {out_dir}"), file = "create.log", append = TRUE)

  # create new seurat object based on input sample names and sample directories
  seurat_obj = load_sample_counts_matrix(opts$sample_name, opts$sample_dir)

  # save right after loading data in case analysis needs to be customized and later steps are disrupted
  saveRDS(seurat_obj, file = "seurat_obj.rds")

  # filter by number of genes and mitochondrial genes percentage (optional parameters)
  seurat_obj = filter_data(seurat_obj, min_genes = opts$min_genes, max_genes = opts$max_genes, max_mt = opts$mt)

  # calculate various variance metrics
  seurat_obj = calculate_variance(seurat_obj)

  saveRDS(seurat_obj, file = "seurat_obj.rds")

} else {

  # all commands besides "create" start with an existing seurat object
  if (file.exists("seurat_obj.rds")) {

    message("loading seurat_obj")
    seurat_obj = readRDS("seurat_obj.rds")

  } else {

    stop("seurat obj does not already exist (run 'create' step first)")

  }

  if (opts$cluster) {

    analysis_step = "cluster"
    message(glue("\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"))

    # calculate clusters
    seurat_obj = calculate_clusters(seurat_obj, num_pcs = as.integer(opts$num_pcs))
    saveRDS(seurat_obj, file = "seurat_obj.rds")

    # consider using BuildClusterTree() to merge clusters (see pbmc33k examples)

  }

  if (opts$identify || opts$diff) {

    # set resolution in the seurat object
    res_val = as.numeric(opts$resolution)
    seurat_obj = set_resolution(seurat_obj, resolution = res_val)

    # use a resolution-specific sub-directory for all output
    res_str = gsub("\\.", "", opts$resolution)
    num_clusters = seurat_obj@ident %>% unique() %>% as.character() %>% length()
    res_label = glue("clust{num_clusters}")
    res_dir = glue("clusters-res{res_str}-{res_label}")
    if (!dir.exists(res_dir)) dir.create(res_dir)
    setwd(res_dir)

    if (opts$identify) {

      analysis_step = "identify"
      message(glue("\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"))

      # create tSNE plot (should already exist in the main directory)
      tsne_filename = glue("tsne.res{res_str}.{res_label}")
      seurat_obj = plot_clusters(seurat_obj, resolution = res_val, plot_filename = tsne_filename)

      # cluster stat tables (number of cells and average expression)
      calculate_cluster_stats(seurat_obj, label = res_label)

      # plot cluster markers
      calculate_cluster_markers(seurat_obj, label = res_label, test = "roc")
      calculate_cluster_markers(seurat_obj, label = res_label, test = "wilcox")
      calculate_cluster_markers(seurat_obj, label = res_label, test = "bimod")

    }

    # cluster stat tables with sample/library info included
    num_samples = seurat_obj@meta.data %>% pull(orig.ident) %>% unique() %>% as.character() %>% length()
    if (num_samples > 1) {

      if (opts$identify) {

        res_label = glue("{res_label}-per-sample")
        seurat_obj@meta.data$clust.sample = paste0(seurat_obj@ident[rownames(seurat_obj@meta.data)], "-",
                                                   seurat_obj@meta.data$orig.ident)
        seurat_obj = SetAllIdent(seurat_obj, id = "clust.sample")
        calculate_cluster_stats(seurat_obj, label = res_label)

      }

      # differential expression
      if (opts$diff) {

          analysis_step = "diff"
          message(glue("\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"))

          calculate_cluster_de_genes(seurat_obj, label = res_label)

      }

    }

  }

}

message(glue("\n\n ========== finished analysis step {analysis_step} for {out_dir} ========== \n\n"))

# delete Rplots.pdf
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")



# end
