#!/usr/bin/env Rscript

'
Description:
  Mark putative doublets in single-cell RNA-seq data stored as a Seurat object using scDblFinder.
  Input is a Seurat object stored as an RDS file.

Usage:
  scrna-doublets-scdblfinder.R <seurat_rds>

Arguments:
  <seurat_rds>     input directory

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
seurat_rds = opts$seurat_rds
out_dir = dirname(seurat_rds)

message("seurat object: ", seurat_rds)
message("output dir: ", out_dir)

# check if the input is valid
if (!file.exists(seurat_rds)) { stop("object file does not exist") }
if (!dir.exists(out_dir)) { stop("output dir does not exist") }

# load libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(glue)
  library(scran)
  library(scDblFinder)
})

# import seurat object
seurat_obj = readRDS(seurat_rds)

# set output directory as working directory
setwd(out_dir)

# check if output exists already
if (file.exists("doublets.scDblFinder.csv.gz")) { stop("output already exists") }

# log to file
write(glue("scDblFinder version: {packageVersion('scDblFinder')}"), file = "create.log", append = TRUE)

Idents(seurat_obj) = "orig.ident"
sce = Seurat::as.SingleCellExperiment(seurat_obj, assay = "RNA")
sce = scran::computeSumFactors(sce, BPPARAM = BiocParallel::MulticoreParam(4))

set.seed(99)
if (all(c("hash.ID", "HTO_classification.global") %in% names(seurat_obj@meta.data))) {
  # hashed multi-sample experiment
  # samples are independent captures, not biological samples, if multiplexed using cell hashes
  if ("library" %in% names(seurat_obj@meta.data)) {
    if ("Doublet" %in% seurat_obj@meta.data$HTO_classification.global) {
      message("library/batch: library with known doublets")
      known_doublets = sce$HTO_classification.global == "Doublet"
    } else {
      message("library/batch: library without known doublets")
      known_doublets = NULL
    }
    doublet_tbl =
      scDblFinder(
        sce, samples = "library", knownDoublets = known_doublets,
        returnType = "table", BPPARAM = BiocParallel::MulticoreParam(4)
      )
  } else {
    stop("hashed multi-sample experiment should have a 'library' metadata column")
  }
} else if (n_distinct(seurat_obj@meta.data$orig.ident) > 1) {
  # multi-sample experiment
  message("library/batch: orig.ident")
  doublet_tbl = scDblFinder(sce, samples = "orig.ident", returnType = "table", BPPARAM = BiocParallel::MulticoreParam(4))
} else {
  message("library/batch: none")
  doublet_tbl = scDblFinder(sce, returnType = "table", BPPARAM = BiocParallel::MulticoreParam(4))
}

# using the "samples" parameter does not return a table (fixed in 1.11.4)
if (class(doublet_tbl) == "SingleCellExperiment") {
  doublet_tbl = colData(doublet_tbl) %>% as.data.frame() %>% dplyr::select(starts_with("scDblFinder"))
  colnames(doublet_tbl) = stringr::str_remove(colnames(doublet_tbl), "scDblFinder.")
  doublet_tbl$type = "real"
}
doublet_tbl = doublet_tbl %>% as_tibble(rownames = "cell") %>% dplyr::filter(type == "real") %>% dplyr::arrange(cell)
write_csv(doublet_tbl, "doublets.scDblFinder.csv.gz")

if (nrow(doublet_tbl) != ncol(seurat_obj)) { stop("doublet table and seurat object are not the same size") }

# add doublet stats to the seurat object
doublet_tbl = doublet_tbl %>% select(cell, doublet_score_scDblFinder = score, doublet_class_scDblFinder = class)
doublet_df = doublet_tbl %>% as.data.frame() %>% column_to_rownames("cell") %>% sample_frac()
seurat_obj = AddMetaData(seurat_obj, doublet_df)
seurat_obj@meta.data$doublet_class_scDblFinder = factor(seurat_obj@meta.data$doublet_class_scDblFinder)

# check doublet rate
num_doublets = table(seurat_obj@meta.data$doublet_class_scDblFinder)[["doublet"]]
doublet_rate = round(num_doublets / ncol(seurat_obj), 3)
message(glue("doublet rate: {doublet_rate}"))
write(glue("num doublets: {num_doublets}"), file = "create.log", append = TRUE)
write(glue("doublet rate: {doublet_rate}"), file = "create.log", append = TRUE)

# dot size for plots
num_cells = ncol(seurat_obj)
pt_size = 1.8
if (num_cells > 1000) pt_size = 1.4
if (num_cells > 5000) pt_size = 1.0
if (num_cells > 10000) pt_size = 0.6
if (num_cells > 50000) pt_size = 0.2

# plot doublet score
featplot_colors = colorRampPalette(c("#d9cfcb", "#d49070", "#ca5528", "#b72600", "#981000", "#730000"))(100)
random_cells = sample(colnames(seurat_obj))
plot_umap =
  FeaturePlot(
    seurat_obj, features = "doublet_score_scDblFinder", reduction = "umap",
    cells = random_cells, pt.size = pt_size, cols = featplot_colors
  ) +
  theme_cowplot() +
  theme(
    plot.background = element_rect(fill = "white"),
    aspect.ratio = 1, plot.title = element_text(hjust = 0.5),
    axis.ticks = element_blank(), axis.text = element_blank()
  )
save_plot("dr.umap.doublet_score_scDblFinder.png", plot = plot_umap, base_height = 6.5, base_width = 8)
Sys.sleep(1)
save_plot("dr.umap.doublet_score_scDblFinder.pdf", plot = plot_umap, base_height = 6.5, base_width = 8)
Sys.sleep(1)

# plot doublet class
plot_umap =
  DimPlot(
    seurat_obj, group.by = "doublet_class_scDblFinder", reduction = "umap",
    cells = random_cells, pt.size = pt_size, cols = c("#E41A1C", "#377EB8")
  ) +
  theme_cowplot() +
  theme(
    plot.background = element_rect(fill = "white"),
    aspect.ratio = 1, plot.title = element_text(hjust = 0.5),
    axis.ticks = element_blank(), axis.text = element_blank()
  )
save_plot("dr.umap.doublet_class_scDblFinder.png", plot = plot_umap, base_height = 6.5, base_width = 8)
Sys.sleep(1)
save_plot("dr.umap.doublet_class_scDblFinder.pdf", plot = plot_umap, base_height = 6.5, base_width = 8)
Sys.sleep(1)

# delete Rplots.pdf
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")

# save
Idents(seurat_obj) = "orig.ident"
saveRDS(seurat_obj, "seurat_obj.rds")



# end
