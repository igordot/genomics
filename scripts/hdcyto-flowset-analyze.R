#!/usr/bin/env Rscript

"
Streamlined analysis of high-dimensional flow or mass cytometry data.
This script converts a flowSet to a SingleCellExperiment and performs basic analysis.
Based on the CATALYST workflow: https://doi.org/10.12688/f1000research.11622.4

The input is a project directory that contains `r-data/flowset.rds`.
The directory can be created with `hdcyto-fcs-flowset.R`, but any flowSet RDS file should work.

Usage:
  hdcyto-flowset-analyze.R <proj_dir>

" -> doc

# output width
options(width = 120)
# print warnings as they occur
options(warn = 1)
# default type for the bitmap devices such as png (should default to "cairo")
options(bitmapType = "cairo")

# retrieve the command-line arguments
library(docopt)
opts <- docopt(doc)

# check inputs
proj_dir <- opts$proj_dir
if (!dir.exists(proj_dir)) {
  stop("dir does not exist: ", proj_dir)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(janitor)
  library(cowplot)
  library(RColorBrewer)
  library(ggsci)
  library(flowCore)
  library(CATALYST)
})

# should consider filtering uninformative markers

flowset_rds <- glue("{proj_dir}/r-data/flowset.rds")
if (!file.exists(flowset_rds)) {
  stop("flowset.rds not found: ", flowset_rds)
}

sce_rds <- glue("{proj_dir}/r-data/sce.rds")

message("loading flowSet: ", flowset_rds)

fs <- readRDS(flowset_rds)

# fs
# str(fs)

# pData(fs)

# fs[[1]]@parameters@data

# antibodies
# colnames(fs)

# dim(exprs(fs[[1]]))

if (n_distinct(fs[[1]]@parameters@data$name) == 1) {
  stop("parameter names are all identical")
}
if (n_distinct(fs[[1]]@parameters@data$desc) == 1) {
  stop("parameter descriptions are all identical")
}

panel <-
  fs[[1]]@parameters@data %>%
  dplyr::filter(!str_detect(name, "TIME|Time")) %>%
  dplyr::filter(!str_detect(desc, "FSC|SSC|Viability|LiveDead")) %>%
  select(fcs_colname = name, desc) %>%
  mutate(antigen = str_remove(desc, " :.*"))

# remove antibodies
# if (max(str_length(bad_abs)) > 0) {
#   panel = panel %>% dplyr::filter(!str_detect(desc, paste(bad_abs, collapse = "|")))
# }
# panel

md <- pData(fs) %>% dplyr::rename(file_name = name)
# md

# remove rownames (cause an error with diffcyt)
rownames(panel) <- NULL
rownames(md) <- NULL

if (!dir.exists(glue("{proj_dir}/qc"))) {
  dir.create(glue("{proj_dir}/qc"))
}

write_csv(panel, glue("{proj_dir}/qc/sce-panel.csv"))
write_csv(md, glue("{proj_dir}/qc/sce-metadata.csv"))

message("generating SingleCellExperiment")

# if data is pre-normalized
# sce = prepData(fs, panel, md, transform = FALSE, FACS = TRUE)

# sce = prepData(fs, panel, md, transform = TRUE, FACS = TRUE)
# sce = prepData(fs, panel, md, transform = TRUE, cofactor = 5, FACS = TRUE)
sce <- prepData(fs, panel, md, transform = TRUE, cofactor = 150, FACS = TRUE)
# sce

# prepData keeps only sample_id, patient_id, and condition
# keep only sample_id and condition
colData(sce) <- colData(sce)[, c("sample_id", "condition")]
head(colData(sce))

md_extra <-
  left_join(
    dplyr::select(as.data.frame(colData(sce)), sample_id),
    dplyr::select(pData(fs), !any_of(c("condition", "name"))),
    by = c("sample_id")
  ) %>%
  select(!sample_id)
if (!all(md_extra$sample_name == colData(sce)$sample_id)) {
  stop("pData and colData mismatch")
}

colData(sce) <- cbind(colData(sce), md_extra)

# assayNames(sce)
# assayNames(sce) = "exprs"
# assayNames(sce)

levels(sce$sample_id)

# n_cells(sce)

# rowData(sce)

colnames(sce) <- make.names(colData(sce)$sample_id, unique = TRUE)

# colData(sce)

# plotCounts(sce, group_by = "sample_id", color_by = "condition") + scale_fill_igv()

# p = plotExprs(sce, color_by = "condition")
# p = plotExprs(sce, color_by = "sample_id")
# p$facet$params$ncol <- 6
# p

# assay(sce)[1:5, 1:3]

# assay(sce, "counts")[1:5, 1:3]

# assay(sce, "exprs")[1:5, 1:3]

message("plotting densities")

# subset to random samples to make plots more readable
samples_subset <- levels(sce$sample_id)
if (length(samples_subset) > 20) {
  set.seed(99)
  samples_subset <- sort(sample(samples_subset, 20))
}
sce_rand <- sce[, sce$sample_id %in% samples_subset]

# subset to random cells to reduce object size
if (ncol(sce_rand) > 100000) {
  set.seed(99)
  sce_rand <- sce_rand[, sample(colnames(sce_rand), 100000)]
}

# density plot for original values
dens_plot <-
  t(assay(sce_rand, "counts")) %>%
  as_tibble(rownames = "cell_id") %>%
  pivot_longer(!cell_id, names_to = "marker") %>%
  group_by(marker) %>%
  # mutate(min_cutoff = quantile(value, 0.01), max_cutoff = quantile(value, 0.99)) %>%
  mutate(zscore = scale(value)) %>%
  ungroup() %>%
  # dplyr::filter(value > min_cutoff, value < max_cutoff) %>%
  dplyr::filter(zscore > -3, zscore < 3) %>%
  left_join(as_tibble(colData(sce_rand), rownames = "cell_id"), by = "cell_id") %>%
  ggplot(aes(x = value, color = sample_id)) +
  geom_density() +
  facet_wrap(vars(marker), scales = "free") +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank()
  ) +
  scale_color_igv()
save_plot(glue("{proj_dir}/qc/exp-marker-density-original.png"), dens_plot, base_width = 12, base_height = 6)
dens_plot

# density plot for normalized values
dens_plot <-
  t(assay(sce_rand, "exprs")) %>%
  as_tibble(rownames = "cell_id") %>%
  pivot_longer(!cell_id, names_to = "marker") %>%
  group_by(marker) %>%
  # mutate(min_cutoff = quantile(value, 0.005), max_cutoff = quantile(value, 0.995)) %>%
  mutate(zscore = scale(value)) %>%
  ungroup() %>%
  # dplyr::filter(value > min_cutoff, value < max_cutoff) %>%
  dplyr::filter(zscore > -3, zscore < 3) %>%
  left_join(as_tibble(colData(sce_rand), rownames = "cell_id"), by = "cell_id") %>%
  ggplot(aes(x = value, color = sample_id)) +
  geom_density() +
  facet_wrap(vars(marker), scales = "free") +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank()
  ) +
  scale_color_igv()
save_plot(glue("{proj_dir}/qc/exp-marker-density-arcsinh.png"), dens_plot, base_width = 12, base_height = 6)
dens_plot

# density plot for original-scaled values
dens_plot <-
  t(scale(assay(sce_rand))) %>%
  as_tibble(rownames = "cell_id") %>%
  pivot_longer(!cell_id, names_to = "marker") %>%
  dplyr::filter(value > -3, value < 3) %>%
  left_join(as_tibble(colData(sce_rand), rownames = "cell_id"), by = "cell_id") %>%
  ggplot(aes(x = value, color = sample_id)) +
  geom_density() +
  facet_wrap(vars(marker), scales = "free") +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank()
  ) +
  scale_color_igv()
save_plot(glue("{proj_dir}/qc/exp-marker-density-scaled.png"), dens_plot, base_width = 12, base_height = 6)
dens_plot

if (all(c("CD4", "CD8") %in% rownames(sce))) {
  cor_plot <-
    t(assay(sce_rand, "exprs")) %>%
    as_tibble(rownames = "cell_id") %>%
    # dplyr::filter(CD3 > 4) %>%
    left_join(as_tibble(colData(sce_rand), rownames = "cell_id"), by = "cell_id") %>%
    ggplot(aes(x = CD4, y = CD8)) +
    # geom_point(size = 0.1, alpha = 0.2) +
    # geom_density_2d(color = "darkred", alpha = 0.8) +
    geom_density_2d_filled(contour_var = "ndensity") +
    facet_wrap(vars(sample_id), scales = "free") +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white"),
      panel.grid.minor = element_blank(),
      aspect.ratio = 1
    ) +
    scale_fill_viridis_d(option = "plasma")
  save_plot(glue("{proj_dir}/qc/exp-marker-cor-CD4-CD8.png"), cor_plot, base_width = 12, base_height = 12)
}

if (all(c("CD3", "CD19") %in% rownames(sce))) {
  cor_plot <-
    t(assay(sce_rand, "exprs")) %>%
    as_tibble(rownames = "cell_id") %>%
    left_join(as_tibble(colData(sce_rand), rownames = "cell_id"), by = "cell_id") %>%
    ggplot(aes(x = CD3, y = CD19)) +
    # geom_point(size = 0.1, alpha = 0.2) +
    # geom_density_2d(color = "darkred", alpha = 0.8) +
    geom_density_2d_filled(contour_var = "ndensity") +
    facet_wrap(vars(sample_id), scales = "free") +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white"),
      panel.grid.minor = element_blank(),
      aspect.ratio = 1
    ) +
    scale_fill_viridis_d(option = "plasma")
  save_plot(glue("{proj_dir}/qc/exp-marker-cor-CD3-CD19.png"), cor_plot, base_width = 12, base_height = 12)
}

if (!dir.exists(glue("{proj_dir}/exp"))) {
  dir.create(glue("{proj_dir}/exp"))
}

pb_plot <- plotPbExprs(sce, features = NULL) + scale_color_igv()
save_plot(glue("{proj_dir}/exp/exp-marker-pseudobulk.png"), pb_plot, base_width = 8, base_height = 6)
pb_plot

# based on CATALYST::.anno_factors
.anno_factors <- function(x, ids, which, type = c("row", "column")) {
  type <- match.arg(type)
  # get non-numeric cell metadata variables
  cd <- colData(x)
  df <- data.frame(cd, check.names = FALSE)
  df <- select_if(df, ~ !is.numeric(.))
  df <- mutate_all(df, ~ droplevels(factor(.x)))

  # store sample matching
  m <- match(ids, df$sample_id)

  # get number of matches per variable
  ns <- split(df, df$sample_id) %>%
    lapply(mutate_all, droplevels) %>%
    lapply(summarize_all, nlevels) %>%
    do.call(what = "rbind")

  # keep only uniquely mapable factors included in 'which'
  keep <- names(which(colMeans(ns) == 1))
  keep <- setdiff(keep, c("sample_id", "cluster_id"))
  if (is.character(which)) {
    keep <- intersect(keep, which)
  }
  if (length(keep) == 0) {
    return(NULL)
  }
  df <- df[m, keep, drop = FALSE]

  # get list of colors for each annotation
  lvls <- lapply(as.list(df), levels)
  nlvls <- vapply(lvls, length, numeric(1))
  pal <- pal_igv("default")(51)
  if (any(nlvls > length(pal))) {
    pal <- colorRampPalette(pal)(max(nlvls))
  }
  names(is) <- is <- colnames(df)
  cols <- lapply(is, function(i) {
    u <- pal[seq_len(nlvls[i])]
    names(u) <- lvls[[i]]
    u
  })

  ComplexHeatmap::HeatmapAnnotation(
    which = type, df = df,
    col = cols, gp = grid::gpar(col = "white")
  )
}

# based on CATALYST::plotExprHeatmap
plotExprHeatmap_ed <- function(x, features = NULL,
                               by = c("sample_id", "cluster_id", "both"), k = "meta20", m = NULL,
                               assay = "exprs", fun = c("median", "mean", "sum"),
                               scale = c("first", "last", "never"), q = 0.01,
                               row_anno = TRUE, col_anno = TRUE,
                               row_clust = TRUE, col_clust = TRUE,
                               row_dend = TRUE, col_dend = TRUE,
                               bars = FALSE, perc = FALSE, bin_anno = FALSE,
                               hm_pal = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
                               k_pal = CATALYST:::.cluster_cols, m_pal = k_pal,
                               distance = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                               linkage = c("average", "ward.D", "single", "complete", "mcquitty", "median", "centroid", "ward.D2")
                               ) {

  # check validity of input arguments
  args <- as.list(environment())
  CATALYST:::.check_args_plotExprHeatmap(args)
  distance <- match.arg(distance)
  linkage <- match.arg(linkage)
  scale <- match.arg(scale)
  fun <- match.arg(fun)
  by <- match.arg(by)

  # subset features of interest
  x <- x[unique(CATALYST:::.get_features(x, features)), ]

  # get specified cluster IDs
  if (by != "sample_id") {
    .check_k(x, k)
    x$cluster_id <- cluster_ids(x, k)
  }
  if (by == "both") {
    by <- c("cluster_id", "sample_id")
  }

  # aggregate to pseudobulks by sample/cluster/both
  # using 'assay' data & 'fun' as summary statistic
  .do_agg <- function() {
    z <- CATALYST:::.agg(x, by, fun, assay)
    if (length(by) == 1) {
      return(z)
    }
    set_rownames(
      do.call("rbind", z),
      levels(x$cluster_id)
    )
  }
  # do 0-1 scaling for each marker trimming
  # lower ('q'%) & upper (1-'q'%) quantiles
  .do_scale <- function() {
    if (scale == "first") {
      z <- assay(x, assay)
      z <- CATALYST:::.scale_exprs(z, 1, q)
      assay(x, assay, FALSE) <- z
      return(x)
    } else {
      CATALYST:::.scale_exprs(z, 1, q)
    }
  }

  # apply one of...
  # - scale & trim then aggregate
  # - aggregate then scale & trim
  # - aggregate only
  z <- switch(scale,
    first = {
      x <- .do_scale()
      .do_agg()
    },
    last = {
      z <- .do_agg()
      .do_scale()
    },
    never = {
      .do_agg()
    }
  )
  if (length(by) == 1) z <- t(z)

  if (scale != "never" && !(assay == "counts" && fun == "sum")) {
    qs <- round(quantile(z, c(0.01, 0.99)) * 5) / 5
    lgd_aes <- list(at = seq(qs[1], qs[2], 0.2))
  } else {
    lgd_aes <- list()
  }
  lgd_aes$title_gp <- grid::gpar(
    fontsize = 10,
    fontface = "bold",
    lineheight = 0.8
  )

  # left-hand side heatmap annotation:
  # non-numeric cell metadata variables
  if (!isFALSE(row_anno)) {
    left_anno <- switch(by[1],
      sample_id = .anno_factors(x, levels(x$sample_id), row_anno, "row"),
      .anno_clusters(x, k, m, k_pal, m_pal)
    )
  } else {
    left_anno <- NULL
  }
  if (!isFALSE(col_anno) && length(by) == 2) {
    top_anno <- .anno_factors(x, levels(x$sample_id), col_anno, "colum")
  } else {
    top_anno <- NULL
  }

  # right-hand side heatmap annotation:
  # labeled barplot of event counts by sample
  if (bars) {
    right_anno <- .anno_counts(x[[by[1]]], perc)
  } else {
    right_anno <- NULL
  }

  # get bin annotation
  if (bin_anno) {
    cell_fun <- function(j, i, x, y, ...) {
      grid.text(
        gp = gpar(fontsize = 8),
        sprintf("%.2f", z[i, j]), x, y
      )
    }
  } else {
    cell_fun <- NULL
  }

  a <- ifelse(assay == "exprs", "expression", assay)
  f <- switch(fun, "median" = "med", fun)
  hm_title <- switch(scale,
    first = sprintf("%s %s\n%s", fun, "scaled", a),
    last = sprintf("%s %s\n%s", "scaled", fun, a),
    never = paste(fun, a, sep = "\n")
  )
  if (length(by) == 2) {
    col_title <- features
  } else if (length(features) == 1 &&
    features %in% c("type", "state")) {
    col_title <- paste0(features, "_markers")
  } else {
    col_title <- ""
  }

  ComplexHeatmap::Heatmap(
    matrix = z,
    name = hm_title,
    col = circlize::colorRamp2(
      seq(min(z), max(z), l = n <- 100),
      colorRampPalette(hm_pal)(n)
    ),
    column_title = col_title,
    column_title_side = ifelse(length(by) == 2, "top", "bottom"),
    cell_fun = cell_fun,
    cluster_rows = row_clust,
    cluster_columns = col_clust,
    show_row_dend = row_dend,
    show_column_dend = col_dend,
    clustering_distance_rows = distance,
    clustering_method_rows = linkage,
    clustering_distance_columns = distance,
    clustering_method_columns = linkage,
    show_row_names = (
      is.null(left_anno) ||
        isTRUE(by == "sample_id")) && !perc,
    row_names_side = ifelse(
      by[1] == "cluster_id" ||
        isFALSE(row_anno) && !row_dend ||
        isFALSE(row_clust),
      "left", "right"
    ),
    top_annotation = top_anno,
    left_annotation = left_anno,
    right_annotation = right_anno,
    rect_gp = grid::gpar(col = "white"),
    heatmap_legend_param = lgd_aes
  )
}

# hm_plot <- plotExprHeatmap(sce, k_pal = pal_igv()(51), m_pal = pal_igv()(51), scale = "last")
hm_plot <- plotExprHeatmap_ed(sce, scale = "last")
png(glue("{proj_dir}/exp/heatmap-samples.png"), width = 9, height = 6, units = "in", res = 300)
hm_plot
dev.off()

# based on CATALYST::plotExprHeatmap
agg_exprs <- function(x, assay = "exprs", fun = c("median", "mean", "sum"),
                      scale = c("first", "last", "never"), q = 0.01) {

  # check validity of input arguments
  scale <- match.arg(scale)
  fun <- match.arg(fun)
  by <- "sample_id"

  # aggregate to pseudobulks by sample/cluster/both
  # using 'assay' data & 'fun' as summary statistic
  .do_agg <- function() {
    z <- CATALYST:::.agg(x, by, fun, assay)
    if (length(by) == 1) {
      return(z)
    }
    set_rownames(
      do.call("rbind", z),
      levels(x$cluster_id)
    )
  }
  # do 0-1 scaling for each marker trimming
  # lower ('q'%) & upper (1-'q'%) quantiles
  .do_scale <- function() {
    if (scale == "first") {
      z <- assay(x, assay)
      z <- CATALYST:::.scale_exprs(z, 1, q)
      assay(x, assay, FALSE) <- z
      return(x)
    } else {
      CATALYST:::.scale_exprs(z, 1, q)
    }
  }

  # apply one of...
  # - scale & trim then aggregate
  # - aggregate then scale & trim
  # - aggregate only
  z <- switch(scale,
    first = {
      x <- .do_scale()
      .do_agg()
    },
    last = {
      z <- .do_agg()
      .do_scale()
    },
    never = {
      .do_agg()
    }
  )
  if (length(by) == 1) z <- t(z)

  if (scale != "never" && !(assay == "counts" && fun == "sum")) {
    qs <- round(quantile(z, c(0.01, 0.99)) * 5) / 5
    lgd_aes <- list(at = seq(qs[1], qs[2], 0.2))
  } else {
    lgd_aes <- list()
  }
  lgd_aes$title_gp <- grid::gpar(
    fontsize = 10,
    fontface = "bold",
    lineheight = 0.8
  )

  as_tibble(t(z), rownames = "marker")
}

message("aggregating expression")

# aggregate assay data is just (no scaling)
agg_exprs(sce, scale = "never") %>% write_csv(glue("{proj_dir}/exp/exp-samples-agg.csv"))
# aggregate assay data first and scale subsequently (range of each marker will be 0-1)
agg_exprs(sce, scale = "last") %>% write_csv(glue("{proj_dir}/exp/exp-samples-agg-scale.csv"))
# scale and trim then aggregate
# agg_exprs(sce, scale = "first") %>% write_csv(glue("{proj_dir}/exp-samples-scale-agg.csv"))

message("running UMAP")

set.seed(99)
sce <- runDR(sce, dr = "UMAP", cells = 10000, n_neighbors = 50, features = NULL)
sce

set.seed(99)
sce_rand <- sce[, sample(colnames(sce))]

# plot phenotypes
for (p in names(colData(sce))) {
  umap_pheno <- plotDR(sce_rand, "UMAP", color_by = p) + theme_classic() + theme(aspect.ratio = 1) + scale_color_igv()
  save_plot(glue("{proj_dir}/dr-umap-pheno-{p}.png"), umap_pheno, base_width = 8, base_height = 6)
}

# plot markers
marker_colors <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
for (m in sort(rownames(sce))) {
  # hcl.colors(10, "reds", rev = TRUE)
  marker_plot <-
    plotDR(sce_rand, "UMAP", color_by = m, a_pal = marker_colors) +
    theme_cowplot() +
    theme(aspect.ratio = 1, strip.background = element_blank())
  save_plot(glue("{proj_dir}/exp/dr-umap-marker-{m}.png"), marker_plot, base_width = 8, base_height = 6)
}

message("generating clusters")

set.seed(99)
sce <- cluster(sce, features = NULL, xdim = 10, ydim = 10, maxK = 20, seed = 99)

# plotExprHeatmap(sce, by = "cluster_id", k_pal = pal_igv()(51), scale = "last")

# plotAbundances(sce, k = "meta20", by = "sample_id")

set.seed(99)
sce_rand <- sce[, sample(colnames(sce))]

umap_meta20 <- plotDR(sce_rand, "UMAP", color_by = "meta20") + theme_classic() + theme(aspect.ratio = 1) + scale_color_igv()
save_plot(glue("{proj_dir}/dr-umap-cluster-meta20.png"), umap_meta20, base_width = 8, base_height = 6)

message("saving SingleCellExperiment: ", sce_rds)
saveRDS(sce, sce_rds)

# summary_pop_tbl = colData(sce) %>% as.data.frame() %>% janitor::tabyl(sample_id, pop)
# write_csv(summary_pop_tbl, (glue("{proj_dir}/pop/summary-pop.csv")))
# summary_pop_tbl

# delete Rplots.pdf
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}

message("analysis complete in: ", proj_dir)
