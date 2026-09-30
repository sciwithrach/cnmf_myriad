#!/usr/bin/env Rscript
# Cluster stability across resolutions with clustree, one plot for each clustering series.
#   clustree.R cluster_table.csv clustree_series.txt
# cluster_table.csv: cells x clusterings (one column per resolution, e.g. topo_clusters_ms_res0.2)
# clustree_series.txt: one tab-separated line per series: name, column prefix, plot title
#   (both written by bin/cluster_compare.py)
# Writes clustree_<name>.png for each series.

suppressPackageStartupMessages({
    library(clustree)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: clustree.R cluster_table.csv clustree_series.txt")

clusters <- read.csv(args[1], header = TRUE, row.names = 1, check.names = FALSE)
series   <- read.delim(args[2], header = FALSE, col.names = c("name", "prefix", "title"), stringsAsFactors = FALSE)

for (i in seq_len(nrow(series))) {
    p <- clustree(clusters, prefix = series$prefix[i]) + ggtitle(series$title[i])
    ggsave(paste0("clustree_", series$name[i], ".png"), p, width = 10, height = 15)
}
