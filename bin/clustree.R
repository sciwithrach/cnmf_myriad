#!/usr/bin/env Rscript
# Cluster stability across resolutions with clustree.
#   clustree.R cluster_table.csv clustree_series.txt
# cluster_table.csv: cells x clusterings (one column per resolution, e.g. topo_clusters_ms_res0.2)
# clustree_series.txt: two lines, the column prefix and the plot title (both written by bin/cluster_compare.py)
# Writes clustree.png.

suppressPackageStartupMessages({
    library(clustree)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: clustree.R cluster_table.csv clustree_series.txt")

series   <- readLines(args[2])
prefix   <- series[1]
title    <- series[2]
clusters <- read.csv(args[1], header = TRUE, row.names = 1, check.names = FALSE)

p <- clustree(clusters, prefix = prefix) + ggtitle(title)
ggsave("clustree.png", p, width = 10, height = 15)
