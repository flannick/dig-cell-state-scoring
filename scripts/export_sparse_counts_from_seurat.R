#!/usr/bin/env Rscript

args <- commandArgs(FALSE)
file_arg <- args[grepl("^--file=", args)]
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "cell_state_de/scripts/export_sparse_counts_from_seurat.R"
target <- file.path(dirname(normalizePath(script_path, mustWork = FALSE)), "export_rank_universe_10x_from_seurat.R")
if (!file.exists(target)) {
  stop("Could not locate export_rank_universe_10x_from_seurat.R", call. = FALSE)
}
source(target, local = new.env(parent = globalenv()))
