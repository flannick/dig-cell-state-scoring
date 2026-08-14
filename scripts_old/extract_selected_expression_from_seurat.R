#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

parse_args <- function(argv) {
  args <- list(
    assay = NA_character_,
    slot = "data",
    layer = NA_character_,
    cell_id_col = "cell_id",
    metadata_cols = NA_character_,
    cell_filter_col = NA_character_,
    cell_filter_values = NA_character_,
    drop_zero = FALSE
  )
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!grepl("^--", key)) {
      stop("Unexpected positional argument: ", key, call. = FALSE)
    }
    name <- sub("^--", "", key)
    if (name == "drop-zero") {
      args$drop_zero <- TRUE
      i <- i + 1
      next
    }
    if (i == length(argv)) {
      stop("Missing value for ", key, call. = FALSE)
    }
    value <- argv[[i + 1]]
    name <- gsub("-", "_", name)
    args[[name]] <- value
    i <- i + 2
  }
  args
}

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  extract_selected_expression_from_seurat.R --rds map.rds --genes genes.txt",
      "    --metadata-out metadata.tsv.gz --expression-out expression.tsv.gz",
      "",
      "Optional:",
      "  --assay RNA --slot data --layer data",
      "  --metadata-cols donor_id,cell_type,disease_group",
      "  --cell-filter-col cell_type --cell-filter-values 'Type A,Type B'",
      "  --drop-zero",
      sep = "\n"
    ),
    "\n"
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("rds", "genes", "metadata_out", "expression_out")
missing <- required[!nzchar(unlist(args[required]))]
if (length(missing) > 0) {
  usage()
  stop("Missing required argument(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

read_gene_list <- function(path) {
  genes <- readLines(path, warn = FALSE)
  genes <- trimws(genes)
  unique(genes[nzchar(genes) & !grepl("^#", genes)])
}

split_csv <- function(value) {
  if (is.na(value) || !nzchar(value)) {
    return(character())
  }
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

write_tsv <- function(df, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(df, con, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

message("Reading Seurat object: ", args$rds)
obj <- readRDS(args$rds)
assay <- args$assay
if (is.na(assay) || !nzchar(assay)) {
  assay <- DefaultAssay(obj)
}
source_label <- if (!is.na(args$layer) && nzchar(args$layer)) {
  paste0("layer=", args$layer)
} else {
  paste0("slot=", args$slot)
}
message("Using assay=", assay, " ", source_label)

genes <- read_gene_list(args$genes)
if (length(genes) == 0) {
  stop("No genes found in ", args$genes, call. = FALSE)
}

meta <- obj@meta.data
meta[[args$cell_id_col]] <- rownames(meta)

if (!is.na(args$cell_filter_col) && nzchar(args$cell_filter_col)) {
  if (!args$cell_filter_col %in% names(meta)) {
    stop("Cell filter column not found in metadata: ", args$cell_filter_col, call. = FALSE)
  }
  filter_values <- split_csv(args$cell_filter_values)
  if (length(filter_values) == 0) {
    stop("--cell-filter-values is required when --cell-filter-col is set", call. = FALSE)
  }
  keep <- meta[[args$cell_filter_col]] %in% filter_values
  meta <- meta[keep, , drop = FALSE]
}

cells <- meta[[args$cell_id_col]]
metadata_cols <- split_csv(args$metadata_cols)
if (length(metadata_cols) > 0) {
  metadata_cols <- unique(c(args$cell_id_col, metadata_cols))
  absent <- setdiff(metadata_cols, names(meta))
  if (length(absent) > 0) {
    warning("Metadata columns not found and omitted: ", paste(absent, collapse = ", "))
  }
  metadata_cols <- intersect(metadata_cols, names(meta))
  meta_out <- meta[, metadata_cols, drop = FALSE]
} else {
  meta_out <- meta
}

message("Writing metadata for ", nrow(meta_out), " cells: ", args$metadata_out)
write_tsv(meta_out, args$metadata_out)

get_assay_matrix <- function(obj, assay, slot, layer) {
  layer_or_slot <- slot
  if (!is.na(layer) && nzchar(layer)) {
    layer_or_slot <- layer
  }
  tryCatch(
    GetAssayData(obj, assay = assay, layer = layer_or_slot),
    error = function(e) GetAssayData(obj, assay = assay, slot = layer_or_slot)
  )
}

expr <- get_assay_matrix(obj, assay, args$slot, args$layer)
present_genes <- intersect(genes, rownames(expr))
missing_genes <- setdiff(genes, present_genes)
if (length(missing_genes) > 0) {
  warning("Genes not present in expression matrix: ", paste(missing_genes, collapse = ", "))
}
if (length(present_genes) == 0) {
  stop("None of the requested genes were present in the expression matrix", call. = FALSE)
}

cells <- intersect(cells, colnames(expr))
if (length(cells) == 0) {
  stop("No selected metadata cells were present in the expression matrix", call. = FALSE)
}

message("Writing selected expression for ", length(present_genes), " genes and ", length(cells), " cells: ", args$expression_out)
con <- gzfile(args$expression_out, open = "wt")
writeLines("cell_id\tgene\texpression", con)

for (gene in present_genes) {
  values <- as.numeric(expr[gene, cells, drop = TRUE])
  if (isTRUE(args$drop_zero)) {
    keep <- values != 0
    values <- values[keep]
    out_cells <- cells[keep]
  } else {
    out_cells <- cells
  }
  if (length(out_cells) == 0) {
    next
  }
  write.table(
    data.frame(cell_id = out_cells, gene = gene, expression = values),
    con,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    na = ""
  )
}

close(con)
message("Done")
