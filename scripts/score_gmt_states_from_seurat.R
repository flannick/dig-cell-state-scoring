#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

parse_args <- function(argv) {
  args <- list(
    rds = NA_character_,
    tarball = NA_character_,
    tar_member = NA_character_,
    assay = NA_character_,
    slot = "data",
    layer = NA_character_,
    cell_id_col = "cell_id",
    metadata_cols = NA_character_,
    cell_filter_col = NA_character_,
    cell_filter_values = NA_character_,
    state_regex = NA_character_,
    min_markers_found = "1",
    wide_out = NA_character_,
    score_method = "ucell",
    max_rank = "1500",
    thresholds_out = NA_character_,
    null_n = "500",
    null_percentile = "0.99",
    null_max_cells = "20000",
    random_seed = "1",
    expression_bins = "20",
    detection_bins = "5",
    rank_block_size = "500"
  )
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!grepl("^--", key)) {
      stop("Unexpected positional argument: ", key, call. = FALSE)
    }
    if (i == length(argv)) {
      stop("Missing value for ", key, call. = FALSE)
    }
    name <- gsub("-", "_", sub("^--", "", key))
    args[[name]] <- argv[[i + 1]]
    i <- i + 2
  }
  args
}

split_csv <- function(value) {
  if (is.na(value) || !nzchar(value)) {
    return(character())
  }
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

read_gmt <- function(path, state_regex) {
  lines <- readLines(path, warn = FALSE)
  rows <- lapply(lines, function(line) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 3) {
      return(NULL)
    }
    list(state = parts[[1]], library = parts[[2]], genes = unique(parts[-c(1, 2)]))
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!is.na(state_regex) && nzchar(state_regex)) {
    rows <- rows[vapply(rows, function(x) grepl(state_regex, x$state), logical(1))]
  }
  rows
}

write_tsv <- function(df, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(df, con, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

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

rank_matrix_for_ucell <- function(expr, max_rank, genes_to_keep = rownames(expr), block_size = 500) {
  genes_to_keep <- intersect(unique(genes_to_keep), rownames(expr))
  out <- matrix(NA_real_, nrow = length(genes_to_keep), ncol = ncol(expr))
  rownames(out) <- genes_to_keep
  colnames(out) <- colnames(expr)
  starts <- seq(1, ncol(expr), by = block_size)
  for (start in starts) {
    end <- min(start + block_size - 1, ncol(expr))
    block_cells <- colnames(expr)[start:end]
    dense <- as.matrix(expr[, block_cells, drop = FALSE])
    ranks <- apply(dense, 2, function(x) {
      r <- rank(-x, ties.method = "average", na.last = "keep")
      r[is.na(r) | r > max_rank] <- max_rank + 1
      r
    })
    if (is.null(dim(ranks))) {
      ranks <- matrix(ranks, ncol = 1)
      rownames(ranks) <- rownames(expr)
      colnames(ranks) <- block_cells
    }
    out[, block_cells] <- ranks[genes_to_keep, , drop = FALSE]
  }
  out
}

score_ucell_from_ranks <- function(ranks, genes, max_rank) {
  present <- intersect(genes, rownames(ranks))
  n <- length(present)
  if (n == 0) {
    return(rep(NA_real_, ncol(ranks)))
  }
  max_u <- n * max_rank - (n * (n + 1)) / 2
  if (max_u <= 0) {
    return(rep(NA_real_, ncol(ranks)))
  }
  rank_sum <- Matrix::colSums(ranks[present, , drop = FALSE])
  u_stat <- rank_sum - (n * (n + 1)) / 2
  score <- 1 - (u_stat / max_u)
  pmin(pmax(as.numeric(score), 0), 1)
}

score_mean_expression <- function(expr, genes) {
  present <- intersect(genes, rownames(expr))
  if (length(present) == 0) {
    return(rep(NA_real_, ncol(expr)))
  }
  as.numeric(Matrix::colMeans(expr[present, , drop = FALSE]))
}

score_gene_set <- function(expr, ranks, genes, score_method, max_rank) {
  if (score_method == "ucell") {
    score_ucell_from_ranks(ranks, genes, max_rank)
  } else if (score_method == "mean_expression") {
    score_mean_expression(expr, genes)
  } else {
    stop("Unsupported --score-method: ", score_method, call. = FALSE)
  }
}

percent_rank <- function(x) {
  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }
  rank(x, ties.method = "average", na.last = "keep") / sum(!is.na(x))
}

make_bins <- function(values, n_bins) {
  values <- as.numeric(values)
  if (length(unique(values[is.finite(values)])) <= 1) {
    return(rep(1L, length(values)))
  }
  probs <- seq(0, 1, length.out = n_bins + 1)
  cuts <- unique(as.numeric(quantile(values, probs = probs, na.rm = TRUE, names = FALSE)))
  if (length(cuts) <= 2) {
    return(rep(1L, length(values)))
  }
  as.integer(cut(values, breaks = cuts, include.lowest = TRUE, labels = FALSE))
}

sample_matched_genes <- function(target_genes, universe, bin_key, rng_exclude) {
  sampled <- character()
  for (gene in target_genes) {
    candidates <- universe[bin_key[universe] == bin_key[[gene]] & !(universe %in% rng_exclude)]
    if (length(candidates) == 0) {
      candidates <- universe[!(universe %in% rng_exclude)]
    }
    if (length(candidates) == 0) {
      candidates <- universe
    }
    sampled <- c(sampled, sample(candidates, 1))
  }
  sampled
}

calibrate_null_thresholds <- function(expr, sets, score_method, max_rank, null_n, null_percentile,
                                      null_max_cells, random_seed, expression_bins, detection_bins,
                                      tissue, cell_type, rank_block_size) {
  set.seed(random_seed)
  n_cells <- ncol(expr)
  if (n_cells > null_max_cells) {
    calibration_cells <- sort(sample(colnames(expr), null_max_cells))
  } else {
    calibration_cells <- colnames(expr)
  }
  cal_expr <- expr[, calibration_cells, drop = FALSE]
  mean_expr <- Matrix::rowMeans(cal_expr)
  detected <- Matrix::rowMeans(cal_expr > 0)
  expr_bin <- make_bins(mean_expr, expression_bins)
  detect_bin <- make_bins(detected, detection_bins)
  bin_key <- paste(expr_bin, detect_bin, sep = ":")
  names(bin_key) <- rownames(expr)
  universe <- rownames(expr)

  rows <- list()
  for (i in seq_along(sets)) {
    set <- sets[[i]]
    present <- intersect(set$genes, rownames(expr))
    marker_coverage <- length(present) / length(set$genes)
    if (length(present) == 0) {
      threshold <- NA_real_
    } else {
      random_sets <- vector("list", null_n)
      for (j in seq_len(null_n)) {
        random_sets[[j]] <- sample_matched_genes(present, universe, bin_key, present)
      }
      cal_ranks <- NULL
      if (score_method == "ucell") {
        cal_ranks <- rank_matrix_for_ucell(
          cal_expr,
          max_rank,
          genes_to_keep = unique(c(present, unlist(random_sets, use.names = FALSE))),
          block_size = rank_block_size
        )
      }
      null_scores <- numeric()
      for (j in seq_len(null_n)) {
        null_scores <- c(null_scores, score_gene_set(cal_expr, cal_ranks, random_sets[[j]], score_method, max_rank))
      }
      threshold <- as.numeric(quantile(null_scores, probs = null_percentile, na.rm = TRUE, names = FALSE))
    }
    rows[[i]] <- data.frame(
      tissue = tissue,
      cell_type = cell_type,
      state = set$state,
      threshold_method = paste0("matched_random_gene_set_null", format(null_percentile, trim = TRUE)),
      threshold_value = threshold,
      null_percentile = null_percentile,
      mixture_boundary = NA_real_,
      n_cells_used = length(calibration_cells),
      marker_coverage_fraction = marker_coverage,
      null_n = null_n,
      score_method = score_method,
      max_rank = max_rank,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("gmt", "scores_out", "metadata_out")
missing <- required[!nzchar(unlist(args[required]))]
if (length(missing) > 0) {
  stop("Missing required argument(s): ", paste(missing, collapse = ", "), call. = FALSE)
}
if ((!is.na(args$rds) && nzchar(args$rds)) == (!is.na(args$tarball) && nzchar(args$tarball))) {
  stop("Provide exactly one of --rds or --tarball", call. = FALSE)
}
if (!is.na(args$tarball) && nzchar(args$tarball) && (is.na(args$tar_member) || !nzchar(args$tar_member))) {
  stop("--tar-member is required when --tarball is used", call. = FALSE)
}

read_input_rds <- function(args) {
  if (!is.na(args$rds) && nzchar(args$rds)) {
    return(readRDS(args$rds))
  }
  cmd <- paste("tar -xOzf", shQuote(args$tarball), shQuote(args$tar_member))
  con <- pipe(cmd, open = "rb")
  on.exit(close(con), add = TRUE)
  readRDS(con)
}

obj <- read_input_rds(args)
assay <- args$assay
if (is.na(assay) || !nzchar(assay)) {
  assay <- DefaultAssay(obj)
}
expr <- get_assay_matrix(obj, assay, args$slot, args$layer)
score_method <- args$score_method
if (!score_method %in% c("ucell", "mean_expression")) {
  stop("--score-method must be ucell or mean_expression", call. = FALSE)
}
max_rank <- as.integer(args$max_rank)
if (is.na(max_rank) || max_rank < 2) {
  stop("--max-rank must be an integer >= 2", call. = FALSE)
}

meta <- obj@meta.data
meta[[args$cell_id_col]] <- rownames(meta)
if (!is.na(args$cell_filter_col) && nzchar(args$cell_filter_col)) {
  if (!args$cell_filter_col %in% names(meta)) {
    stop("Cell filter column not found in metadata: ", args$cell_filter_col, call. = FALSE)
  }
  filter_values <- split_csv(args$cell_filter_values)
  meta <- meta[meta[[args$cell_filter_col]] %in% filter_values, , drop = FALSE]
}
cells <- intersect(meta[[args$cell_id_col]], colnames(expr))
if (length(cells) == 0) {
  stop("No selected metadata cells were present in expression matrix", call. = FALSE)
}
meta <- meta[cells, , drop = FALSE]

metadata_cols <- split_csv(args$metadata_cols)
if (length(metadata_cols) > 0) {
  metadata_cols <- unique(c(args$cell_id_col, metadata_cols))
  metadata_cols <- intersect(metadata_cols, names(meta))
  meta_out <- meta[, metadata_cols, drop = FALSE]
} else {
  meta_out <- meta
}
write_tsv(meta_out, args$metadata_out)

sets <- read_gmt(args$gmt, args$state_regex)
min_markers_found <- as.integer(args$min_markers_found)
if (length(sets) == 0) {
  stop("No GMT marker sets matched the requested filter", call. = FALSE)
}

expr <- expr[, cells, drop = FALSE]
ranks <- NULL
rank_block_size <- as.integer(args$rank_block_size)
if (is.na(rank_block_size) || rank_block_size < 1) {
  stop("--rank-block-size must be a positive integer", call. = FALSE)
}
if (score_method == "ucell") {
  score_genes <- unique(unlist(lapply(sets, function(set) set$genes), use.names = FALSE))
  ranks <- rank_matrix_for_ucell(expr, max_rank, genes_to_keep = score_genes, block_size = rank_block_size)
}
long_con <- gzfile(args$scores_out, open = "wt")
writeLines("cell_id\tstate\tscore\tucell_score\tscore_method\tscore_percentile_within_scope\tn_markers\tmarker_genes_total\tn_markers_found\tmarker_genes_present\tmarker_coverage_fraction\tcell_type_scope\tlibrary", long_con)
wide <- data.frame(cell_id = cells, stringsAsFactors = FALSE)

for (set in sets) {
  present <- intersect(set$genes, rownames(expr))
  if (length(present) < min_markers_found) {
    next
  }
  score <- score_gene_set(expr, ranks, present, score_method, max_rank)
  score_percentile <- percent_rank(score)
  ucell_score <- if (score_method == "ucell") score else NA_real_
  marker_coverage <- length(present) / length(set$genes)
  cell_type_scope <- sub("^pancreas_", "", set$state)
  cell_type_scope <- sub("_cell_.*$", " cell", cell_type_scope)
  cell_type_scope <- gsub("_", " ", cell_type_scope)
  out <- data.frame(
    cell_id = cells,
    state = set$state,
    score = as.numeric(score),
    ucell_score = as.numeric(ucell_score),
    score_method = score_method,
    score_percentile_within_scope = as.numeric(score_percentile),
    n_markers = length(set$genes),
    marker_genes_total = length(set$genes),
    n_markers_found = length(present),
    marker_genes_present = length(present),
    marker_coverage_fraction = marker_coverage,
    cell_type_scope = cell_type_scope,
    library = set$library,
    stringsAsFactors = FALSE
  )
  write.table(out, long_con, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, na = "")
  wide[[set$state]] <- as.numeric(score)
}

close(long_con)
if (!is.na(args$wide_out) && nzchar(args$wide_out)) {
  write_tsv(wide, args$wide_out)
}
if (!is.na(args$thresholds_out) && nzchar(args$thresholds_out)) {
  cell_type <- NA_character_
  if (!is.na(args$cell_filter_col) && nzchar(args$cell_filter_col) && length(split_csv(args$cell_filter_values)) == 1) {
    cell_type <- split_csv(args$cell_filter_values)[[1]]
  }
  thresholds <- calibrate_null_thresholds(
    expr = expr,
    sets = sets,
    score_method = score_method,
    max_rank = max_rank,
    null_n = as.integer(args$null_n),
    null_percentile = as.numeric(args$null_percentile),
    null_max_cells = as.integer(args$null_max_cells),
    random_seed = as.integer(args$random_seed),
    expression_bins = as.integer(args$expression_bins),
    detection_bins = as.integer(args$detection_bins),
    tissue = NA_character_,
    cell_type = cell_type,
    rank_block_size = rank_block_size
  )
  write_tsv(thresholds, args$thresholds_out)
}

message("Wrote scores for ", ncol(wide) - 1, " states and ", length(cells), " cells")
