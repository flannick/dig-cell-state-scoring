#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(Matrix)
  library(Seurat)
})

parse_args <- function(argv) {
  args <- list(
    rds = NA_character_,
    tarball = NA_character_,
    tar_member = NA_character_,
    assay = NA_character_,
    slot = "counts",
    layer = NA_character_,
    membership = NA_character_,
    analysis_types = "cell_type,state,state_association",
    cell_id_col = "cell_id",
    donor_col = "donor_id",
    group_col = "disease_group",
    phenotype_type = "categorical",
    cell_type_col = "cell_type",
    cell_filter_col = NA_character_,
    cell_filter_values = NA_character_,
    treatment_col = NA_character_,
    treatment_values = NA_character_,
    case_values = "T2D",
    control_values = "ND",
    derive_any_aab_col = NA_character_,
    aab_cols = "aab_gada_positive,aab_ia2_positive,aab_znt8_positive,aab_iaa_positive",
    min_cells = "20",
    min_donors = "3",
    min_count = "10"
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

label_values <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  x
}

as_bool <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
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

aggregate_counts <- function(counts, meta, cell_id_col, donor_col, label_col, min_cells) {
  cells <- intersect(meta[[cell_id_col]], colnames(counts))
  meta <- meta[match(cells, meta[[cell_id_col]]), , drop = FALSE]
  meta$sample_id <- paste(meta[[donor_col]], meta[[label_col]], sep = "||")
  n_cells <- table(meta$sample_id)
  keep_samples <- names(n_cells)[n_cells >= min_cells]
  meta <- meta[meta$sample_id %in% keep_samples, , drop = FALSE]
  if (nrow(meta) == 0) {
    return(NULL)
  }
  sample_factor <- factor(meta$sample_id)
  design_mat <- Matrix::sparse.model.matrix(~ 0 + sample_factor)
  colnames(design_mat) <- levels(sample_factor)
  pb <- counts[, meta[[cell_id_col]], drop = FALSE] %*% design_mat
  sample_info <- unique(meta[, c("sample_id", donor_col, label_col), drop = FALSE])
  rownames(sample_info) <- sample_info$sample_id
  sample_info <- sample_info[colnames(pb), , drop = FALSE]
  kept_n_cells <- as.integer(n_cells[colnames(pb)])
  names(kept_n_cells) <- colnames(pb)
  list(counts = pb, samples = sample_info, n_cells = kept_n_cells)
}

aggregate_counts_by_donor <- function(counts, meta, cell_id_col, donor_col, min_cells) {
  cells <- intersect(meta[[cell_id_col]], colnames(counts))
  meta <- meta[match(cells, meta[[cell_id_col]]), , drop = FALSE]
  meta$sample_id <- as.character(meta[[donor_col]])
  n_cells <- table(meta$sample_id)
  keep_samples <- names(n_cells)[n_cells >= min_cells]
  meta <- meta[meta$sample_id %in% keep_samples, , drop = FALSE]
  if (nrow(meta) == 0) {
    return(NULL)
  }
  sample_factor <- factor(meta$sample_id)
  design_mat <- Matrix::sparse.model.matrix(~ 0 + sample_factor)
  colnames(design_mat) <- levels(sample_factor)
  pb <- counts[, meta[[cell_id_col]], drop = FALSE] %*% design_mat
  sample_info <- unique(meta[, c("sample_id", donor_col), drop = FALSE])
  rownames(sample_info) <- sample_info$sample_id
  sample_info <- sample_info[colnames(pb), , drop = FALSE]
  kept_n_cells <- as.integer(n_cells[colnames(pb)])
  names(kept_n_cells) <- colnames(pb)
  list(counts = pb, samples = sample_info, n_cells = kept_n_cells)
}

run_edgeR <- function(pb_counts, sample_info, label_col, donor_col, n_cells, case_values, control_values, min_donors, min_count) {
  sample_info$contrast_group <- ifelse(sample_info[[label_col]] %in% case_values, "case",
    ifelse(sample_info[[label_col]] %in% control_values, "control", NA)
  )
  keep <- !is.na(sample_info$contrast_group)
  sample_info <- sample_info[keep, , drop = FALSE]
  pb_counts <- pb_counts[, rownames(sample_info), drop = FALSE]
  n_cells <- n_cells[colnames(pb_counts)]
  donor_counts <- table(sample_info$contrast_group)
  if (any(!(c("case", "control") %in% names(donor_counts))) ||
      donor_counts[["case"]] < min_donors || donor_counts[["control"]] < min_donors) {
    message("Skipping contrast with donor counts: ", paste(names(donor_counts), donor_counts, sep = "=", collapse = ", "))
    return(NULL)
  }
  y <- DGEList(counts = pb_counts)
  keep_genes <- rowSums(y$counts) >= min_count
  y <- y[keep_genes, , keep.lib.sizes = FALSE]
  if (nrow(y) == 0) {
    message("Skipping contrast because no genes passed min_count=", min_count)
    return(NULL)
  }
  group <- factor(sample_info$contrast_group, levels = c("control", "case"))
  design <- model.matrix(~ group)
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  tab <- topTags(glmQLFTest(fit, coef = 2), n = Inf, sort.by = "none")$table
  data.frame(
    gene = rownames(tab),
    log_fc = tab$logFC,
    pvalue = tab$PValue,
    qvalue = tab$FDR,
    case_donors = length(unique(sample_info[[donor_col]][sample_info$contrast_group == "case"])),
    control_donors = length(unique(sample_info[[donor_col]][sample_info$contrast_group == "control"])),
    case_cells = sum(n_cells[sample_info$contrast_group == "case"]),
    control_cells = sum(n_cells[sample_info$contrast_group == "control"]),
    stringsAsFactors = FALSE
  )
}

run_edgeR_continuous <- function(pb_counts, sample_info, donor_col, n_cells, phenotype, phenotype_name, min_donors, min_count) {
  sample_info$phenotype_value <- as.numeric(phenotype[rownames(sample_info)])
  keep <- is.finite(sample_info$phenotype_value)
  sample_info <- sample_info[keep, , drop = FALSE]
  pb_counts <- pb_counts[, rownames(sample_info), drop = FALSE]
  n_cells <- n_cells[colnames(pb_counts)]
  if (length(unique(sample_info[[donor_col]])) < min_donors) {
    message("Skipping continuous association for ", phenotype_name, " with donor count=", length(unique(sample_info[[donor_col]])))
    return(NULL)
  }
  if (length(unique(sample_info$phenotype_value)) < 3) {
    message("Skipping continuous association for ", phenotype_name, " with fewer than 3 unique values")
    return(NULL)
  }
  y <- DGEList(counts = pb_counts)
  keep_genes <- rowSums(y$counts) >= min_count
  y <- y[keep_genes, , keep.lib.sizes = FALSE]
  if (nrow(y) == 0) {
    message("Skipping continuous association because no genes passed min_count=", min_count)
    return(NULL)
  }
  design <- model.matrix(~ phenotype_value, data = sample_info)
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  tab <- topTags(glmQLFTest(fit, coef = "phenotype_value"), n = Inf, sort.by = "none")$table
  data.frame(
    gene = rownames(tab),
    log_fc = tab$logFC,
    pvalue = tab$PValue,
    qvalue = tab$FDR,
    case_donors = length(unique(sample_info[[donor_col]])),
    control_donors = NA_integer_,
    case_cells = sum(n_cells),
    control_cells = NA_integer_,
    stringsAsFactors = FALSE
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("out")
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

message("Reading input Seurat object")
obj <- read_input_rds(args)
message("Loaded object")
assay <- args$assay
if (is.na(assay) || !nzchar(assay)) {
  assay <- DefaultAssay(obj)
}
message("Extracting assay matrix")
counts <- get_assay_matrix(obj, assay, args$slot, args$layer)
message("Loaded assay matrix with ", nrow(counts), " genes and ", ncol(counts), " cells")
meta <- obj@meta.data
meta[[args$cell_id_col]] <- rownames(meta)

if (!is.na(args$derive_any_aab_col) && nzchar(args$derive_any_aab_col)) {
  aab_cols <- split_csv(args$aab_cols)
  missing_aab_cols <- setdiff(aab_cols, names(meta))
  if (length(missing_aab_cols) > 0) {
    stop("Metadata missing autoantibody column(s): ", paste(missing_aab_cols, collapse = ", "), call. = FALSE)
  }
  aab_frame <- data.frame(lapply(meta[, aab_cols, drop = FALSE], as_bool), check.names = FALSE)
  any_observed <- rowSums(!is.na(meta[, aab_cols, drop = FALSE])) > 0
  meta[[args$derive_any_aab_col]] <- ifelse(!any_observed, NA_character_, ifelse(rowSums(aab_frame) > 0, "positive", "negative"))
}

required_cols <- c(args$cell_id_col, args$donor_col, args$group_col, args$cell_type_col)
missing_cols <- setdiff(required_cols, names(meta))
if (length(missing_cols) > 0) {
  stop("Metadata missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

if (!is.na(args$treatment_col) && nzchar(args$treatment_col)) {
  treatment_values <- split_csv(args$treatment_values)
  meta <- meta[meta[[args$treatment_col]] %in% treatment_values, , drop = FALSE]
}

if (!is.na(args$cell_filter_col) && nzchar(args$cell_filter_col)) {
  if (!args$cell_filter_col %in% names(meta)) {
    stop("Cell filter column not found in metadata: ", args$cell_filter_col, call. = FALSE)
  }
  cell_filter_values <- split_csv(args$cell_filter_values)
  meta <- meta[meta[[args$cell_filter_col]] %in% cell_filter_values, , drop = FALSE]
}

analysis_types <- split_csv(args$analysis_types)
case_values <- split_csv(args$case_values)
control_values <- split_csv(args$control_values)
phenotype_type <- args$phenotype_type
if (!phenotype_type %in% c("categorical", "continuous")) {
  stop("--phenotype-type must be categorical or continuous", call. = FALSE)
}
min_cells <- as.integer(args$min_cells)
min_donors <- as.integer(args$min_donors)
min_count <- as.integer(args$min_count)
all_results <- list()

if ("cell_type" %in% analysis_types) {
  for (ct in sort(unique(meta[[args$cell_type_col]]))) {
    sub <- meta[meta[[args$cell_type_col]] == ct, , drop = FALSE]
    if (phenotype_type == "continuous") {
      phenotype <- tapply(as.numeric(sub[[args$group_col]]), sub[[args$donor_col]], function(x) unique(x[is.finite(x)])[1])
      agg <- aggregate_counts_by_donor(counts, sub, args$cell_id_col, args$donor_col, min_cells)
      if (is.null(agg)) next
      de <- run_edgeR_continuous(agg$counts, agg$samples, args$donor_col, agg$n_cells, phenotype, args$group_col, min_donors, min_count)
    } else {
      sub$de_label <- label_values(sub[[args$group_col]])
      agg <- aggregate_counts(counts, sub, args$cell_id_col, args$donor_col, "de_label", min_cells)
      if (is.null(agg)) next
      de <- run_edgeR(agg$counts, agg$samples, "de_label", args$donor_col, agg$n_cells, case_values, control_values, min_donors, min_count)
    }
    if (is.null(de)) next
    de$analysis_type <- "cell_type"
    de$cell_type <- ct
    de$state <- NA_character_
    de$case <- if (phenotype_type == "continuous") args$group_col else paste(case_values, collapse = ",")
    de$control <- if (phenotype_type == "continuous") "continuous" else paste(control_values, collapse = ",")
    de$test <- if (phenotype_type == "continuous") "edgeR_glmQLF_donor_pseudobulk_continuous" else "edgeR_glmQLF_donor_pseudobulk"
    all_results[[length(all_results) + 1]] <- de
    rm(agg, de)
    gc()
  }
}

if (any(c("state", "state_association") %in% analysis_types)) {
  if (is.na(args$membership) || !nzchar(args$membership)) {
    stop("--membership is required for state analyses", call. = FALSE)
  }
  membership <- read.delim(args$membership, stringsAsFactors = FALSE, check.names = FALSE)
  membership <- membership[as_bool(membership$in_state), , drop = FALSE]
}

if ("state" %in% analysis_types) {
  for (state in sort(unique(membership$state))) {
    cells <- unique(membership$cell_id[membership$state == state])
    sub <- meta[meta[[args$cell_id_col]] %in% cells, , drop = FALSE]
    if (phenotype_type == "continuous") {
      phenotype <- tapply(as.numeric(sub[[args$group_col]]), sub[[args$donor_col]], function(x) unique(x[is.finite(x)])[1])
      agg <- aggregate_counts_by_donor(counts, sub, args$cell_id_col, args$donor_col, min_cells)
      if (is.null(agg)) next
      de <- run_edgeR_continuous(agg$counts, agg$samples, args$donor_col, agg$n_cells, phenotype, args$group_col, min_donors, min_count)
    } else {
      sub$de_label <- label_values(sub[[args$group_col]])
      agg <- aggregate_counts(counts, sub, args$cell_id_col, args$donor_col, "de_label", min_cells)
      if (is.null(agg)) next
      de <- run_edgeR(agg$counts, agg$samples, "de_label", args$donor_col, agg$n_cells, case_values, control_values, min_donors, min_count)
    }
    if (is.null(de)) next
    de$analysis_type <- "state"
    de$cell_type <- if ("cell_type" %in% names(membership)) paste(unique(membership$cell_type[membership$state == state]), collapse = ",") else NA_character_
    de$state <- state
    de$case <- if (phenotype_type == "continuous") args$group_col else paste(case_values, collapse = ",")
    de$control <- if (phenotype_type == "continuous") "continuous" else paste(control_values, collapse = ",")
    de$test <- if (phenotype_type == "continuous") "edgeR_glmQLF_donor_pseudobulk_continuous" else "edgeR_glmQLF_donor_pseudobulk"
    all_results[[length(all_results) + 1]] <- de
    rm(agg, de)
    gc()
  }
}

if ("state_association" %in% analysis_types) {
  for (state in sort(unique(membership$state))) {
    pos_cells <- unique(membership$cell_id[membership$state == state])
    ct_values <- if ("cell_type" %in% names(membership)) unique(membership$cell_type[membership$state == state]) else unique(meta[[args$cell_type_col]])
    sub <- meta[meta[[args$cell_type_col]] %in% ct_values, , drop = FALSE]
    sub$de_label <- ifelse(sub[[args$cell_id_col]] %in% pos_cells, "state_positive", "state_negative")
    agg <- aggregate_counts(counts, sub, args$cell_id_col, args$donor_col, "de_label", min_cells)
    if (is.null(agg)) next
    de <- run_edgeR(agg$counts, agg$samples, "de_label", args$donor_col, agg$n_cells, "state_positive", "state_negative", min_donors, min_count)
    if (is.null(de)) next
    de$analysis_type <- "state_association"
    de$cell_type <- paste(ct_values, collapse = ",")
    de$state <- state
    de$case <- "state_positive"
    de$control <- "state_negative"
    de$test <- "edgeR_glmQLF_donor_pseudobulk"
    all_results[[length(all_results) + 1]] <- de
    rm(agg, de)
    gc()
  }
}

if (length(all_results) == 0) {
  stop("No DE results produced after filters", call. = FALSE)
}
out <- do.call(rbind, all_results)
out <- out[, c("analysis_type", "cell_type", "state", "gene", "log_fc", "pvalue", "qvalue", "case_donors", "control_donors", "case_cells", "control_cells", "case", "control", "test")]
write_tsv(out, args$out)
