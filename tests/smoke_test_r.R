#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

root <- "cell_state_de"
if (!dir.exists(root)) {
  root <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "cell_state_de/tests/smoke_test_r.R"), ".."), mustWork = FALSE)
}

tmp <- tempfile("cell_state_de_r_smoke_")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

counts <- Matrix(
  c(
    10, 0, 8, 0, 2, 0, 2, 0,
    5, 0, 4, 0, 1, 0, 1, 0,
    1, 7, 1, 7, 1, 2, 1, 2
  ),
  nrow = 3,
  byrow = TRUE,
  sparse = TRUE
)
rownames(counts) <- c("G1", "G2", "G3")
colnames(counts) <- paste0("c", 1:8)
obj <- CreateSeuratObject(counts = counts)
obj$donor_id <- rep(c("d1", "d2", "d3", "d4"), each = 2)
obj$disease_group <- rep(c("T2D", "T2D", "ND", "ND"), each = 2)
obj$cell_type <- "Beta"
obj$treatments <- "no_treatment"
obj <- NormalizeData(obj, verbose = FALSE)

rds <- file.path(tmp, "toy.rds")
gmt <- file.path(tmp, "toy.gmt")
scores <- file.path(tmp, "scores.tsv.gz")
wide <- file.path(tmp, "scores_wide.tsv.gz")
thresholds <- file.path(tmp, "thresholds.tsv.gz")
metadata <- file.path(tmp, "metadata.tsv.gz")
de <- file.path(tmp, "de.tsv.gz")
membership <- file.path(tmp, "membership.tsv")
saveRDS(obj, rds)
writeLines("pancreas_beta_cell_state_a\ttoy\tG1\tG2", gmt)
write.table(
  data.frame(cell_id = colnames(counts), state = "pancreas_beta_cell_state_a", in_state = TRUE),
  membership,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

score_script <- file.path(root, "scripts", "score_gmt_states_from_seurat.R")
de_script <- file.path(root, "scripts", "pseudobulk_de_from_seurat.R")

score_status <- system2(
  "/opt/homebrew/bin/Rscript",
  c(
    "--vanilla", score_script,
    "--rds", rds,
    "--gmt", gmt,
    "--scores-out", scores,
    "--wide-out", wide,
    "--thresholds-out", thresholds,
    "--null-n", "5",
    "--null-max-cells", "4",
    "--random-seed", "1",
    "--metadata-out", metadata,
    "--state-regex", "^pancreas_beta_cell_"
  ),
  stdout = TRUE,
  stderr = TRUE
)
if (!is.null(attr(score_status, "status"))) {
  stop("score_gmt_states_from_seurat.R failed", call. = FALSE)
}
stopifnot(file.exists(scores), file.exists(wide), file.exists(metadata), file.exists(thresholds))
score_rows <- read.delim(gzfile(scores), check.names = FALSE)
stopifnot(nrow(score_rows) == 8)
stopifnot("ucell_score" %in% colnames(score_rows))
stopifnot("marker_coverage_fraction" %in% colnames(score_rows))
stopifnot(all(score_rows$score_method == "ucell"))
stopifnot(score_rows$ucell_score[score_rows$cell_id == "c1"] > score_rows$ucell_score[score_rows$cell_id == "c2"])
threshold_rows <- read.delim(gzfile(thresholds), check.names = FALSE)
stopifnot(nrow(threshold_rows) == 1)
stopifnot(threshold_rows$n_cells_used == 4)
stopifnot(!is.na(threshold_rows$threshold_value))

de_status <- system2(
  "/opt/homebrew/bin/Rscript",
  c(
    "--vanilla", de_script,
    "--rds", rds,
    "--membership", membership,
    "--analysis-types", "cell_type,state,state_association",
    "--out", de,
    "--donor-col", "donor_id",
    "--group-col", "disease_group",
    "--cell-type-col", "cell_type",
    "--case-values", "T2D",
    "--control-values", "ND",
    "--min-cells", "1",
    "--min-donors", "1",
    "--min-count", "1"
  ),
  stdout = TRUE,
  stderr = TRUE
)
if (!is.null(attr(de_status, "status"))) {
  stop("pseudobulk_de_from_seurat.R failed", call. = FALSE)
}
stopifnot(file.exists(de))
de_rows <- read.delim(gzfile(de), check.names = FALSE)
stopifnot(nrow(de_rows) > 0)
cat("R smoke tests OK\n")
