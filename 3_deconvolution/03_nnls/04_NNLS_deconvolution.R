############################################
# 04_NNLS_deconvolution.R
# Cell-type deconvolution using NNLS (non-negative least squares)
#
# Inputs:
#  - Signature matrices (genes x celltypes):
#     3_deconvolution/03_nnls/input_signature_matrices/
#       - SignatureMatrix_Zeisel2015_CIBERSORTx.txt
#       - SignatureMatrix_Tasic2018_CIBERSORTx.txt
#       - SignatureMatrix_Yao2021_CIBERSORTx.txt
#
#  - Bulk matrix (genes x samples):
#     data/bulk/CPMcounts_ComBatSeq.txt
#
# Outputs:
#   results/deconvolution/nnls/
#     - NNLS_Zeisel2015_proportions.csv
#     - NNLS_Tasic2018_proportions.csv
#     - NNLS_Yao2021_proportions.csv
############################################

############################################
## 0. Libraries and paths
############################################

suppressPackageStartupMessages({
  library(nnls)
})

# ---- Inputs ----
sig_dir <- file.path("3_deconvolution", "03_nnls", "input_signature_matrices")

sig_zeisel_file <- file.path(sig_dir, "SignatureMatrix_Zeisel2015_CIBERSORTx.txt")
sig_tasic_file  <- file.path(sig_dir, "SignatureMatrix_Tasic2018_CIBERSORTx.txt")
sig_yao_file    <- file.path(sig_dir, "SignatureMatrix_Yao2021_CIBERSORTx.txt")

bulk_file <- file.path("data", "bulk", "CPMcounts_ComBatSeq.txt")

# ---- Outputs ----
out_dir <- file.path("results", "deconvolution", "nnls")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_zeisel <- file.path(out_dir, "NNLS_Zeisel2015_proportions.csv")
out_tasic  <- file.path(out_dir, "NNLS_Tasic2018_proportions.csv")
out_yao    <- file.path(out_dir, "NNLS_Yao2021_proportions.csv")

# ---- Checks ----
min_common_genes <- 1000

############################################
## 1. Helper functions
############################################

quick_range_check <- function(mat, name) {
  cat("\n--- CHECK:", name, "---\n")
  cat("Dimensions: ", paste(dim(mat), collapse = " x "), "\n")
  cat("Range: ", paste(range(mat, na.rm = TRUE), collapse = " to "), "\n")
}

coerce_matrix_numeric <- function(m, label = "matrix") {
  m <- as.matrix(m)
  suppressWarnings(storage.mode(m) <- "numeric")
  if (anyNA(m)) {
    warning("NA values found after numeric coercion in ", label,
            ". Genes with NA will be removed before NNLS.")
  }
  m
}

read_expression_matrix <- function(path, label) {
  if (!file.exists(path)) stop(label, " not found: ", path)
  m <- read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  m <- coerce_matrix_numeric(m, label = paste0(label, ": ", basename(path)))
  m
}

nnls_proportions <- function(sig_mat, bulk_mat) {
  # sig_mat: genes x celltypes
  # bulk_mat: genes x samples
  stopifnot(identical(rownames(sig_mat), rownames(bulk_mat)))
  
  samples   <- colnames(bulk_mat)
  celltypes <- colnames(sig_mat)
  
  props <- matrix(NA_real_, nrow = ncol(bulk_mat), ncol = ncol(sig_mat),
                  dimnames = list(samples, celltypes))
  
  for (i in seq_along(samples)) {
    y <- bulk_mat[, i]
    fit <- nnls::nnls(as.matrix(sig_mat), y)
    x <- fit$x
    
    s <- sum(x, na.rm = TRUE)
    if (is.na(s) || s == 0) {
      props[i, ] <- 0
    } else {
      props[i, ] <- x / s
    }
  }
  
  props
}

write_props_csv <- function(props_mat, out_path) {
  df <- as.data.frame(props_mat)
  df$sample <- rownames(props_mat)
  df <- df[, c("sample", setdiff(colnames(df), "sample")), drop = FALSE]
  write.csv(df, out_path, row.names = FALSE)
}

run_one_reference <- function(signature_mat, bulk_mat, ref_name) {
  cat("\n========== ", ref_name, " ==========\n", sep = "")
  
  common_genes <- intersect(rownames(bulk_mat), rownames(signature_mat))
  common_genes <- unique(common_genes)
  
  cat("Common genes: ", length(common_genes), "\n", sep = "")
  if (length(common_genes) < min_common_genes) {
    warning(ref_name, ": low number of common genes (", length(common_genes),
            "). Check gene identifiers (symbols vs Ensembl) and preprocessing.")
  }
  
  # IMPORTANT: enforce consistent ordering
  common_genes <- sort(common_genes)
  
  bulk_common <- bulk_mat[common_genes, , drop = FALSE]
  sig_common  <- signature_mat[common_genes, , drop = FALSE]
  
  # Remove genes with NA in either matrix
  keep <- complete.cases(bulk_common) & complete.cases(sig_common)
  if (!all(keep)) {
    cat("Removing genes with NA: ", sum(!keep), "\n", sep = "")
    bulk_common <- bulk_common[keep, , drop = FALSE]
    sig_common  <- sig_common[keep, , drop = FALSE]
  }
  
  stopifnot(identical(rownames(bulk_common), rownames(sig_common)))
  
  props <- nnls_proportions(sig_common, bulk_common)
  props
}

############################################
## 2. Load matrices
############################################

signature_zeisel <- read_expression_matrix(sig_zeisel_file, "Signature Zeisel")
signature_tasic  <- read_expression_matrix(sig_tasic_file,  "Signature Tasic")
signature_yao    <- read_expression_matrix(sig_yao_file,    "Signature Yao")

bulk_mat <- read_expression_matrix(bulk_file, "Bulk (CPMcounts_ComBatSeq)")

############################################
## 3. Run NNLS deconvolution
############################################

prop_Z <- run_one_reference(signature_zeisel, bulk_mat, "Zeisel2015")
prop_T <- run_one_reference(signature_tasic,  bulk_mat, "Tasic2018")
prop_Y <- run_one_reference(signature_yao,    bulk_mat, "Yao2021")

############################################
## 4. Save results
############################################

write_props_csv(prop_Z, out_zeisel)
write_props_csv(prop_T, out_tasic)
write_props_csv(prop_Y, out_yao)

message("\nSaved NNLS proportions to:")
message(" - ", out_zeisel)
message(" - ", out_tasic)
message(" - ", out_yao)

############################################
## End of script
############################################
