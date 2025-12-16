############################################
# 01_prepare_rawcounts_matrix.R
# Bulk RNA-seq raw count matrix preparation (featureCounts -> merged gene-level matrix)
# - Merge featureCounts outputs across samples
# - Remove Ensembl version suffix
# - Map Ensembl IDs to gene symbols
# - Collapse duplicated symbols by summing counts
############################################

############################################
## 0. Configuration and libraries
############################################
counts_dir  <- "data/featureCounts"
out_dir     <- "results/bulk"
out_file   <- file.path(out_dir, "rawcounts.csv")
pattern_fc  <- "\\.txt$"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(dplyr)
library(org.Mm.eg.db)
library(AnnotationDbi)

############################################
## 1. Read and merge featureCounts outputs
############################################
read_featurecounts <- function(path) {
  tmp <- read.table(
    path,
    header = TRUE,
    sep = "\t",
    comment.char = "#",
    stringsAsFactors = FALSE
  )
  
  gene_ids <- tmp[["Geneid"]]
  counts   <- tmp[, ncol(tmp)]
  
  sample_name <- gsub(pattern_fc, "", path)
  
  out <- data.frame(counts, row.names = gene_ids)
  colnames(out) <- sample_name
  out
}

fc_files <- list.files(counts_dir, pattern = pattern_fc, full.names = TRUE)
stopifnot(length(fc_files) > 0)
count_list <- lapply(fc_files, read_featurecounts) 
common_genes <- Reduce(intersect, lapply(count_list, rownames)) 
raw_counts <- do.call(cbind, lapply(count_list, `[`, common_genes, , drop = FALSE)) 
ens_ids_nover <- gsub("\\.\\d+$", "", rownames(raw_counts))
rownames(raw_counts) <- ens_ids_nover             

############################################
## 2. Map Ensembl IDs to gene symbols and collapse to gene-level counts
############################################
symbols <- AnnotationDbi::mapIds(
  x = org.Mm.eg.db,
  keys = ens_ids_nover,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
) 

gene_labels <- ifelse(is.na(symbols), ens_ids_nover, symbols) 

df <- as.data.frame(raw_counts, check.names = FALSE) 
df$symbol <- gene_labels
df_collapsed <- df %>%
  dplyr::group_by(symbol) %>%
  dplyr::summarise(
    dplyr::across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

raw_counts_gene <- df_collapsed %>%
  tibble::column_to_rownames(var = "symbol") %>%
  as.matrix()
dim(raw_counts_gene)
head(raw_counts_gene)

############################################
## 3. Save outputs
############################################
write.csv(raw_counts_gene, file = out_file, row.names = TRUE)
