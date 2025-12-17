############################################
# 02_Tasic2018_scRef_forCIBERSORTx.R
# Single-cell RNA-seq reference construction (Tasic et al., 2018)
#
# - Load Tasic 2018 scRNA-seq dataset
# - Map  cell types to 5 major brain cell classes
#   (Neuron, Astrocyte, Microglia, Oligodendrocyte, Endothelial)
# - Extract raw (non-log) gene expression counts
# - Filter genes and cells according to CIBERSORTx requirements
# - Export a gene × cell matrix formatted for CIBERSORTx
#
# Output format:
#   - First row: "GeneSymbol" + one phenotype label per cell
#   - Following rows: gene symbol + raw counts per cell
############################################

############################################
## 0. Configuration and libraries
############################################
out_dir  <- file.path("3_deconvolution", "01_sc_references", "outputs")
out_file <- file.path(out_dir, "Tasic2018_scRef_major5_forCIBERSORTx.tsv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scRNAseq)
  library(SingleCellExperiment)
  library(Seurat)
  library(tidyr)
  library(dplyr)
})

############################################
## 1. Load Tasic 2018 scRNA-seq dataset
############################################
tasic <- TasicBrainData()
message("Tasic: available colData columns:")
print(colnames(colData(tasic)))

############################################
## 2. Convert to Seurat object
############################################
seurat_obj <- as.Seurat(tasic, counts = "counts", data = NULL)
DefaultAssay(seurat_obj) <- "originalexp"
table(seurat_obj$broad_type)

############################################
## 3. Map broad types to major cell types (major5)
############################################
seurat_obj$major5 <- dplyr::case_when(
  seurat_obj$broad_type %in% c("GABA-ergic Neuron", "Glutamatergic Neuron") ~ "Neuron",
  seurat_obj$broad_type %in% c("Astrocyte") ~ "Astrocyte",
  seurat_obj$broad_type %in% c("Oligodendrocyte", "Oligodendrocyte Precursor Cell") ~ "Oligodendrocyte",
  seurat_obj$broad_type %in% c("Microglia") ~ "Microglia",
  seurat_obj$broad_type %in% c("Endothelial Cell") ~ "Endothelial",
  TRUE ~ NA_character_
)

message("Major cell type distribution (including NA):")
print(table(seurat_obj$major5, useNA = "ifany"))

############################################
## 4. Filter cells and enforce CIBERSORTx constraints
############################################
seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[!is.na(seurat_obj$major5)]) # Remove cells without a valid major cell type
tbl <- table(seurat_obj$major5)
print(tbl)
stopifnot(all(tbl >= 3))

############################################
## 5. Extract raw counts and filter genes
############################################
assay_use <- DefaultAssay(seurat_obj)
if ("counts" %in% Layers(seurat_obj, assay = assay_use)) {
  mat <- GetAssayData(seurat_obj, assay = assay_use, layer = "counts") 
} else {
  mat <- GetAssayData(seurat_obj, assay = assay_use, slot = "counts") 
}
mat <- mat[rowSums(mat) > 0, ]
mat <- mat[!grepl("^ERCC-", rownames(mat)), ]
message("Expression matrix dimensions (genes x cells): ",
        paste(dim(mat), collapse = " x "))

############################################
## 6. Prepare phenotype labels
############################################
labels <- seurat_obj$major5[colnames(mat)]
labels <- gsub("\\.", " ", labels)   # CIBERSORTx recommendation
print(table(labels))

############################################
## 7. Ensure valid and unique gene symbols
############################################
genes <- rownames(mat)
genes <- ifelse(is.na(genes) | genes == "", paste0("Gene_", seq_along(genes)), genes)
genes <- make.unique(genes)
rownames(mat) <- genes

############################################
## 8. Write CIBERSORTx reference file
############################################
cat(paste(c("GeneSymbol", labels), collapse = "\t"), file = out_file)
cat("\n", file = out_file, append = TRUE)

write.table(
  cbind(GeneSymbol = rownames(mat), as.matrix(mat)),
  file = out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  append = TRUE
)

message("Saved: ", out_file)

############################################
## 9. Checks
############################################
preview <- readLines(out_file, n = 2)
message("Header preview:")
message(preview[1])

message("First gene row preview:")
message(preview[2])





