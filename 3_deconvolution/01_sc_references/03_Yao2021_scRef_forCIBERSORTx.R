############################################
# 03_Yao2021_scRef_forCIBERSORTx.R
# Single-cell RNA-seq reference construction (Yao et al., 2021)
#
# - Read metadata + raw counts from HDF5 (cells x genes)
# - Keep frontal cortex regions (ACA, AI, MOp, MOs_FRP, PL-ILA-ORB)
# - Subsample cells (to reduce runtime and neuron dominance)
# - Map subclass labels to 5 major brain cell classes
#   (Neuron, Astrocyte, Microglia, Oligodendrocyte, Endothelial)
# - Filter genes (remove all-zero, remove ERCC if present)
# - Export a gene × cell matrix formatted for CIBERSORTx
#
# Output format:
#   - First row: "GeneSymbol" + one phenotype label per cell
#   - Following rows: gene symbol + raw counts per cell
############################################

############################################
## 0. Configuration and libraries
############################################
# Input paths (relative to repository root)
in_dir    <- file.path("3_deconvolution", "01_sc_references", "inputs", "Yao2021")
meta_file <- file.path(in_dir, "metadata.csv")
h5_file   <- file.path(in_dir, "expression_matrix.hdf5")

# Output paths
out_dir  <- file.path("3_deconvolution", "01_sc_references", "outputs")
out_file <- file.path(out_dir, "Yao2021_frontal_SMALL_scRef_major5_forCIBERSORTx.tsv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Frontal cortex region filter (as used in your analysis)
regions_frontals <- c("ACA", "AI", "MOp", "MOs_FRP", "PL-ILA-ORB")

# Subsampling strategy (fast + balanced)
set.seed(123)
max_cells_per_subclass_initial <- 2000  # initial cap per subclass (after frontal filter)
max_neuron_per_subclass        <- 200   # cap per neuronal subclass to avoid neuron dominance
max_oligo_total                <- 1000  # cap total oligodendrocyte-lineage cells
max_microglia_total            <- 250   # optional cap (set to Inf to keep all)
max_endo_total                 <- Inf   # keep all by default
max_astro_total                <- Inf   # keep all by default

suppressPackageStartupMessages({
  library(Seurat)
  library(rhdf5)
  library(dplyr)
})

############################################
## 1. Input checks
############################################

if (!file.exists(meta_file)) stop("metadata.csv not found: ", meta_file)
if (!file.exists(h5_file))   stop("expression_matrix.hdf5 not found: ", h5_file)

message("Metadata: ", meta_file)
message("HDF5:      ", h5_file)
message("Output:    ", out_file)

############################################
## 2. Load metadata + HDF5 index vectors
############################################

meta <- read.csv(meta_file, row.names = 1, stringsAsFactors = FALSE)

if (!"exp_component_name" %in% colnames(meta)) stop("metadata must contain 'exp_component_name'")
if (!"region_label" %in% colnames(meta))       stop("metadata must contain 'region_label'")
if (!"subclass_label" %in% colnames(meta))     stop("metadata must contain 'subclass_label'")

meta$cell_id <- as.character(meta$exp_component_name)

# These datasets must exist in the HDF5 file
genes   <- h5read(h5_file, "/data/gene")
samples <- h5read(h5_file, "/data/samples")

message("HDF5 genes:   ", length(genes))
message("HDF5 samples: ", length(samples))

############################################
## 3. Align metadata with HDF5 samples
############################################

cells_common <- intersect(samples, meta$cell_id)
message("Common cells (samples ∩ metadata): ", length(cells_common))
if (length(cells_common) < 10) stop("Too few common cells. Check metadata cell IDs vs HDF5 samples.")

idx_samples <- match(cells_common, samples)
idx_meta    <- match(cells_common, meta$cell_id)

samples_sub <- samples[idx_samples]
meta_sub    <- meta[idx_meta, , drop = FALSE]

stopifnot(all(samples_sub == meta_sub$cell_id))

############################################
## 4. Restrict to frontal cortex regions
############################################

message("Region counts (aligned set):")
print(table(meta_sub$region_label))

meta_front <- meta_sub[meta_sub$region_label %in% regions_frontals, , drop = FALSE]
message("Cells in frontal regions: ", nrow(meta_front))
if (nrow(meta_front) < 10) stop("Too few frontal cells after filtering. Check region labels / filter list.")

############################################
## 5. Initial subsampling: <= 2000 cells per subclass (frontal only)
############################################

cells_by_subclass <- split(meta_front$cell_id, meta_front$subclass_label)

cells_keep <- unlist(lapply(cells_by_subclass, function(v) {
  sample(v, min(length(v), max_cells_per_subclass_initial))
}), use.names = FALSE)

message("Cells after initial cap (<= ", max_cells_per_subclass_initial, "/subclass): ", length(cells_keep))

# Convert cell IDs to HDF5 row indices (based on 'samples')
idx_cells_keep <- match(cells_keep, samples)
idx_cells_keep <- idx_cells_keep[!is.na(idx_cells_keep)]
message("Cells with valid HDF5 index: ", length(idx_cells_keep))
if (length(idx_cells_keep) < 10) stop("Too few cells with valid HDF5 index. Check matching logic.")

############################################
## 6. Read expression for selected cells (cells x genes)
############################################

expr <- h5read(
  file  = h5_file,
  name  = "/data/counts",
  index = list(idx_cells_keep, NULL)  # rows = cells, cols = genes
)

# Name columns/rows
stopifnot(length(genes) == ncol(expr))
colnames(expr) <- genes

cell_ids_keep <- samples[idx_cells_keep]
stopifnot(length(cell_ids_keep) == nrow(expr))
rownames(expr) <- cell_ids_keep

message("Expression loaded (cells x genes): ", paste(dim(expr), collapse = " x "))

############################################
## 7. Align metadata to expression rows
############################################

meta_keep <- meta_front[match(rownames(expr), meta_front$cell_id), , drop = FALSE]
stopifnot(all(rownames(expr) == meta_keep$cell_id))

############################################
## 8. Map subclasses to major5 phenotypes
############################################

sub <- meta_keep$subclass_label

major5 <- dplyr::case_when(
  # Astro
  sub %in% c("Astro") ~ "Astrocyte",
  # Microglia
  sub %in% c("Micro-PVM") ~ "Microglia",
  # Oligodendrocyte lineage
  sub %in% c("Oligo", "OPC") ~ "Oligodendrocyte",
  # Endothelial / mural-associated (includes SMC-Peri + VLMC)
  sub %in% c("Endo", "SMC-Peri", "VLMC") ~ "Endothelial",
  # Neurons (excitatory + inhibitory) and named subclasses
  sub %in% c(
    "Glutamatergic", "GABAergic",
    "L2 IT ENTl", "L2/3 IT CTX", "L2/3 IT PPP", "L2/3 IT RHP",
    "L4 RSP-ACA", "L4/5 IT CTX", "L5 IT CTX", "L5 PT CTX",
    "L5/6 IT TPE-ENT", "L5/6 NP CTX",
    "L6 CT CTX", "L6 IT CTX", "L6 IT ENTl",
    "L6b CTX", "L6b/CT ENT",
    "CT SUB", "CA2-IG-FC", "Car3", "Meis2",
    "NP PPP", "NP SUB",
    "Sst", "Sst Chodl", "Pvalb", "Vip", "Sncg", "CR", "Lamp5"
  ) ~ "Neuron",
  TRUE ~ NA_character_
)

message("major5 distribution (with NA):")
print(table(major5, useNA = "ifany"))

# Remove NA
keep <- !is.na(major5)
expr <- expr[keep, , drop = FALSE]
meta_keep <- meta_keep[keep, , drop = FALSE]
major5 <- major5[keep]

# Enforce >= 3 cells per phenotype (CIBERSORTx recommendation)
tab <- table(major5)
valid_types <- names(tab[tab >= 3])
keep2 <- major5 %in% valid_types

expr <- expr[keep2, , drop = FALSE]
meta_keep <- meta_keep[keep2, , drop = FALSE]
major5 <- major5[keep2]

message("major5 distribution after >=3 filter:")
print(table(major5))

############################################
## 9. Build SMALL reference (balanced subsampling)
############################################

idx <- seq_len(nrow(expr))

# Neurons: cap per subclass
idx_neuron <- which(major5 == "Neuron")
idx_neuron_keep <- integer(0)
if (length(idx_neuron) > 0) {
  sub_neuron <- meta_keep$subclass_label[idx_neuron]
  idx_by_subclass <- split(idx_neuron, sub_neuron)
  idx_neuron_keep <- unlist(lapply(idx_by_subclass, function(v) {
    sample(v, min(length(v), max_neuron_per_subclass))
  }), use.names = FALSE)
}

# Oligodendrocyte lineage: cap total
idx_oligo <- which(major5 == "Oligodendrocyte")
idx_oligo_keep <- if (length(idx_oligo) > 0) sample(idx_oligo, min(length(idx_oligo), max_oligo_total)) else integer(0)

# Microglia: optional cap total
idx_micro <- which(major5 == "Microglia")
cap_micro <- if (is.finite(max_microglia_total)) max_microglia_total else length(idx_micro)
idx_micro_keep <- if (length(idx_micro) > 0) sample(idx_micro, min(length(idx_micro), cap_micro)) else integer(0)

# Endothelial: optional cap total
idx_endo <- which(major5 == "Endothelial")
cap_endo <- if (is.finite(max_endo_total)) max_endo_total else length(idx_endo)
idx_endo_keep <- if (length(idx_endo) > 0) sample(idx_endo, min(length(idx_endo), cap_endo)) else integer(0)

# Astrocyte: optional cap total
idx_astro <- which(major5 == "Astrocyte")
cap_astro <- if (is.finite(max_astro_total)) max_astro_total else length(idx_astro)
idx_astro_keep <- if (length(idx_astro) > 0) sample(idx_astro, min(length(idx_astro), cap_astro)) else integer(0)

# Combine
idx_keep <- unique(c(idx_neuron_keep, idx_oligo_keep, idx_micro_keep, idx_endo_keep, idx_astro_keep))

expr_small <- expr[idx_keep, , drop = FALSE]
major5_small <- major5[idx_keep]

message("SMALL reference cells: ", nrow(expr_small))
message("SMALL major5 distribution:")
print(table(major5_small))

############################################
## 10. Filter genes and prepare gene×cell matrix
############################################

# Remove all-zero genes
gene_sums <- colSums(expr_small)
expr_small <- expr_small[, gene_sums > 0, drop = FALSE]

# Remove ERCC if present
expr_small <- expr_small[, !grepl("^ERCC-", colnames(expr_small)), drop = FALSE]

message("After gene filtering (cells x genes): ", paste(dim(expr_small), collapse = " x "))

# Transpose to genes x cells
mat <- t(expr_small)

# Gene name hygiene
gene_names <- rownames(mat)
gene_names <- ifelse(is.na(gene_names) | gene_names == "",
                     paste0("Gene_", seq_along(gene_names)),
                     gene_names)
rownames(mat) <- make.unique(gene_names)

# Labels for header
labels <- as.character(major5_small)
labels <- gsub("\\.", " ", labels)

# Warn if any phenotype has <3 cells after subsampling
lab_tab <- table(labels)
if (!all(lab_tab >= 3)) {
  message("WARNING: some phenotypes have <3 cells after subsampling:")
  print(lab_tab)
}

############################################
## 11. Write CIBERSORTx reference file
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

message("CIBERSORTx reference written to:")
message(out_file)

############################################
## 12. Quick sanity checks
############################################

preview <- readLines(out_file, n = 2)
message("Header preview:")
message(preview[1])

message("First gene row preview:")
message(preview[2])
