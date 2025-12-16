############################################
# 02_combatseq_pcr_normalization.R
# Batch correction of bulk raw counts using ComBat-Seq:
# - Loads raw count matrix (gene-level)
# - Builds sample metadata (batch, cell type, treatment) from sample names
# - PCA pre/post (logCPM) for QC
# - Applies ComBat-Seq preserving biological structure (treatment x celltype)
# - Saves corrected/CPM-corrected/logCPM-corrected matrices + PCA plots
############################################

############################################
## 0. Configuration and libraries
############################################
in_dir   <- "results/bulk"
out_dir  <- "results/bulk"
plot_dir <- "results/plots/combatseq"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
in_csv <- file.path(in_dir, "rawcounts.csv")

library(sva)  
library(dplyr)
library(ggplot2)
library(forcats)

############################################
## 1. Load raw counts
############################################
rawcounts_df <- read.csv (in_csv, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
colnames(rawcounts_df)[1] <- "Gene" 
rawcounts_df[, -1] <- lapply(rawcounts_df[, -1, drop = FALSE], function(x) as.numeric(x)) 

rawcounts <- as.matrix(rawcounts_df[, -1, drop = FALSE])
rownames(rawcounts) <- rawcounts_df$Gene

stopifnot(is.matrix(rawcounts))
stopifnot(nrow(rawcounts) > 0, ncol(rawcounts) > 0)
stopifnot(all(colSums(rawcounts) > 0))

############################################
## 2. Build coldata (batch, celltype, treatment)
############################################
sample_names <- colnames(rawcounts)

batch <- ifelse(grepl("_2_", sample_names), "2", "3")
celltype <- ifelse(grepl("CD45", sample_names), "Microglia",
                   ifelse(grepl("THY1", sample_names), "Neurons", "Astrocytes"))
treatment <- ifelse(grepl("VEH", sample_names), "VEH",
                    ifelse(grepl("RISPE", sample_names), "RISPE", "CLOZA"))

coldata <- data.frame(
  sample    = sample_names,
  treatment = factor(treatment, levels = c("VEH", "RISPE", "CLOZA")),
  batch     = factor(batch),
  celltype  = factor(celltype, levels = c("Microglia", "Astrocytes", "Neurons"))
)
rownames(coldata) <- sample_names

############################################
## 3. Pre-ComBat-Seq PCA (logCPM)
############################################
calc_logcpm_pca <- function(count_matrix, coldata) {
  
  stopifnot(all(colnames(count_matrix) %in% rownames(coldata)))
  
  lib_sizes <- colSums(count_matrix)
  cpm <- t(t(count_matrix) / lib_sizes * 1e6)
  logcpm <- log2(cpm + 1)
  
  var_genes <- apply(logcpm, 1, var)
  logcpm_var <- logcpm[var_genes > 0, , drop = FALSE]
  
  pca <- prcomp(t(logcpm_var), scale. = TRUE)
  percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
  
  pca_df <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  pca_df$sample <- rownames(pca_df)
  
  meta <- coldata[rownames(pca_df), , drop = FALSE]
  
  pca_df <- cbind(pca_df, meta) %>%
    mutate(
      treatment = fct_recode(treatment, RIS = "RISPE", CLZ = "CLOZA"),
      batch_label = ifelse(batch == "2", "Batch 1", "Batch 2")
    )
  
  list(pca = pca, pca_df = pca_df, percentVar = percentVar)
}

pre <- calc_logcpm_pca(rawcounts, coldata)

############################################
## 4. Apply ComBat-Seq
############################################
group_vec <- interaction(coldata$treatment, coldata$celltype, drop = TRUE)
combatseqcounts <- ComBat_seq(
  counts = rawcounts,
  batch  = coldata$batch,
  group  = group_vec
)
stopifnot(all(colSums(combatseqcounts) > 0))

############################################
## 5. Post-ComBat-Seq PCA (logCPM)
############################################
post <- calc_logcpm_pca(combatseqcounts, coldata)

############################################
## 6. PCA plots (all samples)
############################################
common_xlim <- c(-200, 200)
common_ylim <- c(-150, 150)

plot_pca <- function(pca_df, percentVar, title_text) {
  hull <- pca_df %>%
    group_by(celltype) %>%
    slice(chull(PC1, PC2))
  
  ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_point(size = 4, aes(shape = treatment), color = "grey40") +
    geom_text(aes(label = batch_label), color = "black", vjust = -1.2, size = 2.3) +
    geom_polygon(data = hull, alpha = 0.2, aes(fill = celltype, colour = celltype)) +
    scale_color_manual(values = c("VEH" = "#F8766D", "RIS" = "#00BA38", "CLZ" = "#619CFF")) +
    labs(
      title = title_text,
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance"),
      color = "Treatment") +
    coord_cartesian(xlim = common_xlim, ylim = common_ylim) +
    theme_minimal(base_size = 14)
}

p_pre  <- plot_pca(pre$pca_df,  pre$percentVar,  "PCA: pre-ComBat-Seq")
p_post <- plot_pca(post$pca_df, post$percentVar, "PCA: post-ComBat-Seq")

ggsave(filename = file.path(plot_dir, "PCA_all_pre_ComBatSeq.pdf"),  plot = p_pre,  width = 6, height = 5, dpi = 300)
ggsave(filename = file.path(plot_dir, "PCA_all_post_ComBatSeq.pdf"), plot = p_post, width = 6, height = 5, dpi = 300)

############################################
## 7. Save corrected COUNTS matrix
############################################
combatseqcounts_df <- data.frame(Gene = rownames(combatseqcounts), combatseqcounts, check.names = FALSE)
write.csv(combatseqcounts_df, file = file.path(out_dir, "rawcounts_ComBatSeq.csv"), row.names = FALSE)

############################################
## 8. Save CPM and logCPM matrices
############################################
lib_sizes_combat <- colSums(combatseqcounts)
CPM_combatseq <- t(t(combatseqcounts) / lib_sizes_combat * 1e6)
logCPM_combatseq <- log2(CPM_combatseq + 1)

write.table(data.frame(Gene = rownames(CPM_combatseq), CPM_combatseq, check.names = FALSE),
  file = file.path(out_dir, "CPMcounts_ComBatSeq.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(data.frame(Gene = rownames(logCPM_combatseq), logCPM_combatseq, check.names = FALSE),
  file = file.path(out_dir, "logCPMcounts_ComBatSeq.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


