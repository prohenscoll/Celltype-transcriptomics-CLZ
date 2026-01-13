############################################
# 01_DESeq2_GSEA_allCellTypes.R
#
# Differential expression (DESeq2) + Functional Enrichment GSEA (clusterProfiler)
#
# Process:
# 1. Load ComBat-Seq corrected counts
# 2. Loop through cell types (Microglia, Astrocyte, Neuron)
# 3. Perform DESeq2 (CLZ vs VEH and RIS vs VEH)
# 4. Run GSEA using clusterProfiler (GO Biological Process)
# 5. Export tables and plots
############################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(openxlsx)
  library(readr)
  library(tibble)
})

############################################
## 0. Configuration & Paths
############################################
input_file <- "results/bulk/rawcounts_ComBatSeq.csv"
out_dir    <- "results/data_analysis"
plot_dir   <- "results/plots/differential_expression"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Define cell types and their corresponding sample patterns
cell_types <- list(
  Microglia = "CD45",
  Astrocyte = "HEPA",
  Neuron    = "THY1"
)

############################################
## 1. Load Data
############################################
raw_counts <- read.csv(input_file, row.names = 1, check.names = FALSE)

# Build metadata dynamically from column names
metadata <- data.frame(SampleID = colnames(raw_counts)) %>%
  mutate(
    CellType  = case_when(
      grepl("CD45", SampleID) ~ "Microglia",
      grepl("HEPA", SampleID) ~ "Astrocyte",
      grepl("THY1", SampleID) ~ "Neuron"
    ),
    Treatment = case_when(
      grepl("VEH", SampleID)   ~ "Vehicle",
      grepl("CLOZA", SampleID) ~ "Clozapine",
      grepl("RISPE", SampleID) ~ "Risperidone"
    )
  ) %>%
  column_to_rownames("SampleID")

metadata$Treatment <- factor(metadata$Treatment, levels = c("Vehicle", "Clozapine", "Risperidone"))

############################################
## 2. Main Analysis Loop
############################################
for (ct in names(cell_types)) {
  message("\n>>> Processing Cell Type: ", ct)
  
  # 2.1 Subset data for current cell type
  ct_samples <- rownames(metadata[metadata$CellType == ct, ])
  ct_counts  <- raw_counts[, ct_samples]
  ct_meta    <- metadata[ct_samples, ]
  
  # 2.2 DESeq2 Pipeline
  dds <- DESeqDataSetFromMatrix(countData = round(ct_counts),
                                colData = ct_meta,
                                design = ~ Treatment)
  dds <- DESeq(dds)
  
  # 2.3 Comparisons: Clozapine and Risperidone vs Vehicle
  comparisons <- list(Clozapine = "Clozapine", Risperidone = "Risperidone")
  
  for (comp_name in names(comparisons)) {
    message("Running comparison: ", comp_name, " vs Vehicle")
    
    res <- results(dds, contrast = c("Treatment", comp_name, "Vehicle"))
    res_df <- as.data.frame(res) %>% 
      rownames_to_column("gene") %>%
      arrange(pvalue)
    
    # Save DEGs table
    write.csv(res_df, file.path(out_dir, paste0("DEGs_", ct, "_", comp_name, ".csv")), row.names = FALSE)
    
    # 2.4 Volcano Plot
    p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 0.5), alpha = 0.5) +
      scale_color_manual(values = c("black", "red")) +
      geom_text_repel(data = head(res_df, 10), aes(label = gene)) +
      labs(title = paste("Volcano Plot:", ct, "-", comp_name),
           subtitle = "Significant genes highlighted in red") +
      theme_minimal() + theme(legend.position = "none")
    
    ggsave(file.path(plot_dir, paste0("Volcano_", ct, "_", comp_name, ".pdf")), p_volcano, width = 7, height = 6)
    
    # 2.5 GSEA Analysis
    message("Running GSEA for: ", ct, " - ", comp_name)
    
    # Create gene list ranked by Wald Statistic
    gene_list <- res_df$stat
    names(gene_list) <- res_df$gene
    gene_list <- sort(gene_list, decreasing = TRUE)
    
    gsea_res <- gseGO(geneList     = gene_list,
                      OrgDb        = org.Mm.eg.db,
                      keyType      = "SYMBOL",
                      ont          = "BP",
                      minGSSize    = 10,
                      maxGSSize    = 500,
                      pvalueCutoff = 0.05,
                      seed         = 123)
    
    if (nrow(gsea_res) > 0) {
      # Save GSEA results
      write.xlsx(as.data.frame(gsea_res), file.path(out_dir, paste0("GSEA_", ct, "_", comp_name, ".xlsx")))
      
      # GSEA Dotplot
      p_gsea <- dotplot(gsea_res, showCategory = 15) + 
        labs(title = paste("GSEA BP:", ct, "-", comp_name))
      
      ggsave(file.path(plot_dir, paste0("GSEA_Dotplot_", ct, "_", comp_name, ".pdf")), p_gsea, width = 9, height = 7)
    } else {
      message("No significant GO terms found for ", ct, " ", comp_name)
    }
  }
}

message("\n✅ Data Analysis Completed. Results in: ", out_dir)