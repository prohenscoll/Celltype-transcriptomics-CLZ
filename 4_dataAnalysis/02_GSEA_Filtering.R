############################################
# 02_GSEA_Filtering_and_Cytoscape.R
# Refinement of GSEA results:
# 1. Significance & Size filtering
# 2. NES Top 40% Selection
# 3. GO Depth analysis (removing too general/specific terms)
# 4. Cytoscape-ready export
############################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(GO.db)
  library(ggplot2)
  library(openxlsx)
  library(tidyr)
})

############################################
## 0. Configuration & Paths
############################################
input_dir  <- "results/data_analysis"
output_dir <- "results/data_analysis/filtered_gsea"
plot_dir   <- "results/plots/gsea_refinement"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Define the files to process based on outputs from script 01
gsea_files <- list.files(input_dir, pattern = "^GSEA_.*\\.xlsx$", full.names = TRUE)

############################################
## 1. Helper Functions
############################################

# Function to get GO term depth (number of ancestors)
get_go_depth <- function(go_id) {
  anc <- tryCatch(GOBPANCESTOR[[go_id]], error = function(e) NA)
  if (is.null(anc) || all(is.na(anc))) return(NA)
  return(length(anc))
}

# Core function to filter GSEA results
refine_gsea <- function(file_path) {
  base_name <- tools::file_path_sans_ext(basename(file_path))
  message("Refining: ", base_name)
  
  # A. Load data
  df <- read_excel(file_path)
  if(nrow(df) == 0) return(NULL)
  
  # B. Filter by Strong Significance & GeneSet Size
  df_filt <- df %>%
    filter(p.adjust < 0.005) %>%
    filter(setSize > 15, setSize < 300)
  
  if(nrow(df_filt) == 0) {
    message("No terms passed initial filters for ", base_name)
    return(NULL)
  }
  
  # C. Filter by Top 40% Absolute NES (Impact)
  topN <- ceiling(0.4 * nrow(df_filt))
  df_filt <- df_filt %>%
    mutate(abs_NES = abs(NES)) %>%
    arrange(desc(abs_NES)) %>%
    slice_head(n = topN) %>%
    select(-abs_NES)
  
  # D. Depth Analysis (Specificity)
  df_filt$depth <- sapply(df_filt$ID, get_go_depth)
  
  # Create a QC Plot for Depth distribution
  p_depth <- ggplot(df_filt, aes(x = depth)) +
    geom_histogram(binwidth = 1, fill = "#619CFF", color = "white") +
    labs(title = paste("GO Depth Distribution:", base_name), x = "Depth", y = "Count") +
    theme_minimal()
  ggsave(file.path(plot_dir, paste0("Depth_Dist_", base_name, ".pdf")), p_depth, width = 6, height = 4)
  
  # Filter by depth: remove too general (<=4) and too specific (>40)
  df_filt <- df_filt %>%
    filter(!is.na(depth), depth > 4, depth <= 40)
  
  # E. Prepare for Cytoscape (EnrichmentMap format)
  df_export <- df_filt %>%
    select(ID, Description, pvalue, p.adjust, NES, core_enrichment) %>%
    rename(GeneSet = ID,
           FDR_q_value = p.adjust,
           Genes = core_enrichment) %>%
    mutate(Genes = gsub("/", ",", Genes)) # Cytoscape friendly
  
  # F. Save results
  write.table(df_export, 
              file = file.path(output_dir, paste0(base_name, "_Cytoscape.txt")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(df_export)
}

############################################
## 2. Execution
############################################

# Apply the refinement to all GSEA files found
refined_list <- lapply(gsea_files, refine_gsea)

message("\ GSEA Refinement completed.")
message("Filtered files ready for Cytoscape in: ", output_dir)