############################################
# 05_NNLSandCIBERSORTx_proportionsanalysis.R
#
# - Read proportions from 3 references (Zeisel/Tasic/Yao)
# - Compute mean ± SD across references per sample & cell type
# - Plot stacked barplot
#
# Works for BOTH:
#  - NNLS outputs (CSV wide)
#  - CIBERSORTx outputs (CSV/TXT wide; may include extra columns like P-value)
############################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tibble)
})

############################################
## 0) Paths
############################################

# ---- NNLS outputs ----
nnls_dir <- file.path("results", "deconvolution", "nnls")
nnls_files <- list(
  Zeisel = file.path(nnls_dir, "NNLS_Zeisel2015_proportions.csv"),
  Tasic  = file.path(nnls_dir, "NNLS_Tasic2018_proportions.csv"),
  Yao    = file.path(nnls_dir, "NNLS_Yao2021_proportions.csv")
)

# ---- CIBERSORTx outputs ----
# Put your exported proportions here (recommended)
cib_dir <- file.path("results", "deconvolution", "cibersortx")
cib_files <- list(
  Zeisel = file.path(cib_dir, "CIBERSORTx_Zeisel2015_proportions.csv"),
  Tasic  = file.path(cib_dir, "CIBERSORTx_Tasic2018_proportions.csv"),
  Yao    = file.path(cib_dir, "CIBERSORTx_Yao2021_proportions.csv")
)

# ---- Output plot dir ----
plot_dir <- file.path("results", "deconvolution", "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

############################################
## 1) Manual sample order + short labels
############################################

sample_order <- c(
  "counts_S586_01_2_CD45_VEH_S1", "counts_S586_10_3_CD45_VEH_S10",
  "counts_S586_07_2_CD45_CLOZA_S7", "counts_S586_16_3_CD45_CLOZA_S16",
  "counts_S586_04_2_CD45_RISPE_S4", "counts_S586_13_3_CD45_RISPE_S13",
  "counts_S586_03_2_HEPA_VEH_S3", "counts_S586_12_3_HEPA_VEH_S12",
  "counts_S586_09_2_HEPA_CLOZA_S9", "counts_S586_18_3_HEPA_CLOZA_S18",
  "counts_S586_06_2_HEPA_RISPE_S6", "counts_S586_15_3_HEPA_RISPE_S15",
  "counts_S586_02_2_THY1_VEH_S2", "counts_S586_11_3_THY1_VEH_S11",
  "counts_S586_08_2_THY1_CLOZA_S8", "counts_S586_17_3_THY1_CLOZA_S17",
  "counts_S586_05_2_THY1_RISPE_S5", "counts_S586_14_3_THY1_RISPE_S14"
)

rename_samples <- c(
  "counts_S586_01_2_CD45_VEH_S1"    = "S1",
  "counts_S586_10_3_CD45_VEH_S10"   = "S2",
  "counts_S586_07_2_CD45_CLOZA_S7"  = "S3",
  "counts_S586_16_3_CD45_CLOZA_S16" = "S4",
  "counts_S586_04_2_CD45_RISPE_S4"  = "S5",
  "counts_S586_13_3_CD45_RISPE_S13" = "S6",
  "counts_S586_03_2_HEPA_VEH_S3"    = "S7",
  "counts_S586_12_3_HEPA_VEH_S12"   = "S8",
  "counts_S586_09_2_HEPA_CLOZA_S9"  = "S9",
  "counts_S586_18_3_HEPA_CLOZA_S18" = "S10",
  "counts_S586_06_2_HEPA_RISPE_S6"  = "S11",
  "counts_S586_15_3_HEPA_RISPE_S15" = "S12",
  "counts_S586_02_2_THY1_VEH_S2"    = "S13",
  "counts_S586_11_3_THY1_VEH_S11"   = "S14",
  "counts_S586_08_2_THY1_CLOZA_S8"  = "S15",
  "counts_S586_17_3_THY1_CLOZA_S17" = "S16",
  "counts_S586_05_2_THY1_RISPE_S5"  = "S17",
  "counts_S586_14_3_THY1_RISPE_S14" = "S18"
)

apply_manual_sample_order <- function(df_long, sample_order, rename_samples) {
  df_long <- df_long %>% mutate(sample = as.character(sample))
  
  missing_in_df <- setdiff(sample_order, unique(df_long$sample))
  if (length(missing_in_df) > 0) {
    warning("Some samples in 'sample_order' are missing from data (will be ignored):\n",
            paste(missing_in_df, collapse = "\n"))
  }
  
  present_order <- sample_order[sample_order %in% unique(df_long$sample)]
  
  df_long <- df_long %>%
    mutate(
      sample = factor(sample, levels = present_order),
      short_label = dplyr::recode(as.character(sample), !!!rename_samples)
    )
  
  # Ensure x-axis order is exactly S1..S18 in the manual order
  df_long$short_label <- factor(df_long$short_label, levels = unname(rename_samples[present_order]))
  
  df_long
}

############################################
## 2) Reading + cleaning proportions
############################################

assert_files_exist <- function(files_named_list, label) {
  missing <- files_named_list[!file.exists(unlist(files_named_list))]
  if (length(missing) > 0) {
    stop(
      label, " - missing file(s):\n",
      paste0(" - ", unlist(missing), collapse = "\n")
    )
  }
}

read_props_wide <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) {
    df <- read.delim(path, check.names = FALSE)
  } else {
    df <- read.csv(path, check.names = FALSE)
  }
  
  if (!"sample" %in% colnames(df)) {
    colnames(df)[1] <- "sample"
  }
  df$sample <- as.character(df$sample)
  df
}

clean_cibersortx_cols <- function(df) {
  exclude_cols <- c( "P-value", "Correlation", "RMSE", "Absolute score (sig.score)"
                     )
  keep <- setdiff(colnames(df), exclude_cols)
  df[, keep, drop = FALSE]
}

coerce_numeric_props <- function(df) {
  for (cn in setdiff(colnames(df), "sample")) {
    df[[cn]] <- suppressWarnings(as.numeric(df[[cn]]))
  }
  df
}

wide_to_long <- function(df, value_name) {
  df %>%
    pivot_longer(
      cols = -sample,
      names_to = "cell_type",
      values_to = value_name
    )
}

############################################
## 3) Core: average across 3 references
############################################

average_across_references <- function(files_named_list,
                                      method_name,
                                      is_cibersortx = FALSE,
                                      celltype_order = c("Astrocyte","Microglia","Neuron","Endothelial","Oligodendrocyte"),
                                      min_common_samples = 1) {
  
  assert_files_exist(files_named_list, method_name)
  
  wide_list <- lapply(files_named_list, function(p) {
    df <- read_props_wide(p)
    if (is_cibersortx) df <- clean_cibersortx_cols(df)
    df <- coerce_numeric_props(df)
    df
  })
  
  common_samples <- Reduce(intersect, lapply(wide_list, `[[`, "sample"))
  if (length(common_samples) < min_common_samples) {
    stop(method_name, ": too few common samples across references: ", length(common_samples))
  }
  
  long_list <- Map(
    function(df, ref) {
      df %>%
        filter(sample %in% common_samples) %>%
        wide_to_long(value_name = ref)
    },
    wide_list,
    names(files_named_list)
  )
  
  merged <- Reduce(
    function(x, y) inner_join(x, y, by = c("sample", "cell_type")),
    long_list
  )
  
  ref_cols <- names(files_named_list)
  
  merged <- merged %>%
    rowwise() %>%
    mutate(
      mean_prop = mean(c_across(all_of(ref_cols)), na.rm = TRUE),
      sd_prop   = sd(c_across(all_of(ref_cols)), na.rm = TRUE)
    ) %>%
    ungroup()
  
  merged$cell_type <- factor(
    merged$cell_type,
    levels = intersect(celltype_order, unique(merged$cell_type))
  )
  
  merged
}

############################################
## 4) Plotting barplot
############################################

plot_barplot <- function(df_long, title, out_pdf, width = 24, height = 8, colors = NULL) {
  p <- ggplot(df_long, aes(x = short_label, y = mean_prop, fill = cell_type)) +
    geom_bar(stat = "identity", width = 0.9, alpha = 0.95) +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(angle = 60, hjust = 0.5, size = 10),
      plot.title   = element_text(hjust = 0.5, size = 22),
      axis.text.y  = element_text(size = 14),
      axis.title.y = element_text(size = 16),
      legend.text  = element_text(size = 14),
      legend.title = element_text(size = 16)
    ) +
    labs(title = title, x = "", y = "Mean proportion", fill = "Cell type")
  
  if (!is.null(colors)) p <- p + scale_fill_manual(values = colors)
  
  ggsave(out_pdf, plot = p, width = width, height = height)
  p
}

############################################
## 5) Colors
############################################

colors_celltypes <- c(
  "Astrocyte"       = "#F8766D",
  "Endothelial"     = "#FDB515",
  "Microglia"       = "#00BA38",
  "Neuron"          = "#619CFF",
  "Oligodendrocyte" = "#FF99FF"
)

############################################
## 6) NNLS: average + manual order + barplot
############################################

nnls_merged <- average_across_references(
  nnls_files,
  method_name = "NNLS",
  is_cibersortx = FALSE
)

nnls_merged <- apply_manual_sample_order(nnls_merged, sample_order, rename_samples)

plot_barplot(
  nnls_merged,
  title   = "NNLS estimated cell-type proportions (mean across Zeisel/Tasic/Yao)",
  out_pdf = file.path(plot_dir, "Barplot_NNLS_mean.pdf"),
  colors  = colors_celltypes
)

############################################
## 7) CIBERSORTx: average + manual order + barplot
############################################

if (all(file.exists(unlist(cib_files)))) {
  
  cib_merged <- average_across_references(
    cib_files,
    method_name = "CIBERSORTx",
    is_cibersortx = TRUE
  )
  
  cib_merged <- apply_manual_sample_order(cib_merged, sample_order, rename_samples)
  
  plot_barplot(
    cib_merged,
    title   = "CIBERSORTx estimated cell-type proportions (mean across Zeisel/Tasic/Yao)",
    out_pdf = file.path(plot_dir, "Barplot_CIBERSORTx_mean.pdf"),
    colors  = colors_celltypes
  )
  
} else {
  message("\n[INFO] CIBERSORTx files not found -> skipping CIBERSORTx barplot.")
  message("Expected:")
  message(paste0(" - ", unlist(cib_files), collapse = "\n"))
}

message("\nDone. Barplots saved in: ", plot_dir)
