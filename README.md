# Celltype-transcriptomics-CLZ

**Cell type specific transcriptional responses to clozapine in microglia, astrocytes, and neurons to investigate treatment-resistant schizophrenia**

This repository contains the full analysis pipeline used to characterize gene expression of microglia, astrocytes and neurons obtained from fontal cortices of mice treated with clozapine (CLZ). The study integrates bulk RNA-seq data and single-cell RNA-seq reference datasets and applies cell-type deconvolution to validate enrichment, followed by differential gene expression and gene set enrichment analyses.

---

## 📚 Project Overview

The goal of this project is to investigate the biological mechanisms underlying treatment-resistant schizophrenia (TRS), and to do that we characterized cell type specific transcriptional responses to clozapine (CLZ) in microglia, astrocytes, and neurons from treated-mice frontal cortices. The study uses bulk RNA-seq data from immunopanned cell populations and applies differential gene expression and gene set enrichment analyses to identify molecular programs associated with CLZ response. Immunopanned cell populations are validated with single-cell RNA-seq reference datasets and applying cell-type deconvolution.

- **Cell types**: microglia, astrocytes or neurons.
- **Treatment**: vehicle (VEH), clozapine (CLZ) or risperidone (RIS).
- **Analyses**: Cell-type deconvolution, Differential Gene Expression (DEG) analysis, Gense Set Enrichment Analysis (GSEA).
  
---

## 🧪 Analysis Pipeline

The analysis is structured across sequential scripts executed in the following order:

| Script | Description |
|--------|-------------|
| `0_setup/` | Installing comand-line tools, check versions, create environments|
| `1_rna_preprocessing/` | Preprocessing of raw FASTQ files, including: quality control, adapter trimming, alignment to the mouse reference genome (GRCm39), gene-level quantification to obtain raw counts |
| `2_deconvolution_preparation.R` | Normalization of count data and preparation of single-cell reference matrices for signature matrix generation in CIBERSORTx |
| `3_deconvolution.R` | Cell type deconvolution using non-negative least squares (NNLS) and CIBERSORTx outputs, including integration of fraction estimates and generation of plots |
| `4_DataAnalysis.R` | DEG analysis and GSEA|

---

## 🛠 Setup (required tools and packages)

### Command-line tools: 

- Homebrew (macOS)
- FastQC
- cutadapt
- HISAT2
- samtools
- subread (featureCounts)

A one-time setup script is provided in `setup/`.

### R packages:

- **Core RNA-seq analysis**:  
  `DESeq2`, `sva`, `BiocGenerics`

- **Deconvolution**:  
  `nnls`

- **Single-cell reference data and handling**:  
  `scRNAseq`, `SingleCellExperiment`, `rhdf5`, `biomaRt`

- **Functional enrichment and GO analysis**:  
  `clusterProfiler`, `org.Mm.eg.db`, `AnnotationDbi`, `GO.db`
- **Data manipulation**:  
  `dplyr`, `tidyr`, `readr`, `tibble`, `stringr`, `purrr`, `reshape2`, `tidyverse`

- **Visualization**:  
  `ggplot2`, `pheatmap`, `RColorBrewer`, `ggrepel`, `VennDiagram`, `ggvenn`, `grid`

- **Input / Output utilities**:  
  `openxlsx`
  
You can install missing packages using:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2","sva","BiocGenerics","clusterProfiler","org.Mm.eg.db", "AnnotationDbi", "GO.db", "scRNAseq", "SingleCellExperiment", "rhdf5", "biomaRt"))
install.packages(c("nnls", "dplyr", "tidyr", "readr", "tibble", "stringr", "purrr", "reshape2", "tidyverse", "ggplot2", "pheatmap", "RColorBrewer", "ggrepel", "VennDiagram", "ggvenn", "grid", "openxlsx"))
```

---

## ⚙️ Reproducibility

To reproduce the full analysis:

1. Clone this repository:
   ```bash
   git clone https://github.com/prohenscoll/Celltype-transcriptomics-CLZ.git
   ```

2. Open R or RStudio and set the working directory to the repo folder.

3. Run the scripts in the order specified above.

*Note: Raw data are in Gene Omnibus under the code GSEXXXX.*

---

## 📄 License

This code is released under the

---

## 📬 Contact

For questions, suggestions, or collaboration opportunities, please contact:

**Sergi Mas Herrero**  
Accredited Researcher · CIBERSAM
Associate Professor · University of Barcelona  
📧 sergimash@ub.edu

and 

**Lucía Prohens Coll**  
Collaborative Researcher ·
PhD Candidate · University of Barcelona  
📧 prohens@ub.edu
