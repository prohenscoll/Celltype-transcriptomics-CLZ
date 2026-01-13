# Celltype-transcriptomics-CLZ

**Cell type specific transcriptional responses to clozapine in microglia, astrocytes, and neurons to investigate treatment-resistant schizophrenia**

This repository contains the full analysis pipeline used to characterize gene expression profiles of microglia, astrocytes and neurons obtained from the frontal cortices of mice treated with clozapine (CLZ). The study integrates bulk RNA-seq data and single-cell RNA-seq reference datasets and applies cell-type deconvolution to validate enrichment, followed by differential gene expression and gene set enrichment analyses.

---

## 📚 Project Overview

The goal of this project is to investigate the biological mechanisms underlying treatment-resistant schizophrenia (TRS). To achieve this, we characterized cell-type specific transcriptional responses to clozapine (CLZ) in microglia, astrocytes, and neurons from treated mice. The study uses bulk RNA-seq data from immunopanned cell populations and applies:
· **Differential gene expression (DEG) analysis**: to identify genes uniquely regulated by treatment
· **Gene Set Enrichment Analysis (GSEA)**: to identify molecular programs associated with CLZ response. 
· **Validation**: Immunopanned populations are validated with single-cell RNA-seq reference datasets and cell-type deconvolution.

**Study Design**:
- **Cell types**: Microglia (CD45+), Astrocytes (HEPA+) and Neurons (THY1+).
- **Treatment**: Vehicle (VEH), Clozapine (CLZ) or Risperidone (RIS).
- **Analyses**: Cell-type deconvolution, DEG analysis and GSEA/GO.
  
---

## 🧪 Analysis Pipeline

The analysis is structured into sequential modules:

| Folder | Description |
|--------|-------------|
| `0_setup/` | Installation of required command-line tools, version checks, and environment setup|
| `1_rna_preprocessing/` | Raw FASTQ files preprocessing: quality control (QC), adapter trimming, alignment to the mouse genome (GRCm39), and gene-level quantification|
| `2_bulk_matrix_preparation/` | Generation of a clean bulk gene expression matrix, handling batch effects (e.g., ComBat-seq) and normalization (CPM) for downstream analyses|
| `3_deconvolution/` | Preparation of scRNA-seq reference matrices and cell-type deconvolution using non-negative least squares (NNLS) and CIBERSORTx, including integration of cell-fraction estimates and visualization|
| `4_dataAnalysis/` | Core statistical analysis: DEG (DESeq2), GSEA (clusterProfiler) and GO term filtering|

---

## 🛠 Setup (required tools and packages)

### Command-line tools: 
The preprocessing pipeline requires the following tools (installable via Homebrew on macOS):
- FastQC
- cutadapt
- HISAT2
- samtools
- subread (featureCounts)

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
  
The analysis depends on several Bioconductor and CRAN packages. You can install all dependencies by running:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
# Bioconductor packages
bioc_pkgs <- c("DESeq2", "sva", "BiocGenerics", "clusterProfiler", "org.Mm.eg.db", 
               "AnnotationDbi", "GO.db", "scRNAseq", "SingleCellExperiment", 
               "rhdf5", "biomaRt")
# CRAN packages
cran_pkgs <- c("nnls", "tidyverse", "pheatmap", "RColorBrewer", "ggrepel", 
               "VennDiagram", "ggvenn", "openxlsx", "reshape2")
BiocManager::install(bioc_pkgs)
install.packages(cran_pkgs)
```
---

## ⚙️ Reproducibility

To reproduce the analysis:

1. Clone the repository:
   ```bash
   git clone https://github.com/prohenscoll/Celltype-transcriptomics-CLZ.git
   ```

2. Set up the environment: Run the scripts in 0_setup/ to ensure all command-line tools are available.

3. Data Access: Raw data are available in the NCBI Gene Expression Omnibus under accession number GSE315204.

4. CIBERSORTx: Please note that deconvolution via CIBERSORTx was performed using the online platform (https://cibersortx.stanford.edu/) according to the developers' guidelines.
   
5. Run the scripts in the order specified above.

---

## 📄 License

This project is licensed under the MIT License. See the LICENSE file for details

---

## 📬 Contact

For questions, suggestions, or collaboration opportunities, please contact:

**Sergi Mas Herrero**  
Accredited Researcher · CIBERSAM · IDIBAPS
Associate Professor · University of Barcelona  
📧 sergimash@ub.edu

**Natalia Rodriguez Ferret**  
Accredited Researcher · CIBERSAM · IDIBAPS
Associate Professor · University of Barcelona  
📧 nrodriguezfe@ub.edu

**Lucía Prohens Coll**  
Collaborative Researcher ·
PhD Candidate · University of Barcelona  
📧 prohens@ub.edu
