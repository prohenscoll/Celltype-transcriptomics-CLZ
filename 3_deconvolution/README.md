## 3_deconvolution/

Preparation of single-cell reference matrices and cell-type deconvolution using **NNLS** and **CIBERSORTx**, including integration of cell-fraction estimates and visualization.

### Overview
This module is divided into three parts:

1) **Single-cell reference preparation (CIBERSORTx signature matrices)**  
   Generate reference expression matrices from three published scRNA-seq datasets:
   - Zeisel 2015
   - Tasic 2018
   - Yao 2021

2) **Deconvolution with NNLS**  
   Estimate cell-type proportions in bulk samples using non-negative least squares.

3) **Proportions analysis**  
   Summarize, compare, and visualize estimated cell fractions from CIBERSORTx and NNLS.

---

### Inputs
- Bulk raw counts matrix (gene-level):  
  `results/bulk/rawcounts.csv` (generated in `2_bulk_matrix/`)
- (Optional) bulk metadata (sample annotations): treatment, cell type, batch (if needed)

### Outputs
- Single-cell reference matrices (formatted for CIBERSORTx)
- CIBERSORTx signature matrices (one per reference)
- NNLS proportion estimates
- Plots and tables comparing fractions across methods and references

---

## Part 1 — Single-cell references for CIBERSORTx signature matrices

Scripts:
- `01_scReferenceFromZeisel2015.R`
- `02_scReferenceFromTasic2018.R`
- `03_scReferenceFromYao2021.R`

Each script:
- downloads/loads the scRNA-seq dataset (reference)
- harmonizes gene identifiers (symbols/Ensembl) to match bulk
- aggregates cells into cell-type profiles (e.g., average expression per cell type)
- exports a CIBERSORTx-ready reference matrix and annotation files (if needed)

Expected outputs (example):
- `results/deconvolution/references/Zeisel2015_reference.tsv`
- `results/deconvolution/references/Tasic2018_reference.tsv`
- `results/deconvolution/references/Yao2021_reference.tsv`

> **Note:** CIBERSORTx signature matrix generation is performed externally using the CIBERSORTx web/standalone workflow. This repository provides the input reference matrices and downstream analyses.

Run:
```bash
Rscript 01_scReferenceFromZeisel2015.R
Rscript 02_scReferenceFromTasic2018.R
Rscript 03_scReferenceFromYao2021.R

## Part 2 — NNLS DECONVOLUTION
