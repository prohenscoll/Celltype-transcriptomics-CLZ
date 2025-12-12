# Celltype-transcriptomics-CLZ for neurobiology of TRS

**Cell type specific transcriptional responses to clozapine in microglia, astrocytes, and neurons to investigate treatment-resistant schizophrenia**

This repository contains the full analysis pipeline used to characterize gene expression of microglia, astrocytes and neurons obtained from fontal cortices of mice treated with clozapine (CLZ). The study integrates bulk RNA-seq data and single-cell RNA-seq reference datasets and applies cell-type deconvolution to validate enrichment, followed by differential gene expression and gene set enrichment analyses.

---

## 📚 Project Overview

The goal of this project is to investigate the biological mechanisms underlying treatment-resistant schizophrenia (TRS), and to do that we characterized cell type specific transcriptional responses to clozapine (CLZ) in microglia, astrocytes, and neurons from treated-mice frontal cortices. The study uses bulk RNA-seq data from immunopanned cell populations and applies differential gene expression and gene set enrichment analyses to identify molecular programs associated with CLZ response. Immunopanned cell populations are validated with single-cell RNA-seq reference datasets and applying cell-type deconvolution.

- **Model system**: Mouse frontal cortex
- **Cell types**: Microglia, astrocytes, and neurons
- **Treatment**: Clozapine (CLZ) versus vehicle (and risperidone, where indicated)
- **Data**: Bulk RNA-seq integrated with single-cell RNA-seq reference datasets
- **Analyses**: Cell type deconvolution, differential gene expression, and gene set enrichment analysis

---

## 🧪 Analysis Pipeline

The analysis is structured across 6 scripts, executed in sequential order:

| Script | Description |
|--------|-------------|
| `1_data_preprocessing_and_descriptives.R` | Initial data cleaning, exclusion criteria, and exploratory analyses |
| `2_data_imputation.R` | Multiple imputation of missing values using MICE |
| `3_variable_selection_LASSO.R` | Feature selection via logistic LASSO with 10-fold cross-validation |
| `4_model_fitting.R` | Logistic regression and ML models fitted on imputed data (clinical, genetic, combined) |
| `5_model_evaluation.R` | Model calibration, discrimination (AUC), bootstrap-based comparisons and variable contribution|
| `6_decision_curve_analysis.R` | Evaluation of clinical utility using decision curve analysis (DCA) |

---

## 📦 Required Packages

Key R packages used include:

- `openxlsx`, `mice`, `VIM`, `naniar`
- `glmnet`, `dplyr`
- `rms`, `caret`, `pROC`, `boot`, `gbm`
- `ggplot2`, `ggsignif`, `fastshap`, `shapviz`
- `dcurves`, `gtsummary` 

You can install missing packages using:

```r
packages <- c("openxlsx", "mice", "VIM", "naniar", "glmnet", "dplyr", "rms", "caret", "pROC", "boot", "gbm", "ggplot2", "ggsignif", "fastshap", "shapviz", "dcurves", "gtsummary")
install.packages(setdiff(packages, rownames(installed.packages())))
```

---

## ⚙️ Reproducibility

To reproduce the full analysis:

1. Clone this repository:
   ```bash
   git clone https://github.com/laurajuliamelis/EOR-prediction-FEP.git
   ```

2. Open R or RStudio and set the working directory to the repo folder.

3. Run the scripts in the order specified above.

*Note: Raw data are not included in this repository due to privacy restrictions.*

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
