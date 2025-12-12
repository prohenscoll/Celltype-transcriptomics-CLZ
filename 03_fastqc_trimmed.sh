#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Step 3: Quality control of trimmed FASTQ files with FastQC
# Input:  data/trimmed_fastq/*_trimmed.fastq.gz
# Output: results/fastqc/trimmed/
# ============================================================

TRIMMED_FASTQ_DIR="data/trimmed_fastq"
OUTPUT_DIR="results/fastqc/trimmed"
THREADS=4

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "Running FastQC on trimmed FASTQ files..."

fastqc \
  --threads "${THREADS}" \
  --outdir "${OUTPUT_DIR}" \
  "${TRIMMED_FASTQ_DIR}"/*.fastq.gz

echo "FastQC on trimmed files completed."
echo "Results saved in: ${OUTPUT_DIR}"
