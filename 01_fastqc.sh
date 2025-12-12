#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Step 1: QC of raw FASTQ files with FastQC
# Input:  data/raw_fastq/*.fastq.gz
# Output: results/fastqc/raw/
# ============================================================

RAW_FASTQ_DIR="data/raw_fastq"
OUTPUT_DIR="results/fastqc/raw"
THREADS=4

# Check that FastQC is installed
command -v fastqc >/dev/null 2>&1 || {
  echo "ERROR: FastQC is not installed or not in PATH."
  echo "Install on Mac with: brew install fastqc"
  exit 1
}

# Create output folder
mkdir -p "${OUTPUT_DIR}"

# Run FastQC on all FASTQ.gz files
fastqc \
  --threads "${THREADS}" \
  --outdir "${OUTPUT_DIR}" \
  "${RAW_FASTQ_DIR}"/*.fastq.gz

echo "FastQC done. Results in: ${OUTPUT_DIR}"
