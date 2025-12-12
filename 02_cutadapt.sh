#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Step 2: Adapter trimming with cutadapt (paired-end RNA-seq)
# Library prep: QIAseq FastSelect RNA Library Kit
# Adapters: Illumina TruSeq
# Input:  data/raw_fastq/*_R1_001.fastq.gz
# Output: data/trimmed_fastq/
# ============================================================

RAW_FASTQ_DIR="data/raw_fastq"
TRIMMED_FASTQ_DIR="data/trimmed_fastq"
THREADS=4

# Illumina TruSeq adapter sequences
ADAPTER_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# Create output directory
mkdir -p "${TRIMMED_FASTQ_DIR}"

echo "Running cutadapt for paired-end FASTQ files..."

for R1 in "${RAW_FASTQ_DIR}"/*_R1_001.fastq.gz; do

  R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"

  if [[ ! -f "${R2}" ]]; then
    echo "WARNING: No matching R2 file for ${R1}. Skipping."
    continue
  fi

  SAMPLE_NAME=$(basename "${R1}" _R1_001.fastq.gz)

  OUT_R1="${TRIMMED_FASTQ_DIR}/${SAMPLE_NAME}_R1_trimmed.fastq.gz"
  OUT_R2="${TRIMMED_FASTQ_DIR}/${SAMPLE_NAME}_R2_trimmed.fastq.gz"

  echo "Trimming adapters for sample: ${SAMPLE_NAME}"

  cutadapt \
    -j "${THREADS}" \
    -a "${ADAPTER_R1}" \
    -A "${ADAPTER_R2}" \
    -m 20 \
    -o "${OUT_R1}" \
    -p "${OUT_R2}" \
    "${R1}" "${R2}"

done

echo "Cutadapt trimming completed successfully."

