#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Step 5: Paired-end alignment with HISAT2 (sequential)
# Input:  data/trimmed_fastq/*_R1_trimmed.fastq.gz
# Output: results/alignment/*.sorted.bam (+ .bai + .flagstat.txt)
# Requires: hisat2, samtools
# ============================================================

INDEX_PREFIX="data/genome_index/mm39"
INPUT_DIR="data/trimmed_fastq"
OUTPUT_DIR="results/alignment"
THREADS=4

mkdir -p "${OUTPUT_DIR}"

echo "Starting HISAT2 alignment for all samples..."
echo "Input:  ${INPUT_DIR}"
echo "Index:  ${INDEX_PREFIX}"
echo "Output: ${OUTPUT_DIR}"
echo ""

for R1 in "${INPUT_DIR}"/*_R1_trimmed.fastq.gz; do
  SAMPLE=$(basename "${R1}" _R1_trimmed.fastq.gz)
  R2="${INPUT_DIR}/${SAMPLE}_R2_trimmed.fastq.gz"

  if [[ ! -f "${R2}" ]]; then
    echo "WARNING: Missing R2 for sample ${SAMPLE}. Skipping."
    continue
  fi

  echo "Aligning sample: ${SAMPLE}"

  # Align and directly create sorted BAM (avoid huge SAM files)
  hisat2 \
    -p "${THREADS}" \
    --dta \
    -x "${INDEX_PREFIX}" \
    -1 "${R1}" \
    -2 "${R2}" \
    | samtools view -bS - \
    | samtools sort -o "${OUTPUT_DIR}/${SAMPLE}.sorted.bam"

  # Index BAM (useful for IGV, downstream tools)
  samtools index "${OUTPUT_DIR}/${SAMPLE}.sorted.bam"

  # Alignment summary
  samtools flagstat "${OUTPUT_DIR}/${SAMPLE}.sorted.bam" \
    > "${OUTPUT_DIR}/${SAMPLE}.flagstat.txt"

  echo "Done: ${SAMPLE}"
  echo ""
done

echo "All alignments completed."
