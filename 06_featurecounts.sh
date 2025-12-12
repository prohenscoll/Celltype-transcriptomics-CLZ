#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Step 6: Gene-level quantification with featureCounts (paired-end)
# Input:  results/alignment/*.sorted.bam
# Output: results/counts/gene_counts.txt (+ .summary)
# Params (as used in the project): -p -O -t exon -g gene_id
# ============================================================

BAM_DIR="results/alignment"
COUNT_DIR="results/counts"
GTF="data/genome/Mus_musculus.GRCm39.110.gtf"
THREADS=4

mkdir -p "${COUNT_DIR}"

# Collect BAM files
BAMS=( "${BAM_DIR}"/*.sorted.bam )

if [[ ${#BAMS[@]} -eq 0 ]]; then
  echo "ERROR: No BAM files found in ${BAM_DIR} (*.sorted.bam)."
  exit 1
fi

if [[ ! -f "${GTF}" ]]; then
  echo "ERROR: GTF file not found: ${GTF}"
  exit 1
fi

echo "Running featureCounts..."
echo "BAM dir: ${BAM_DIR}"
echo "GTF:     ${GTF}"
echo "Output:  ${COUNT_DIR}/gene_counts.txt"
echo ""

featureCounts \
  -T "${THREADS}" \
  -p \
  -O \
  -t exon \
  -g gene_id \
  -a "${GTF}" \
  -o "${COUNT_DIR}/gene_counts.txt" \
  "${BAMS[@]}"

echo ""
echo "featureCounts completed."
echo "Counts:   ${COUNT_DIR}/gene_counts.txt"
echo "Summary:  ${COUNT_DIR}/gene_counts.txt.summary"
