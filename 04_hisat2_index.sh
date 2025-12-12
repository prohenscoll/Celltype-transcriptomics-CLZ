#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Step 4: Download mouse reference genome and build HISAT2 index
# Genome: Mus musculus GRCm39 (mm39)
# ============================================================

GENOME_DIR="data/genome"
INDEX_DIR="data/genome_index"
THREADS=4

mkdir -p "${GENOME_DIR}"
mkdir -p "${INDEX_DIR}"

cd "${GENOME_DIR}"

echo "Downloading mouse reference genome (GRCm39)..."
wget -c ftp://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip -f Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

echo "Downloading gene annotation (GTF)..."
wget -c ftp://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz
gunzip -f Mus_musculus.GRCm39.110.gtf.gz

echo "Building HISAT2 index..."
hisat2-build -p "${THREADS}" \
  Mus_musculus.GRCm39.dna.primary_assembly.fa \
  "${INDEX_DIR}/mm39"

echo "Genome indexing completed."
