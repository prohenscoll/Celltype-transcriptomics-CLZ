#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# One-time setup for RNA-seq preprocessing tools (macOS)
# Installs command-line tools needed for steps 01–05:
# FastQC, cutadapt, HISAT2, samtools, featureCounts (subread)
# ============================================================

echo "=== RNA-seq preprocessing setup (macOS) ==="

# 1) Check Homebrew
if ! command -v brew >/dev/null 2>&1; then
  echo "ERROR: Homebrew is not installed."
  echo "Install it with this command, then re-run this script:"
  echo '/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"'
  exit 1
fi

echo "Homebrew found: $(brew --version | head -n 1)"

# 2) Update Homebrew
echo "Updating Homebrew..."
brew update

# 3) Install tools
echo "Installing required tools..."
brew install fastqc
brew install cutadapt
brew install hisat2
brew install samtools
if command -v conda >/dev/null 2>&1; then
  conda install -c bioconda -y subread
else
  echo "ERROR: conda is not available. Install Miniconda/Anaconda first, or install subread via Homebrew (may require Xcode CLT)."
  exit 1
fi

# 4) Verify installations (print versions)
echo ""
echo "=== Versions ==="
fastqc --version || fastqc -v
cutadapt --version
hisat2 --version | head -n 1
samtools --version | head -n 1
featureCounts -v

echo ""
echo "Setup completed successfully."
echo "Next: run scripts/1_rna_preprocessing/01_fastqc.sh"
