#!/bin/bash

# SILVA Database Download Script
# Downloads and prepares SILVA database for QIIME2

echo "=== Downloading SILVA Database ==="

# Create database directory
mkdir -p silva
cd silva

# Download SILVA files (example URLs - update with actual download links)
echo "1. Downloading SILVA sequences..."
wget -O silva-138-99-seqs.qza "https://data.qiime2.org/2023.9/common/silva-138-99-seqs.qza"

echo "2. Downloading SILVA taxonomy..."
wget -O silva-138-99-tax.qza "https://data.qiime2.org/2023.9/common/silva-138-99-tax.qza"

echo "3. Extracting V3-V4 specific sequences..."
qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer ACTCCTACGGGAGGCAGCAG \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --o-reads silva-138-99-v3v4-seqs.qza

echo "4. Training classifier..."
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-99-v3v4-seqs.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
  --o-classifier silva-138-99-v3v4-classifier.qza

echo "=== SILVA Database Setup Completed ==="
echo "Classifier ready: silva/silva-138-99-v3v4-classifier.qza"
