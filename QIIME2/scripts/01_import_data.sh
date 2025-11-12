#!/bin/bash

# Step 1: Import sequencing data into QIIME2

echo "=== Importing Sequencing Data ==="

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.csv \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33

# Create visualization
qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization visual/demux-summary.qzv

echo "Data import completed: paired-end-demux.qza"
echo "Visualization: visual/demux-summary.qzv"
