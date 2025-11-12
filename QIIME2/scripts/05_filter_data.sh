#!/bin/bash

# QIIME2 Data Filtering Script
# Step 5: Filter low-frequency features and remove contaminants

set -e  # Exit on error

echo "=================================================="
echo "QIIME2 STEP 5: DATA FILTERING"
echo "Started: $(date)"
echo "=================================================="

# Check if input files exist
if [[ ! -f "table.qza" ]]; then
    echo "❌ ERROR: table.qza not found!"
    echo "Please run 03_dada2_denoise.sh first"
    exit 1
fi

if [[ ! -f "taxonomy.qza" ]]; then
    echo "❌ ERROR: taxonomy.qza not found!"
    echo "Please run 04_taxonomic_classification.sh first"
    exit 1
fi

# Check if QIIME2 is available
if ! command -v qiime &> /dev/null; then
    echo "❌ ERROR: QIIME2 not found in PATH"
    echo "Please activate QIIME2 environment:"
    echo "conda activate qiime2-amplicon-2023.9"
    exit 1
fi

echo "⚙️  Filtering Parameters:"
echo "  - Minimum frequency: 1"
echo "  - Include: Bacteria only"
echo "  - Exclude: Mitochondria, Chloroplast"

echo "🔧 Step 5.1: Filter low-frequency features..."

# Filter features with low frequency
qiime feature-table filter-features \
  --i-table table.qza \
  --p-min-frequency 1 \
  --o-filtered-table table-filtered-freq.qza

if [[ $? -eq 0 ]]; then
    echo "✅ Low-frequency filtering completed: table-filtered-freq.qza"
else
    echo "❌ Low-frequency filtering failed"
    exit 1
fi

echo "🔧 Step 5.2: Filter taxonomy (remove mitochondria/chloroplast)..."

# Filter table based on taxonomy
qiime taxa filter-table \
  --i-table table-filtered-freq.qza \
  --i-taxonomy taxonomy.qza \
  --p-include "d__Bacteria" \
  --p-exclude "mitochondria,chloroplast" \
  --o-filtered-table table-filtered-taxa.qza

if [[ $? -eq 0 ]]; then
    echo "✅ Taxonomy filtering completed: table-filtered-taxa.qza"
else
    echo "❌ Taxonomy filtering failed"
    exit 1
fi

echo "🔧 Step 5.3: Filter representative sequences..."

# Filter representative sequences to match filtered table
qiime feature-table filter-seqs \
  --i-data rep-seqs.qza \
  --i-table table-filtered-taxa.qza \
  --o-filtered-data rep-seqs-filtered.qza

if [[ $? -eq 0 ]]; then
    echo "✅ Sequence filtering completed: rep-seqs-filtered.qza"
else
    echo "❌ Sequence filtering failed"
    exit 1
fi

echo "🔧 Step 5.4: Re-classify filtered sequences..."

# Re-classify filtered sequences
qiime feature-classifier classify-sklearn \
  --i-classifier database/silva-138-99-v3v4-classifier.qza \
  --i-reads rep-seqs-filtered.qza \
  --o-classification taxonomy-filtered.qza

if [[ $? -eq 0 ]]; then
    echo "✅ Re-classification completed: taxonomy-filtered.qza"
else
    echo "❌ Re-classification failed"
    exit 1
fi

echo "📊 Creating filtered data visualizations..."

# Create summary of filtered table
qiime feature-table summarize \
  --i-table table-filtered-taxa.qza \
  --o-visualization visual/table-filtered.qzv

# Create taxonomy barplot for filtered data
qiime taxa barplot \
  --i-table table-filtered-taxa.qza \
  --i-taxonomy taxonomy-filtered.qza \
  --m-metadata-file metadata.txt \
  --o-visualization visual/taxa-bar-plots-filtered.qzv

echo "📈 Filtering Statistics:"

# Compare before and after filtering
echo "🔍 Feature count comparison:"
echo "  - Original features: $(qiime feature-table summarize --i-table table.qza --o-visualization /dev/null 2>&1 | grep "Feature count" | cut -d: -f2)"
echo "  - Filtered features: $(qiime feature-table summarize --i-table table-filtered-taxa.qza --o-visualization /dev/null 2>&1 | grep "Feature count" | cut -d: -f2)"

# Log filtering details
echo "📝 Logging filtering details..."
qiime tools inspect table-filtered-taxa.qza > logs/filtering_inspection.txt

echo "=================================================="
echo "DATA FILTERING COMPLETED SUCCESSFULLY"
echo "Finished: $(date)"
echo "=================================================="
echo ""
echo "📋 OUTPUT FILES:"
echo "  - table-filtered-taxa.qza        : Filtered feature table"
echo "  - rep-seqs-filtered.qza          : Filtered representative sequences"
echo "  - taxonomy-filtered.qza          : Filtered taxonomy"
echo "  - visual/table-filtered.qzv      : Filtered table visualization"
echo "  - visual/taxa-bar-plots-filtered.qzv : Filtered taxonomy bar plots"
echo ""
echo "🔍 NEXT STEP: Run 06_phylogenetic_tree.sh"
echo "=================================================="
