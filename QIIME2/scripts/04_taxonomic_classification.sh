#!/bin/bash

# QIIME2 Taxonomic Classification Script
# Step 4: Assign taxonomy using SILVA database

set -e  # Exit on error

echo "=================================================="
echo "QIIME2 STEP 4: TAXONOMIC CLASSIFICATION"
echo "Started: $(date)"
echo "=================================================="

# Check if input files exist
if [[ ! -f "rep-seqs.qza" ]]; then
    echo "❌ ERROR: rep-seqs.qza not found!"
    echo "Please run 03_dada2_denoise.sh first"
    exit 1
fi

if [[ ! -f "database/silva-138-99-v3v4-classifier.qza" ]]; then
    echo "❌ ERROR: SILVA classifier not found!"
    echo "Please download and train classifier first:"
    echo "Run: database/download_silva.sh"
    exit 1
fi

# Check if QIIME2 is available
if ! command -v qiime &> /dev/null; then
    echo "❌ ERROR: QIIME2 not found in PATH"
    echo "Please activate QIIME2 environment:"
    echo "conda activate qiime2-amplicon-2023.9"
    exit 1
fi

echo "🗄️  Using classifier: database/silva-138-99-v3v4-classifier.qza"
echo "🧬 Input sequences: rep-seqs.qza"

echo "🔬 Starting taxonomic classification..."

# Classify sequences using pre-trained classifier
time qiime feature-classifier classify-sklearn \
  --i-classifier database/silva-138-99-v3v4-classifier.qza \
  --i-reads rep-seqs.qza \
  --p-confidence 0.7 \
  --p-n-jobs 8 \
  --o-classification taxonomy.qza

if [[ $? -eq 0 ]]; then
    echo "✅ Taxonomic classification completed: taxonomy.qza"
else
    echo "❌ Taxonomic classification failed"
    exit 1
fi

echo "📊 Creating taxonomy visualization..."

# Create taxonomy visualization
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization visual/taxonomy.qzv

# Create interactive taxonomy barplot
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file metadata.txt \
  --o-visualization visual/taxa-bar-plots.qzv

echo "📈 Taxonomy Summary:"

# Export taxonomy for manual inspection
echo "📝 Exporting taxonomy table..."
qiime tools export \
  --input-path taxonomy.qza \
  --output-path taxonomy-export

# Show taxonomy overview
echo "🔍 Top taxonomic assignments:"
head -20 taxonomy-export/taxonomy.tsv

echo "📊 Classification confidence statistics:"
# This would require additional processing to get confidence distribution

# Log classification details
echo "📝 Logging classification details..."
qiime tools inspect taxonomy.qza > logs/taxonomy_inspection.txt

echo "=================================================="
echo "TAXONOMIC CLASSIFICATION COMPLETED SUCCESSFULLY"
echo "Finished: $(date)"
echo "=================================================="
echo ""
echo "📋 OUTPUT FILES:"
echo "  - taxonomy.qza                   : Taxonomic assignments"
echo "  - visual/taxonomy.qzv            : Taxonomy table visualization"
echo "  - visual/taxa-bar-plots.qzv      : Interactive bar plots"
echo "  - taxonomy-export/taxonomy.tsv   : Exported taxonomy table"
echo "  - logs/taxonomy_inspection.txt   : Taxonomy inspection log"
echo ""
echo "🔍 NEXT STEP: Run 05_filter_data.sh"
echo "=================================================="
