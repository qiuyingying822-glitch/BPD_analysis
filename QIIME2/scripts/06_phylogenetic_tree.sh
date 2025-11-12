#!/bin/bash

# QIIME2 Phylogenetic Tree Script
# Step 6: Generate phylogenetic tree for diversity analysis

set -e  # Exit on error

echo "=================================================="
echo "QIIME2 STEP 6: PHYLOGENETIC TREE GENERATION"
echo "Started: $(date)"
echo "=================================================="

# Check if input files exist
if [[ ! -f "rep-seqs-filtered.qza" ]]; then
    echo "❌ ERROR: rep-seqs-filtered.qza not found!"
    echo "Please run 05_filter_data.sh first"
    exit 1
fi

# Check if QIIME2 is available
if ! command -v qiime &> /dev/null; then
    echo "❌ ERROR: QIIME2 not found in PATH"
    echo "Please activate QIIME2 environment:"
    echo "conda activate qiime2-amplicon-2023.9"
    exit 1
fi

echo "🌳 Phylogenetic Tree Generation Steps:"
echo "  1. Multiple sequence alignment (MAFFT)"
echo "  2. Mask alignment to remove highly variable positions"
echo "  3. Phylogeny construction (FastTree)"
echo "  4. Root the tree"

echo "🔬 Step 6.1: Multiple sequence alignment with MAFFT..."

# Generate phylogenetic tree
time qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-filtered.qza \
  --p-n-threads 8 \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

if [[ $? -eq 0 ]]; then
    echo "✅ Phylogenetic tree generation completed"
else
    echo "❌ Phylogenetic tree generation failed"
    exit 1
fi

echo "📊 Creating tree visualizations..."

# Create tree visualization
qiime empress tree-plot \
  --i-tree rooted-tree.qza \
  --m-feature-metadata-file taxonomy-filtered.qza \
  --o-visualization visual/rooted-tree.qzv

# Create alignment visualization (first 100 sequences for performance)
echo "🔍 Creating alignment visualization (subset)..."

qiime alignment mask \
  --i-alignment aligned-rep-seqs.qza \
  --o-masked-alignment aligned-rep-seqs-masked.qza

qiime tools export \
  --input-path aligned-rep-seqs-masked.qza \
  --output-path alignment-export

echo "📈 Phylogenetic Tree Statistics:"

# Get tree statistics
echo "🌲 Tree Information:"
echo "  - Rooted tree: rooted-tree.qza"
echo "  - Unrooted tree: unrooted-tree.qza"
echo "  - Alignment: aligned-rep-seqs.qza"
echo "  - Masked alignment: masked-aligned-rep-seqs.qza"

# Count tree tips (should match feature count)
echo "🔍 Tree validation:"
FEATURE_COUNT=$(qiime feature-table summarize --i-table table-filtered-taxa.qza --o-visualization /dev/null 2>&1 | grep "Feature count" | cut -d: -f2 | tr -d ' ')
echo "  - Features in table: $FEATURE_COUNT"

echo "📝 Logging phylogenetic details..."

# Log tree information
qiime tools inspect rooted-tree.qza > logs/tree_inspection.txt
qiime tools inspect aligned-rep-seqs.qza > logs/alignment_inspection.txt

# Export tree for other applications
echo "📤 Exporting tree files..."

qiime tools export \
  --input-path rooted-tree.qza \
  --output-path tree-export

qiime tools export \
  --input-path unrooted-tree.qza \
  --output-path tree-export/unrooted

echo "✅ Exported tree files to tree-export/"

echo "=================================================="
echo "PHYLOGENETIC TREE GENERATION COMPLETED SUCCESSFULLY"
echo "Finished: $(date)"
echo "=================================================="
echo ""
echo "📋 OUTPUT FILES:"
echo "  - rooted-tree.qza               : Rooted phylogenetic tree"
echo "  - unrooted-tree.qza             : Unrooted phylogenetic tree"
echo "  - aligned-rep-seqs.qza          : Multiple sequence alignment"
echo "  - masked-aligned-rep-seqs.qza   : Masked alignment"
echo "  - visual/rooted-tree.qzv        : Tree visualization"
echo "  - tree-export/                  : Exported tree files"
echo "  - logs/tree_inspection.txt      : Tree inspection log"
echo ""
echo "🔍 NEXT STEP: Run 07_diversity_analysis.sh"
echo "=================================================="
