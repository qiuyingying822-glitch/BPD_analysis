#!/bin/bash

# QIIME2 Diversity Analysis Script
# Step 7: Alpha and beta diversity analysis

set -e  # Exit on error

echo "=================================================="
echo "QIIME2 STEP 7: DIVERSITY ANALYSIS"
echo "Started: $(date)"
echo "=================================================="

# Check if input files exist
if [[ ! -f "table-filtered-taxa.qza" ]]; then
    echo "❌ ERROR: table-filtered-taxa.qza not found!"
    echo "Please run 05_filter_data.sh first"
    exit 1
fi

if [[ ! -f "rooted-tree.qza" ]]; then
    echo "❌ ERROR: rooted-tree.qza not found!"
    echo "Please run 06_phylogenetic_tree.sh first"
    exit 1
fi

if [[ ! -f "metadata.txt" ]]; then
    echo "❌ ERROR: metadata.txt not found!"
    echo "Please create metadata file for diversity analysis"
    exit 1
fi

# Check if QIIME2 is available
if ! command -v qiime &> /dev/null; then
    echo "❌ ERROR: QIIME2 not found in PATH"
    echo "Please activate QIIME2 environment:"
    echo "conda activate qiime2-amplicon-2023.9"
    exit 1
fi

# Diversity analysis parameters
SAMPLING_DEPTH=22000

echo "⚙️  Diversity Analysis Parameters:"
echo "  - Sampling depth: $SAMPLING_DEPTH"
echo "  - Phylogenetic tree: rooted-tree.qza"
echo "  - Metadata: metadata.txt"

echo "📊 Step 7.1: Core diversity metrics analysis..."

# Core diversity metrics (alpha and beta diversity)
time qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table-filtered-taxa.qza \
  --p-sampling-depth $SAMPLING_DEPTH \
  --m-metadata-file metadata.txt \
  --p-n-jobs 8 \
  --output-dir core-metrics-results

if [[ $? -eq 0 ]]; then
    echo "✅ Core diversity metrics completed"
else
    echo "❌ Core diversity metrics failed"
    echo "⚠️  Trying with lower sampling depth..."
    
    # Try with automatic sampling depth
    qiime diversity core-metrics-phylogenetic \
      --i-phylogeny rooted-tree.qza \
      --i-table table-filtered-taxa.qza \
      --p-sampling-depth 1000 \
      --m-metadata-file metadata.txt \
      --p-n-jobs 8 \
      --output-dir core-metrics-results
fi

echo "📈 Step 7.2: Alpha diversity analysis..."

# Alpha rarefaction analysis
echo "📊 Generating alpha rarefaction curves..."

qiime diversity alpha-rarefaction \
  --i-table table-filtered-taxa.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth $SAMPLING_DEPTH \
  --m-metadata-file metadata.txt \
  --o-visualization visual/alpha-rarefaction.qzv

echo "📊 Step 7.3: Alpha group significance testing..."

# Alpha diversity group significance
for metric in observed_features shannon faith_pd; do
    echo "🔍 Testing group significance for: $metric"
    
    qiime diversity alpha-group-significance \
      --i-alpha-diversity core-metrics-results/${metric}_vector.qza \
      --m-metadata-file metadata.txt \
      --o-visualization visual/alpha-group-significance-${metric}.qzv
done

echo "📊 Step 7.4: Beta diversity group significance..."

# Beta diversity group significance
for metric in bray_curtis weighted_unifrac unweighted_unifrac; do
    echo "🔍 Testing group significance for: $metric"
    
    qiime diversity beta-group-significance \
      --i-distance-matrix core-metrics-results/${metric}_distance_matrix.qza \
      --m-metadata-file metadata.txt \
      --m-metadata-column BPD_Group \
      --p-pairwise \
      --o-visualization visual/beta-group-significance-${metric}.qzv
done

echo "📊 Step 7.5: PCoA visualization..."

# Create Emperor PCoA plots
for metric in bray_curtis weighted_unifrac unweighted_unifrac; do
    echo "🎨 Creating PCoA plot for: $metric"
    
    qiime emperor plot \
      --i-pcoa core-metrics-results/${metric}_pcoa_results.qza \
      --m-metadata-file metadata.txt \
      --o-visualization visual/emperor-${metric}.qzv
done

echo "📈 Diversity Analysis Summary:"

# Display diversity metrics summary
echo "🔍 Core Metrics Generated:"
ls core-metrics-results/*.qza | xargs -n 1 basename

echo "📊 Alpha Diversity Metrics:"
echo "  - observed_features: Species richness"
echo "  - shannon: Shannon diversity index"
echo "  - faith_pd: Faith's phylogenetic diversity"
echo "  - pielou_e: Pielou's evenness"

echo "📊 Beta Diversity Metrics:"
echo "  - bray_curtis: Bray-Curtis dissimilarity"
echo "  - weighted_unifrac: Weighted UniFrac"
echo "  - unweighted_unifrac: Unweighted UniFrac"

echo "📝 Logging diversity analysis details..."

# Log diversity results
qiime tools inspect core-metrics-results/ > logs/diversity_inspection.txt

echo "=================================================="
echo "DIVERSITY ANALYSIS COMPLETED SUCCESSFULLY"
echo "Finished: $(date)"
echo "=================================================="
echo ""
echo "📋 OUTPUT FILES:"
echo "  - core-metrics-results/          : Core diversity metrics"
echo "  - visual/alpha-rarefaction.qzv   : Alpha rarefaction curves"
echo "  - visual/alpha-group-significance-*.qzv : Alpha diversity tests"
echo "  - visual/beta-group-significance-*.qzv  : Beta diversity tests"
echo "  - visual/emperor-*.qzv           : PCoA plots"
echo "  - logs/diversity_inspection.txt  : Diversity analysis log"
echo ""
echo "🔍 NEXT STEP: Run 08_export_data.sh"
echo "=================================================="
