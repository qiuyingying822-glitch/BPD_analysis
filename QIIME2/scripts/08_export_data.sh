#!/bin/bash

# QIIME2 Data Export Script
# Step 8: Export results for downstream analysis

set -e  # Exit on error

echo "=================================================="
echo "QIIME2 STEP 8: DATA EXPORT"
echo "Started: $(date)"
echo "=================================================="

# Check if input files exist
if [[ ! -f "table-filtered-taxa.qza" ]]; then
    echo "❌ ERROR: table-filtered-taxa.qza not found!"
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

echo "📤 Exporting QIIME2 results for downstream analysis..."

# Create export directories
mkdir -p export/feature_table
mkdir -p export/rep_seqs
mkdir -p export/taxonomy
mkdir -p export/diversity
mkdir -p export/rarefied

echo "📊 Step 8.1: Export feature table and taxonomy..."

# Export main feature table (non-rarefied)
qiime tools export \
  --input-path table-filtered-taxa.qza \
  --output-path export/feature_table

# Export representative sequences
qiime tools export \
  --input-path rep-seqs-filtered.qza \
  --output-path export/rep_seqs

# Export taxonomy
qiime tools export \
  --input-path taxonomy-filtered.qza \
  --output-path export/taxonomy

echo "🔧 Step 8.2: Process taxonomy file..."

# Process taxonomy file header
sed -i -e '1 s/Feature/#Feature/' -e '1 s/Taxon/taxonomy/' export/taxonomy/taxonomy.tsv

echo "🔧 Step 8.3: Add taxonomy to BIOM file..."

# Add taxonomy to BIOM file
biom add-metadata \
  -i export/feature_table/feature-table.biom \
  -o export/feature_table/feature_table_taxonomy.biom \
  --observation-metadata-fp export/taxonomy/taxonomy.tsv \
  --sc-separated taxonomy

echo "🔧 Step 8.4: Convert BIOM to TSV..."

# Convert BIOM to TSV with taxonomy
biom convert \
  -i export/feature_table/feature_table_taxonomy.biom \
  -o export/feature_table/feature_table_taxonomy.txt \
  --to-tsv \
  --header-key taxonomy

echo "📊 Step 8.5: Generate rarefied table..."

# Create rarefied table for specific analyses
qiime feature-table rarefy \
  --i-table table-filtered-taxa.qza \
  --p-sampling-depth 22000 \
  --o-rarefied-table export/rarefied/table-rarefied.qza

# Export rarefied table
qiime tools export \
  --input-path export/rarefied/table-rarefied.qza \
  --output-path export/rarefied

echo "🔧 Step 8.6: Filter sequences for rarefied table..."

# Filter sequences to match rarefied table
qiime feature-table filter-seqs \
  --i-data rep-seqs-filtered.qza \
  --i-table export/rarefied/table-rarefied.qza \
  --o-filtered-data export/rarefied/rep-seqs-rarefied.qza

echo "🔧 Step 8.7: Classify rarefied sequences..."

# Classify rarefied sequences
qiime feature-classifier classify-sklearn \
  --i-classifier database/silva-138-99-v3v4-classifier.qza \
  --i-reads export/rarefied/rep-seqs-rarefied.qza \
  --o-classification export/rarefied/taxonomy-rarefied.qza

# Export rarefied taxonomy
qiime tools export \
  --input-path export/rarefied/taxonomy-rarefied.qza \
  --output-path export/rarefied/taxonomy

echo "🔧 Step 8.8: Process rarefied taxonomy..."

# Process rarefied taxonomy file
sed -i -e '1 s/Feature/#Feature/' -e '1 s/Taxon/taxonomy/' export/rarefied/taxonomy/taxonomy.tsv

# Add taxonomy to rarefied BIOM
biom add-metadata \
  -i export/rarefied/feature-table.biom \
  -o export/rarefied/feature_table_taxonomy.biom \
  --observation-metadata-fp export/rarefied/taxonomy/taxonomy.tsv \
  --sc-separated taxonomy

# Convert rarefied BIOM to TSV
biom convert \
  -i export/rarefied/feature_table_taxonomy.biom \
  -o export/rarefied/feature_table_taxonomy.txt \
  --to-tsv \
  --header-key taxonomy

echo "📊 Step 8.9: Export diversity results..."

# Export diversity results if they exist
if [[ -d "core-metrics-results" ]]; then
    echo "📈 Exporting diversity metrics..."
    
    # Export alpha diversity vectors
    for file in core-metrics-results/*_vector.qza; do
        if [[ -f "$file" ]]; then
            base=$(basename "$file" .qza)
            qiime tools export \
              --input-path "$file" \
              --output-path "export/diversity/${base}"
        fi
    done
    
    # Export beta diversity matrices
    for file in core-metrics-results/*_distance_matrix.qza; do
        if [[ -f "$file" ]]; then
            base=$(basename "$file" .qza)
            qiime tools export \
              --input-path "$file" \
              --output-path "export/diversity/${base}"
        fi
    done
fi

echo "📈 Export Summary:"

# Display export summary
echo "📁 Exported Files Structure:"
find export/ -type f -name "*.tsv" -o -name "*.txt" -o -name "*.biom" | sort

echo "🔢 File Statistics:"
echo "  - Features in main table: $(wc -l < export/feature_table/feature_table_taxonomy.txt | awk '{print $1-1}')"
echo "  - Features in rarefied table: $(wc -l < export/rarefied/feature_table_taxonomy.txt | awk '{print $1-1}')"
echo "  - Samples in main table: $(head -1 export/feature_table/feature_table_taxonomy.txt | wc -w)"
echo "  - Samples in rarefied table: $(head -1 export/rarefied/feature_table_taxonomy.txt | wc -w)"

echo "📝 Logging export details..."

# Create export manifest
echo "QIIME2 Data Export Manifest" > logs/export_manifest.txt
echo "Generated: $(date)" >> logs/export_manifest.txt
echo "============================" >> logs/export_manifest.txt
find export/ -type f >> logs/export_manifest.txt

echo "=================================================="
echo "DATA EXPORT COMPLETED SUCCESSFULLY"
echo "Finished: $(date)"
echo "=================================================="
echo ""
echo "📋 EXPORTED DATA:"
echo "  - export/feature_table/     : Main feature tables with taxonomy"
echo "  - export/rep_seqs/          : Representative sequences (FASTA)"
echo "  - export/taxonomy/          : Taxonomy assignments"
echo "  - export/rarefied/          : Rarefied data for diversity analysis"
echo "  - export/diversity/         : Diversity metrics (if available)"
echo ""
echo "🎉 QIIME2 ANALYSIS PIPELINE COMPLETED!"
echo "Next: Use exported data for SparCC network analysis"
echo "Run: ../sparcc_analysis/run_sparcc_pipeline.sh"
echo "=================================================="
