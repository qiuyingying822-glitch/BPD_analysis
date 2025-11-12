#!/bin/bash

# QIIME2 DADA2 Denoising Script
# Step 3: Quality filtering, denoising, and ASV calling

set -e  # Exit on error

echo "=================================================="
echo "QIIME2 STEP 3: DADA2 DENOISING"
echo "Started: $(date)"
echo "=================================================="

# Check if input file exists
if [[ ! -f "16s_paired_end_primer.qza" ]]; then
    echo "❌ ERROR: 16s_paired_end_primer.qza not found!"
    echo "Please run 02_trim_primers.sh first"
    exit 1
fi

# Check if QIIME2 is available
if ! command -v qiime &> /dev/null; then
    echo "❌ ERROR: QIIME2 not found in PATH"
    echo "Please activate QIIME2 environment:"
    echo "conda activate qiime2-amplicon-2023.9"
    exit 1
fi

# DADA2 parameters from config or defaults
TRUNC_LEN_F=270    # Forward read truncation length
TRUNC_LEN_R=218    # Reverse read truncation length
MAX_EE_F=2.0       # Forward read max expected errors
MAX_EE_R=2.0       # Reverse read max expected errors
TRUNC_Q=2          # Truncation quality score

echo "⚙️  DADA2 Parameters:"
echo "  - Forward truncation: $TRUNC_LEN_F bp"
echo "  - Reverse truncation: $TRUNC_LEN_R bp"
echo "  - Max expected errors (F/R): $MAX_EE_F / $MAX_EE_R"
echo "  - Truncation quality: $TRUNC_Q"
echo "  - Overlap length: $((TRUNC_LEN_F + TRUNC_LEN_R - 468)) bp"

echo "🔬 Starting DADA2 denoising (this may take a while)..."

# Run DADA2 denoising with timing
time qiime dada2 denoise-paired \
  --i-demultiplexed-seqs 16s_paired_end_primer.qza \
  --p-trunc-len-f $TRUNC_LEN_F \
  --p-trunc-len-r $TRUNC_LEN_R \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-max-ee-f $MAX_EE_F \
  --p-max-ee-r $MAX_EE_R \
  --p-trunc-q $TRUNC_Q \
  --p-chimera-method consensus \
  --p-n-threads 8 \
  --p-n-reads-learn 1000000 \
  --p-hashed-feature-ids True \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats denoising-stats.qza

if [[ $? -eq 0 ]]; then
    echo "✅ DADA2 denoising completed successfully"
else
    echo "❌ DADA2 denoising failed"
    exit 1
fi

echo "📊 Creating DADA2 result visualizations..."

# Create feature table summary
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization visual/table.qzv

# Create representative sequences summary
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization visual/rep-seqs.qzv

# Create denoising statistics summary
qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization visual/denoising-stats.qzv

echo "📈 DADA2 Statistics Summary:"

# Extract basic statistics
echo "🔍 Feature table information:"
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization visual/table-summary.qzv 2>&1 | grep -E "(Sample count|Feature count|Total frequency)"

# Log denoising results
echo "📝 Logging denoising details..."
qiime tools inspect table.qza > logs/table_inspection.txt
qiime tools inspect rep-seqs.qza > logs/rep_seqs_inspection.txt

echo "=================================================="
echo "DADA2 DENOISING COMPLETED SUCCESSFULLY"
echo "Finished: $(date)"
echo "=================================================="
echo ""
echo "📋 OUTPUT FILES:"
echo "  - table.qza                      : Feature table (ASV abundances)"
echo "  - rep-seqs.qza                   : Representative sequences (ASVs)"
echo "  - denoising-stats.qza            : Denoising statistics"
echo "  - visual/table.qzv               : Feature table visualization"
echo "  - visual/rep-seqs.qzv            : Representative sequences visualization"
echo "  - visual/denoising-stats.qzv     : Denoising stats visualization"
echo ""
echo "🔍 NEXT STEP: Run 04_taxonomic_classification.sh"
echo "=================================================="
