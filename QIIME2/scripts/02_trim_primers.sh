#!/bin/bash

# QIIME2 Primer Trimming Script
# Step 2: Remove V3-V4 primer sequences (338F/806R)

set -e  # Exit on error

echo "=================================================="
echo "QIIME2 STEP 2: PRIMER TRIMMING"
echo "Started: $(date)"
echo "=================================================="

# Check if input file exists
if [[ ! -f "paired-end-demux.qza" ]]; then
    echo "❌ ERROR: paired-end-demux.qza not found!"
    echo "Please run 01_import_data.sh first"
    exit 1
fi

# Check if QIIME2 is available
if ! command -v qiime &> /dev/null; then
    echo "❌ ERROR: QIIME2 not found in PATH"
    echo "Please activate QIIME2 environment:"
    echo "conda activate qiime2-amplicon-2023.9"
    exit 1
fi

# Primer sequences for V3-V4 region
FORWARD_PRIMER="ACTCCTACGGGAGGCAGCAG"  # 338F
REVERSE_PRIMER="GGACTACHVGGGTWTCTAAT"  # 806R

echo "🔬 Primer sequences:"
echo "  - Forward (338F): $FORWARD_PRIMER"
echo "  - Reverse (806R): $REVERSE_PRIMER"

echo "✂️  Trimming primer sequences..."

# Trim primers using cutadapt
qiime cutadapt trim-paired \
  --p-cores 8 \
  --i-demultiplexed-sequences paired-end-demux.qza \
  --p-front-f $FORWARD_PRIMER \
  --p-front-r $REVERSE_PRIMER \
  --p-error-rate 0.1 \
  --p-indels 1 \
  --p-minimum-length 100 \
  --o-trimmed-sequences 16s_paired_end_primer.qza

if [[ $? -eq 0 ]]; then
    echo "✅ Successfully trimmed primers: 16s_paired_end_primer.qza"
else
    echo "❌ Failed to trim primers"
    exit 1
fi

echo "📊 Creating primer trimming summary..."

# Create visualization of trimmed sequences
qiime demux summarize \
  --i-data 16s_paired_end_primer.qza \
  --o-visualization visual/16s_paired_end_primer.qzv

if [[ $? -eq 0 ]]; then
    echo "✅ Primer trimming summary created: visual/16s_paired_end_primer.qzv"
else
    echo "❌ Failed to create primer trimming summary"
    exit 1
fi

# Compare sequence counts before and after trimming
echo "📈 Sequence statistics:"
echo "  - Input file: paired-end-demux.qza"
echo "  - Output file: 16s_paired_end_primer.qza"

# Log trimming information
echo "📝 Logging trimming details..."
qiime tools inspect 16s_paired_end_primer.qza > logs/trimming_inspection.txt

echo "=================================================="
echo "PRIMER TRIMMING COMPLETED SUCCESSFULLY"
echo "Finished: $(date)"
echo "=================================================="
echo ""
echo "📋 OUTPUT FILES:"
echo "  - 16s_paired_end_primer.qza     : Primer-trimmed sequences"
echo "  - visual/16s_paired_end_primer.qzv : Trimming summary visualization"
echo "  - logs/trimming_inspection.txt  : Trimming inspection log"
echo ""
echo "🔍 NEXT STEP: Run 03_dada2_denoise.sh"
echo "=================================================="
