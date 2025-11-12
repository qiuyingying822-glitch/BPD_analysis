#!/bin/bash

# QIIME2 16S Analysis Pipeline for BPD Microbiome Data
# Complete workflow from raw sequences to analyzed data

set -e  # Exit on error

echo "=== Starting QIIME2 16S Analysis Pipeline ==="

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/config.yaml"

# Load configuration
if [[ -f "$CONFIG_FILE" ]]; then
    echo "Loading configuration from: $CONFIG_FILE"
    # Parse config (simplified)
    FORWARD_PRIMER=$(grep "forward:" "$CONFIG_FILE" | cut -d' ' -f2)
    REVERSE_PRIMER=$(grep "reverse:" "$CONFIG_FILE" | cut -d' ' -f2)
    TRUNC_LEN_F=$(grep "trunc_len_f:" "$CONFIG_FILE" | cut -d' ' -f2)
    TRUNC_LEN_R=$(grep "trunc_len_r:" "$CONFIG_FILE" | cut -d' ' -f2)
    SAMPLING_DEPTH=$(grep "sampling_depth:" "$CONFIG_FILE" | cut -d' ' -f2)
else
    echo "Using default parameters"
    FORWARD_PRIMER="ACTCCTACGGGAGGCAGCAG"
    REVERSE_PRIMER="GGACTACHVGGGTWTCTAAT"
    TRUNC_LEN_F=270
    TRUNC_LEN_R=218
    SAMPLING_DEPTH=22000
fi

# Create output directories
mkdir -p logs
mkdir -p visual
mkdir -p rarefy
mkdir -p feature_table
mkdir -p rep_seqs
mkdir -p taxonomy

echo "1. Activating QIIME2 environment..."
# conda activate qiime2-amplicon-2023.9

echo "2. Running analysis steps..."

# Step 1: Import data
echo "Step 1: Importing sequencing data..."
./scripts/01_import_data.sh

# Step 2: Trim primers
echo "Step 2: Trimming primers..."
./scripts/02_trim_primers.sh

# Step 3: DADA2 denoising
echo "Step 3: DADA2 denoising..."
./scripts/03_dada2_denoise.sh

# Step 4: Taxonomic classification
echo "Step 4: Taxonomic classification..."
./scripts/04_taxonomic_classification.sh

# Step 5: Data filtering
echo "Step 5: Filtering data..."
./scripts/05_filter_data.sh

# Step 6: Phylogenetic tree
echo "Step 6: Generating phylogenetic tree..."
./scripts/06_phylogenetic_tree.sh

# Step 7: Diversity analysis
echo "Step 7: Diversity analysis..."
./scripts/07_diversity_analysis.sh

# Step 8: Export data
echo "Step 8: Exporting results..."
./scripts/08_export_data.sh

echo "=== QIIME2 Analysis Pipeline Completed ==="
echo "Results available in:"
echo "  - visual/          : QIIME2 visualizations"
echo "  - feature_table/   : Feature tables"
echo "  - taxonomy/        : Taxonomy assignments"
echo "  - rarefy/          : Rarefied data"
echo "  - logs/            : Analysis logs"
