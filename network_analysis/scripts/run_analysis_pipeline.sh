#!/bin/bash

# BPD Network Analysis Pipeline
# This script runs the complete analysis pipeline

echo "=== Starting BPD Network Analysis Pipeline ==="

# Create necessary directories
mkdir -p network
mkdir -p sensitivity_analysis
mkdir -p gephi

# Run sensitivity analysis
echo "1. Running sensitivity analysis..."
Rscript network_analysis/sensitivity_analysis/analyze_low_correlation.R

# Run keystone taxa identification
echo "2. Identifying keystone taxa..."
Rscript network_analysis/keystone_analysis/identify_keystone_taxa_enhanced.R

# Run sensitivity analysis for keystone taxa
echo "3. Running keystone taxa sensitivity analysis..."
Rscript network_analysis/keystone_analysis/identify_keystone_taxa_sensitivity.R

echo "=== Analysis Pipeline Completed ==="
echo "Results saved in:"
echo "  - network/"
echo "  - sensitivity_analysis/"
echo "  - gephi/"
