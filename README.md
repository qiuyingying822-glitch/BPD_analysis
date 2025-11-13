# BPD Microbiome Network Analysis

A comprehensive computational pipeline for analyzing the relationship between gut microbial networks and Bronchopulmonary Dysplasia (BPD) severity in preterm infants.

## 📋 Project Overview

This repository contains a complete analytical workflow for studying the gut microbiome's role in BPD development through network analysis approaches.

## Project Structure
BPD_analysis/
├── QIIME2/
│   ├── README.md
│   ├── run_qiime2_pipeline.sh
│   ├── scripts/
│   ├── config.yaml
│   └── manifest.csv
├── network_analysis/
│   ├── README.md
│   ├── sparcc_analysis/
│   │   ├── run_sparcc_pipeline.sh
│   │   ├── process_sparcc_results.R
│   │   ├── install_sparcc.sh
│   │   └── config.yaml
│   ├── keystone_analysis/
│   │   ├── identify_keystone_taxa_enhanced.R
│   │   ├── identify_keystone_taxa_sensitivity.R
│   │   └── export_gephi_files.R
│   └── sensitivity_analysis/
│       ├── analyze_low_correlation.R
│       └── method_comparison.R
└── regression_analysis/
    ├── README.md
    ├── network_metrics_model/
    ├── keystone_multivariate/
    └── README.md

## 🔄 Analysis Workflow

### Phase 1: 16S Data Processing (`qiime2_analysis/`)
- Raw sequence processing and quality control
- ASV identification and taxonomic classification
- Diversity analysis and data export

### Phase 2: Microbial Network Analysis (`network_analysis/`)
- **Network Construction** (`sparcc_analysis/`): Build co-occurrence networks using SparCC
- **Topology Analysis** (`keystone_analysis/`): Identify keystone taxa and network properties
- **Method Validation** (`sensitivity_analysis/`): Compare methods and test robustness

### Phase 3: Statistical Modeling (`regression_analysis/`)
- Predict BPD severity using network metrics
- Analyze keystone taxa associations with clinical outcomes

## 🚀 Quick Start

### 1. Sequential Analysis
```bash
# 1. Process 16S data
cd qiime2_analysis
./run_qiime2_pipeline.sh

# 2. Build and analyze networks
cd ../network_analysis/sparcc_analysis
./run_sparcc_pipeline.sh

cd ../keystone_analysis
Rscript identify_keystone_taxa_enhanced.R

# 3. Statistical modeling
cd ../../regression_analysis/network_metrics_model
Rscript network_metrics_regression.R

###2. Individual Module Usage
#Each module can be run independently - see the respective README files for detailed instructions.

##📁 Module Details
##🔬 QIIME2 Analysis
Input: Raw FASTQ sequences

Output: ASV tables, taxonomic assignments

See qiime2_analysis/README.md

## 🌐 Network Analysis (Integrated)
SparCC Networks: Correlation-based network construction

Keystone Identification: Topological analysis and key species detection

Sensitivity Tests: Method comparisons and robustness validation

See network_analysis/README.md

## 📊 Statistical Modeling
Network metrics as predictors for BPD severity

Multivariate analysis of keystone taxa

See regression_analysis/README.md

## ⚙️ Configuration
Each module contains its own config.yaml for parameter customization:

Analysis parameters

File paths and formats

Statistical thresholds

Visualization settings

## 📊 Output
Standardized outputs across all modules:

Processed data tables

Statistical results

Publication-ready visualizations

Reproducibility logs

## 🔧 Dependencies
See individual module READMEs for specific requirements:

QIIME2 (2023.9+)

R (4.0.0+) with tidyverse, igraph, MASS

Python (3.6+) with SparCC, numpy, scipy

## 📚 Citation
Please cite the relevant methodological papers and this repository if used in your research.

## 🤝 Contributing
Issues and pull requests are welcome for improvements to the analysis pipeline.
