# QIIME2 16S rRNA Analysis Pipeline

A complete QIIME2 workflow for 16S rRNA sequencing data analysis of BPD (Bronchopulmonary Dysplasia) microbiome data.

## 📁 Project Structure

```
qiime2_analysis/
├── run_qiime2_pipeline.sh          # 🚀 Complete analysis pipeline
├── config.yaml                     # ⚙️ Analysis parameters
├── manifest.csv                    # 📋 Sample manifest template
├── metadata.txt                    # 📊 Metadata template
├── scripts/                        # 🔧 Individual analysis steps
├── database/                       # 🗄️ Reference databases
```

## 🚀 Quick Start

### 1. Setup Environment
```bash
# Activate QIIME2 environment
conda activate qiime2-amplicon-2023.9

# Download SILVA database
chmod +x database/download_silva.sh
./database/download_silva.sh
```

### 2. Prepare Input Files
- `manifest.csv` - Sample sequencing files
- `metadata.txt` - Sample metadata

### 3. Run Complete Pipeline
```bash
chmod +x run_qiime2_pipeline.sh
./run_qiime2_pipeline.sh
```

### 4. Run Individual Steps
```bash
# Import data
./scripts/01_import_data.sh

# Trim primers
./scripts/02_trim_primers.sh

# DADA2 denoising
./scripts/03_dada2_denoise.sh

# Taxonomic classification
./scripts/04_taxonomic_classification.sh

# Data filtering
./scripts/05_filter_data.sh

# Phylogenetic analysis
./scripts/06_phylogenetic_tree.sh

# Diversity analysis
./scripts/07_diversity_analysis.sh

# Export results
./scripts/08_export_data.sh
```

## 📊 Input Requirements

### Manifest File (`manifest.csv`)
```csv
sample-id,absolute-filepath,direction
sample1,/path/to/sample1_R1.fastq.gz,forward
sample1,/path/to/sample1_R2.fastq.gz,reverse
```

### Metadata File (`metadata.txt`)
```txt
#SampleID	BPD_Group	Gestational_Age	Birth_Weight
sample1	Mild	28	1200
sample2	Severe	26	800
```

## 🔧 Analysis Steps

1. **Data Import** - Import paired-end sequences
2. **Primer Trimming** - Remove V3-V4 primers (338F/806R)
3. **DADA2 Denoising** - Quality filtering and ASV calling
4. **Taxonomic Classification** - SILVA database annotation
5. **Data Filtering** - Remove mitochondria/chloroplast
6. **Phylogenetic Tree** - Phylogenetic diversity analysis
7. **Diversity Analysis** - Alpha and beta diversity
8. **Data Export** - Export tables for downstream analysis

## ⚙️ Configuration

Edit `config.yaml` to adjust analysis parameters:

```yaml
primers:
  forward: "ACTCCTACGGGAGGCAGCAG"
  reverse: "GGACTACHVGGGTWTCTAAT"
  
dada2:
  trunc_len_f: 270
  trunc_len_r: 218
  
filtering:
  min_frequency: 1
  sampling_depth: 22000
```

## 📈 Output Files

- **Feature Tables**: ASV abundance tables
- **Taxonomy**: Taxonomic classifications
- **Phylogenetic Trees**: For diversity analysis
- **Diversity Metrics**: Alpha and beta diversity results
- **Visualizations**: Interactive QIIME2 views

## 🗄️ Database Requirements

- SILVA 138 database for 16S rRNA classification
- V3-V4 region specific classifier

## 🔍 Quality Control

- Sequence quality plots
- Denoising statistics
- Feature table summaries
- Rarefaction curves

## 📚 References

- QIIME2 Documentation: https://docs.qiime2.org
- SILVA Database: https://www.arb-silva.de
- Bolyen E, et al. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology. 2019.
