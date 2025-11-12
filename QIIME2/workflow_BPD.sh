#!/bin/bash
# QIIME2 workflow for BPD 16S analysis
# Version: 2023.9
# Author: [Your Name]
# ------------------------------------------

# 设置工作目录
wd=~/BPD
cd $wd
conda activate qiime2-amplicon-2023.9

# 1. Import data
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.csv \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33

# 2. Primer trimming
qiime cutadapt trim-paired \
  --p-cores 8 \
  --i-demultiplexed-sequences paired-end-demux.qza \
  --p-front-f ACTCCTACGGGAGGCAGCAG \
  --p-front-r GGACTACHVGGGTWTCTAAT \
  --o-trimmed-sequences 16s_paired_end_primer.qza

# 3. Visualize sequence quality
qiime demux summarize \
  --i-data 16s_paired_end_primer.qza \
  --o-visualization 16s_paired_end_primer.qzv

# 4. DADA2 denoising
time qiime dada2 denoise-paired \
  --i-demultiplexed-seqs 16s_paired_end_primer.qza \
  --p-trunc-len-f 270 \
  --p-trunc-len-r 218 \
  --o-representative-sequences representative-sequences.qza \
  --o-table table.qza \
  --o-denoising-stats denoising-stats.qza

# 5. Taxonomic classification (SILVA 138)
qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer ACTCCTACGGGAGGCAGCAG \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --o-reads 99-v3v4-seqs.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads 99-v3v4-seqs.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
  --o-classifier silva-138-99-v3v4-classifier.qza

qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-v3v4-classifier.qza \
  --i-reads representative-sequences.qza \
  --o-classification taxonomy_train_silva.qza

# 6. Filter and rarefaction
qiime feature-table filter-features \
  --i-table table.qza \
  --p-min-frequency 1 \
  --o-filtered-table table_filter_freq_1.qza

qiime taxa filter-table \
  --i-table table_filter_freq_1.qza \
  --i-taxonomy taxonomy_train_silva.qza \
  --p-include "d__Bacteria; p_" \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table table-with-bacteria.qza

qiime feature-table rarefy \
  --i-table table-with-bacteria.qza \
  --p-sampling-depth 22000 \
  --o-rarefied-table ./rarefy/table-rarefy-22000.qza

# 7. Phylogenetic tree
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences representative-sequences.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

# 8. Alpha rarefaction
qiime diversity alpha-rarefaction \
  --i-table table-with-bacteria.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 30000 \
  --m-metadata-file metadata.txt \
  --o-visualization alpha-rarefaction.qzv

# 9. Export results
qiime tools export --input-path table-with-bacteria.qza --output-path feature_table
qiime tools export --input-path representative-sequences.qza --output-path rep_seqs
qiime tools export --input-path taxonomy_train_silva.qza --output-path taxonomy
