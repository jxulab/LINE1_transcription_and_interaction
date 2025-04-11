# Analysis of LINE1 chrRNA-seq and L1CaP-C Data

This repository contains the code used to analyze LINE1 chrRNA-seq and L1CaP-C data for the study:

> *"LINE-1 locus transcription nucleates oncogenic chromatin architecture"*

The analyses are implemented in R and Bash and are organized into seven main scripts, each targeting a specific aspect of the study.

---

## 📁 Repository Structure

### 1. `Clustering_LINE1_from_chrRNA-seq.R`
- Performs clustering of all full-length LINE1 elements across 10 cell lines.
- Generates heatmaps for all full-length LINE1s or specific LINE1 subfamilies.

### 2. `Identify_expression_threshold.R`
- Defines a TPM threshold to distinguish transcriptional noise from true expression. LINE1 loci with TPM above the threshold are considered expressed.

### 3. `motif_analysis.R`
- Performs TP motif enrichment analysis on the 5'UTR promoter sequences of expressed vs. non-expressed LINE1 loci across ten cell lines.

### 4. `Calculate_significance_of_interactions.R`
- Custom pipeline for identifying significant intra- and inter-chromosomal interactions.

### 5. `HILL_analysis.R`
- Identifies Highly Interacting LINE1 Loci (HILLs).
- Generates boxplots comparing expression of HILLs and HILL-interacting genes.

### 6. `HMM_feature_analysis.R`
- Identifies ChromHMM features enriched at expressed vs. non-expressed LINE1 loci.
- Evaluates enrichment of ChromHMM features at LINE1-interacting regions.

### 7. `Chromunity_analysis.R`
- Performs higher-order synergy (Chromunity) analysis across three cell lines and LINE-1 ASO knockdown data.
- Compares HILLs with detected 'Chromunity' events.

---

## 📦 Dependencies

Make sure the following R packages are installed before running the scripts:

bedtoolsr,dplyr,data.table,corrplot, Hmisc, cowplot, gplots, RColorBrewer, ggplot2, stringr,GenomicRanges,rtracklayer,skitools,chromunity,MASS,magrittr,gUtils,readr, circlize
