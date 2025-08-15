# Tick Population Genomics Pipeline

**Population genomics analysis of *Ixodes scapularis* across the Midwest**

## Overview

This repository contains a comprehensive 17-step bioinformatics pipeline for analyzing population structure in blacklegged ticks (*Ixodes scapularis*) using whole-genome sequencing data. The pipeline processes 61 tick samples from Nebraska, Iowa, and Kansas through quality control, variant calling, filtering, and population genetic analysis.

## Pipeline Workflow

### 🔵 Preprocessing (Steps 1-8)
1. **Concatenate Reads** - Combine lane files
2. **Fastp QC** - Quality control and trimming  
3. **FastQC** - Quality assessment
4. **BWA Alignment** - Map to reference genome
5. **SAM to BAM** - Convert and sort alignments
6. **Add Read Groups** - GATK metadata preparation
7. **Remove Duplicates** - PCR artifact removal
8. **Summary Statistics** - Comprehensive QC report

### 🟢 Variant Calling (Steps 9-11)
9. **HaplotypeCaller** - Individual variant discovery
10. **GenomicsDB** - Consolidate variant data
11. **Joint Genotyping** - Population variant calling

### 🟠 Filtering (Steps 12-15)
12. **Select SNPs** - Extract SNPs only
13. **Hard Filtering** - Apply GATK quality filters
14. **Select PASS** - Keep high-quality variants
15. **Population Filters** - MAF and missingness filtering

### 🟣 Analysis (Steps 16-17)
16. **PLINK + PCA** - Population structure analysis
17. **Visualization** - Create publication plots

## Key Features

- **Comprehensive Documentation**: Detailed methodology for each pipeline step
- **GATK Best Practices**: Rigorous variant calling and quality control
- **Population Optimization**: MAF ≥5% and missingness ≤30% filtering
- **Scalable Processing**: Optimized for HPC cluster environments
- **Publication Ready**: High-quality visualizations and analysis

## Dataset

- **Samples**: 61 *Ixodes scapularis* ticks
- **Regions**: Nebraska, Iowa, Kansas  
- **Sequencing**: NovaSeq whole-genome sequencing
- **Data Volume**: >1TB genomic data

## Technical Specifications

- **Platform**: SLURM-based HPC cluster
- **Key Tools**: GATK4, BWA-MEM, PLINK, VCFtools
- **Languages**: Bash scripting, Python visualization
- **Memory**: Optimized for large-scale population genomics

## Repository Contents

- `scripts/` - All 17 pipeline scripts
- `pipeline-documentation.md` - Comprehensive step-by-step analysis
- `pca_generation.py` - Population structure visualization
- `README.md` - This overview

## Quick Start

1. Ensure all dependencies are installed (GATK4, BWA, PLINK, etc.)
2. Modify file paths in scripts for your environment
3. Run scripts sequentially (Steps 1-17)
4. See `pipeline-documentation.md` for detailed guidance


---

**Input:** Raw FASTQ files → **Output:** Population structure insights
