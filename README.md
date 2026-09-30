# Salivary Microbiome Smoking Analysis

## Overview

This repository contains the R-based analysis workflow used for the
forensic assessment of smoking-associated salivary microbiome profiles.

The workflow includes taxonomic profiling, prevalence analysis,
microbial relative abundance, diversity analysis, ordination,
statistical testing, differential abundance analysis, and visualization.

## Study Focus

The analysis investigates bacterial community profiles in saliva in
relation to smoking status.

The workflow includes comparisons among:

- Current smokers
- Former smokers
- Never smokers

## Analysis Workflow

The R script includes the following major analytical steps:

1. Import and processing of Kraken2/Bracken taxonomic results
2. Construction of microbial abundance tables
3. Taxonomic profiling
4. Relative abundance analysis
5. Prevalence analysis
6. Alpha diversity analysis
7. Beta diversity analysis
8. PCoA ordination
9. PERMANOVA
10. Differential abundance analysis
11. Visualization of microbial patterns

## Main R Packages

The analysis uses R packages including:

- `phyloseq`
- `vegan`
- `DESeq2`
- `ggplot2`
- `dplyr`
- `tidyr`

## Repository Contents

| File | Description |
|------|-------------|
| `Kraken2_R_analysis.R` | Complete R analysis workflow for the study |

## Data

Raw sequencing data and participant-level metadata are **not included**
in this repository.

Only the R analysis workflow is provided for reproducibility and
demonstration purposes.

## Reproducibility

The script is provided as a complete analytical workflow. Users wishing
to reproduce the analysis should provide their own appropriately
formatted Kraken2/Bracken results and metadata.

## Research Context

This work forms part of research investigating the potential application
of salivary microbiome profiles in forensic identification.

## Author

**Sreekutti S**

National Forensic Sciences University  
Gandhinagar, Gujarat, India
