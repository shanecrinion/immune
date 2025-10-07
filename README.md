## Transcriptomic Analysis of Immune Dysregulation in Schizophrenia 

### Overview

This repository documents the analytical framework for studying immune dysregulation in schizophrenia (SZ) using transcriptomic data, as presented in Chapter 3 of my PhD thesis. The study examines differential gene expression and immune cell-type composition in SZ patients.

Poster presented at the World Congress of Psychiatric Genetics 2024 in Singapore: https://docs.google.com/presentation/d/15D2L4XgWHvQBBKx9b9_kwarNQo5bRBKDXp3PJaWysyA/edit?usp=drive_link

### Workflow

#### Step 1: RNA-Seq Data Processing

- **Quality Control (QC):** Checking read quality with FastQC and MultiQC.

- **Read Alignment:** Mapping reads to the reference genome using HISAT2.

- **Normalization & Batch Correction:** Correcting for technical biases using DESeq2 and edgeR.

#### Step 2: Differential Gene Expression (DGE) Analysis

- **Case-Control Comparison:** Identifying up/downregulated genes.

- **Statistical Significance Testing:** Adjusting for multiple comparisons (FDR correction).

#### Step 3: Functional & Pathway Analysis

- **Gene Ontology (GO) & Pathway Enrichment:** Investigating biological implications.

- **Cell-Type Deconvolution:**  Estimating immune cell proportions within bulk RNA-seq samples.

- **Visualization:** Generating volcano plots, heatmaps, and PCA plots.

#### Required Software & Dependencies

- FastQC, MultiQC (QC)

- HISAT2 (Read alignment)

- DESeq2, edgeR (Differential expression analysis)

- CIBERSORT (Cellular deconvolution)

- R & Bioconductor (Statistical analysis, visualization)


