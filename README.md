## RNA-seq Analysis to Investigate Immune Function and Stress Response in an Irish Schizophrenia Case-Control Study

### Overview

Immune factors are known to play a role in SZ, however their mechanisms remain largely unexplored. Blood based transcriptomic studies have valuable potential to investigate their complex interplay. This research, as part of the Immune Response & Social Cognition in Schizophrenia (iRELATE) project seeks to disentangle the complex genetic architecture of stress and schizophrenia (SZ) by investigating the presence of differentially expressed genes, enrichment for genes associated with disrupted cells, tissues and pathways, constructing new gene networks and their association with clinical and psychiatric phenotypes and examining cell type proportions

I received funding valued at €30,000 to generate RNA-sequencing data from the iRELATE study, which I used to complete this analysis. RNA-seq was performed for blood-based samples for 50 SZ cases and 50 healthy controls, before and after a stress test. The steps in this RNA-seq analysis workflow outlined in this repository include cloud computing set-up, sample QC, sequence alignment, post-alignment QC, read counting, differential expression analysis (DEA), gene enrichment analysis, weighted gene co-expression network analysis (WGCNA) and cellular deconvolution.

DEA identified 2,545 DEGs associated with SZ status and 770 DEGs associated with stress response. DEGs for SZ indicated immune responses, while DEGs for stress-test response were biologically aligned. Gene network construction identified one network was enriched for positive symptoms of SZ. One constructed module for SZ was found to be significantly enriched for CNS demyelination processes. Cell type proportions reveal a significant enrichment for T-cells CD8 associated with SZ. This study provides evidence for the relationship between immune response and SZ, and identifies biologically relevant genes, pathways and gene networks. This project also demonstrates that blood-based RNA-seq studies can lead to significant findings on the role of immune function in psychiatry and reveal the underpinning of its pathophysiology.

:speech_balloon: [Poster from World Congress of Psychiatric Genetics](media/WCPG_2024.pdf) — Singapore, 2024

### Project Description

This repository documents the analytical workflow for an RNA-seq analysis I developed to investigate immune dysregulation using blood-based transcriptomic data from the Immune Response and Social Cognition Study (iRELATE). This project was completed in part fulfillment of my PhD thesis titled '*Multiomic analysis of neuropsychiatric disorders*'. I completed a PhD in Genomics Data Science as a student of the SFI Centre for Research Training in Genomics Data Science.

#### Graphical Overview (experimental design and analytical workflow).

<div style="text-align: center; width="100%">
  <img src="media/experimental_design.png" alt="Overview of experimental design and analytical workflow." style="display: block; margin: 0 auto;">
</div>

**A.** *In the Patient Assessment Protocol, blood samples and neuropsychiatric tests are taken before and after the Trier Social Stress Test (TSST).* **B.** *In Sample Preparation & Library Generation, total RNA is isolated, polyadenylated, and sequenced using complementary DNA (cDNA). In Bioinformatics Analysis, DEA is used to identify differentially expressed genes and explore the biological pathophysiology linked to SZ and or stress response.*

#### RNA-Seq analysis workflow included exploratory analysis such as differential expression, gene-set enrichment and cellular deconvolution.

<div style="text-align: center; width="100%">
  <img src="media/rnaseq_workflow.png" alt="RNA-Seq analysis workflow." style="display: block; margin: 0 auto;">
</div>

*A number of bioinformatics tools were used in this project. Scripts on analysis methods are available in their respective subdirectories for each analysis step however post-alignment analysis were carried out primarily in R (using the tools DESeq2 for DEA, WGCNA for gene network construction and bMIND for cellular deconvolution).*

## Installation

-   Read quality was assessed before and after sequence alignment using FastQC (<https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>) and MultiQC (<https://seqera.io/multiqc/>).

-   Alignment was performed using HISAT2 and the pre-built human reference genome build 37 (<https://daehwankimlab.github.io/hisat2/>).

-   Post-alignment QC was performed using samtools with FastQC and MultiQC.

-   Read counting was measured using the Subread toolkit called featureCounts (<https://subread.sourceforge.net/featureCounts.html>).

The main R packages used for the RNA-seq analysis were DESeq2 (version 1.42.1), gage and WGCNA. These, and many other packages used in the individual analyses, are installed in R using BiocManager:

```         
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "gage", "WGCNA"))
```

## Procedure 

1.  Sequence Alignment | :file_folder: [seq_alignment/](seq_alignment)
    1.  Transfer files from external supplier to local server (<code>1_download_remaining_files.sh</code>).
    2.  Trim adaptor sequences (<code>2_trim_all.sh</code>).
    3.  Perform paired-end sequence alignment to GRCh37 (<code>3_run_hisat2.sh</code>).
    4.  Post-alignment QC using samtools including sorting, converting to BAM, indexing, remove PCR duplicates and generate QC metrics (<code>4_postprocessing.sh</code>).


2. Differential expression analysis | :file_folder: [diff_express_analysis/](diff_express_analysis)
    1. Format and merge phenotype and read count data for DESeq2, generate data object (<code>1_import_makeDDS.R</code>)
    2. QC plots to investigate the distribution and profile of read counts. (<code>Data_QC_plots.md</code>)
    3. Observe results of DEA for cases vs controls, create volcano plots and extract extract top genes (<code>2_cases_vs_controls/DEA_CasesVsControls.md</code>)
    4. Observe results of DEA for pre- (T1) and post-stress test (T3), create volcano plots and extract extract top genes (<code>3_tp1_vs_tp3/DEA_T1vsT3.md</code>)

3. Gene set analysis | :file_folder: [gene_set_analysis/](gene_set_analysis)
    1. GSA including library loading, Entrez-Ensembl ID mapping, pathway analysis using gage to investigate up and down regulated pathways in cases vs controls (<code>1_gene_set_analysis_cases_vs_controls.md</code>)
    2. Same as above for Timepoint 1 vs Timepoint 3 (<code>2_gene_set_analysis_T1vsT3.md</code>).
    
4. Weighted Gene Co-expression Network Analysis (WGCNA) |  :file_folder: [wgcna/](wgcna)
    1. Construction of gene networks to identify latent factors based on gene expression profiles. Script to construct networks, investigate association with case/control status, stress test time-point and interaction term between both. The significant modules are then extracted and subsets of the data are created to assess for association with 1) clinical measurements 2) psychiatric measurments (cases only) (<code>wgcna.md</code>).

5. Cellular deconvolution | :file_folder: [cellular_deconv/](cellular_deconv)
    1. This was performed to predict cell-type proportions on the blood-based data using a single-cell immune cell reference. This script imports the bulk and reference data to calculate cell type proportions using MIND, plot fractions by sample and case/control status and statistical tests to identify significant differences between cases and controls <code>deconv.md</code>. 
  

## Acknowledgements 

This work has emanated from research supported by a research grant from Science
Foundation Ireland (SFI) under Grant number 18/CRT/6214 and is co-funded under the European Regional Development Fund and by Genuity Science.


