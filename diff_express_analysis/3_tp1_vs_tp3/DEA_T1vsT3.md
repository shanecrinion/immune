DEA: Timepoint 1 vs Timepoint 3
================
Shane Crinion
2025-10-20

Libraries:

``` r
library(DESeq2)
library(EnhancedVolcano)
library(apeglm)
library(ashr)
```

Import DESeq2 object:

``` r
# ├── import using dds 
dds <- readRDS('../dds_covarIncluded.rds')
```

Import gene names for labelling:

``` r
# import gene names
gene_info <- readRDS('~/Documents/dev/dev2/remake/geneIDs.rds')
```

Extract the results table for the T1 vs T3 analysis:

``` r
results <- results(dds, alpha=0.05, test='Wald', name='time.point_T3_vs_T1')
results$ensembl <- row.names(results)
```

Annotate the results with gene names:

``` r
results_w_gene_name <- merge(as.data.frame(results),
                             gene_info, 
                             by.x='ensembl',
                             by.y='ensembl_gene_id',
                             all.x=T, all.y=F)
#write.csv(subset(results_w_gene_name, padj<0.05), 'time.point_T3_vs_T1_padj05.csv')
```

### 2. Visualisations

We will generate a number of visualisations to understand the underlying
patterns of the data.

We will make MA plots, volcano plots, heatmaps, plot counts, and
clustering dendrograms

### 2.1. MA-plot

MA-plots are used to show the log2 fold change (M) between T1 and T3
over the average expression level (A) of genes across all samples. Blue
indicates the adjusted p-value is less than 0.05.

``` r
plotMA(results, alpha=0.05, ylim=c(-1.5,1.5))
```

![](DEA_T1vsT3_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

The above is based on a normal prior distribution, centered on zero and
with a scale that is fit to the data. Shrunken log2 fold changes can
remove the noise associated with low read counts and create improved
visualisations. Here we shrink log2 fold change values using the ASHR
and APEGLM methods to improve visualisations.

First, we generate the shrunk results using the adaptive t prior
shrinkage estimator (apeglm) and the adaptive shrinkage estimator
(ashr):

``` r
resultsASHR <- lfcShrink(dds, res = results,
                         coef='time.point_T3_vs_T1',
                         type="ashr")
```

    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041

``` r
resultsApeglm <- lfcShrink(dds, res = results,
                           coef='time.point_T3_vs_T1',
                           type="apeglm")
```

    ## using 'apeglm' for LFC shrinkage. If used in published research, please cite:
    ##     Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
    ##     sequence count data: removing the noise and preserving large differences.
    ##     Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

Next, we compare the normal and shrunken estimators.

``` r
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-1.5,1.5)
plotMA(results, xlim=xlim, ylim=ylim, main="normal")
plotMA(resultsASHR, xlim=xlim, ylim=ylim, main="ashr")
plotMA(resultsApeglm, xlim=xlim, ylim=ylim, main="apeglm")
```

![](DEA_T1vsT3_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

*Interpretation:* The effect change estimators do generate improved
visualisation of the gene expression changes linked to schizophrenia.
The apeglm method appears to represent these changes most clearly. Blue
points indicate genes that are significantly differentially expressed.
There are more significantly upregulated genes than significantly
downregulated genes between cases and controls. Many genes are far from
the centre line, indicating that expression of these genes is likely
linked to the schizophrenia. The distribution of points is generally
evenly spread, indicating that there does not appear to be any technical
artifacts or batch effects.

### 2.2. Volcano plots

We generate volcano plots to compare log fold change to statistical
significance. Here we can see individual genes and may identify genes
that have been previously linked to schizophrenia or the immune
response.

First we investigated IL genes due to the immune aspect of the project.

``` r
#BiocManager::install('EnhancedVolcano')
  EnhancedVolcano(results_w_gene_name,
    lab = results_w_gene_name$external_gene_name,
    selectLab = results_w_gene_name[grep('^IL', results_w_gene_name$external_gene_name),
                                    'external_gene_name'],
    FCcutoff = 0.3,
    x = 'log2FoldChange', subtitle= NULL, 
    caption = NULL, title = NULL,
    xlim=c(-1.362953,1.745288), # min, max 2fc
    y = 'pvalue',
    ylim=c(0,8.5)) # min, max log10 pvalue + .5
```

![](DEA_T1vsT3_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
#ggsave('vp_IL.png', dpi=1200)
#getwd()
```

*Interpretation* The upper right quadrant shows that there are a number
of genes that are signifcantly upregulated. We added labels to all genes
beginning in “IL” to investigate immune related genes however, we did
not find that any of these had significant differential expression. We
will investigate what genes do show as significant in the next plot.

``` r
#BiocManager::install('EnhancedVolcano')
  EnhancedVolcano(results_w_gene_name,
    lab = results_w_gene_name$external_gene_name,
  #  selectLab = results_w_gene_name[grep('^IL', results_w_gene_name$external_gene_name),
   #                                 'external_gene_name'],
    FCcutoff = 0.5,
    labSize = 4,
    x = 'log2FoldChange', subtitle= NULL, 
    caption = NULL, title = NULL,
    xlim=c(-1.362953,1.745288), # min, max 2fc
    y = 'pvalue',
  drawConnectors = T,
    ylim=c(0,8.5)) # min, max log10 pvalue + .5
```

![](DEA_T1vsT3_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
#ggsave('volcanoPlot_SZstatusVsControls.png', dpi=1200)
```

A number of significantly up and down regulated genes were detected
between cases and controls. We take a further look at these below to
investigate the number of genes above padj/log2FoldChange thresholds:

``` r
sum(results_w_gene_name$padj < 0.05, na.rm = T) # 770
```

    ## [1] 770

``` r
sum(results_w_gene_name$padj<0.05 & abs(results_w_gene_name$log2FoldChange) > 0.5,na.rm=T) # 44
```

    ## [1] 44

``` r
sum(results_w_gene_name$padj<0.05 & abs(results_w_gene_name$log2FoldChange) > 0.3,na.rm=T) # 273
```

    ## [1] 273

Let’s take a look at the top genes:

``` r
head(results_w_gene_name[order(results_w_gene_name$padj, results_w_gene_name$pvalue),])
```

    ##               ensembl  baseMean log2FoldChange      lfcSE      stat
    ## 12507 ENSG00000168439 283.25839     -0.2400624 0.04007602 -5.990176
    ## 29348 ENSG00000226891 140.22433      0.5091531 0.09401439  5.415694
    ## 17800 ENSG00000197451 307.97036     -0.1661974 0.03072446 -5.409286
    ## 12615 ENSG00000168906 614.32340     -0.1866110 0.03473660 -5.372172
    ## 26426 ENSG00000223336  23.32137      0.5750913 0.10834607  5.307911
    ## 9578  ENSG00000152022 430.96991     -0.1985895 0.03756207 -5.286970
    ##             pvalue         padj external_gene_name
    ## 12507 2.096135e-09 4.471476e-05              STIP1
    ## 29348 6.105160e-08 4.148746e-04          LINC01359
    ## 17800 6.327636e-08 4.148746e-04            HNRNPAB
    ## 12615 7.779385e-08 4.148746e-04              MAT2A
    ## 26426 1.108886e-07 4.421377e-04            RNU2-6P
    ## 9578  1.243590e-07 4.421377e-04               <NA>
