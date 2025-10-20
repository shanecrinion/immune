Data QC plots
================
2025-10-20

Here we will generate a number of plots to investigate the count data.
However, while we use raw counts and discrete distributions to explore
differential expression, transformed versions of the count data are
better for some downstream analyses such as clustering and heatmaps.

The report is based on the guidance
[here](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-transformations-and-visualization)
and uses descriptions from this tutorial for the transformation.

### 1. Library and Data Import

``` r
library(DESeq2)
library("vsn")
library("pheatmap")
library("RColorBrewer")
library(tidyverse)
```

Import DESeq2 object:

``` r
# ├── import using dds 
dds <- readRDS('dds_covarIncluded.rds')
dds
```

    ## class: DESeqDataSet 
    ## dim: 63677 180 
    ## metadata(1): version
    ## assays(4): counts mu H cooks
    ## rownames(63677): ENSG00000223972 ENSG00000227232 ... ENSG00000210195
    ##   ENSG00000210196
    ## rowData names(50): baseMean baseVar ... deviance maxCooks
    ## colnames(180): PSR120 PSR089 ... PSR072 PSR057
    ## colData names(11): SubjectID tobacco ... Plate sizeFactor

### 2. Extract Transformed Data

We will variance stabilizing transformations (VST) which incorporates a
prior on the sample differences. This transforms data on the log2 scale
which has been normalized with respect to library size or other
normalization factors.

The point of this transformation is to remove the dependence of the
variance on the mean, particularly the high variance of the logarithm of
count data when the mean is low. VST uses the experiment-wide trend of
variance over mean, in order to transform the data to remove the
experiment-wide trend. Note that we do not require or desire that all
the genes have exactly the same variance after transformation. Indeed,
in a figure below, you will see that after the transformation the genes
with the same mean do not have exactly the same standard deviations, but
that the experiment-wide trend has flattened. It is those genes with row
variance above the trend which will allow us to cluster samples into
interesting groups.

Here we generate the transformed counts using vst and extract a matrix
of normalised values using the assay function.

``` r
vsd <- vst(dds, blind=FALSE)
assay(vsd)[1:3,1:3]
```

    ##                   PSR120   PSR089   PSR114
    ## ENSG00000223972 6.850766 6.696634 6.205485
    ## ENSG00000227232 8.936643 8.703745 8.870793
    ## ENSG00000243485 3.667966 3.549275 3.167936

**Variance stabilizing transformation** Above, we used a parametric fit
for the dispersion. In this case, the closed-form expression for the
variance stabilizing transformation is used by the vst function. If a
local fit is used (option fitType=“locfit” to estimateDispersions) a
numerical integration is used instead. The transformed data should be
approximated variance stabilized and also includes correction for size
factors or normalization factors. The transformed data is on the log2
scale for large counts.

### 2.1 Effects of transformations on the variance

First we will generate normal transformed counts:

``` r
# BiocManager::install('vsn')
# this gives log2(n + 1)
ntd <- normTransform(dds)
```

Explore the transformations to determine how they change dependent on
the normalisation

``` r
table(counts(dds) == assay(dds), useNA = 'ifany') # conforming that these are the same
```

    ## 
    ##     TRUE 
    ## 11461860

``` r
assay(ntd)[1:6,1:6]; counts(dds, normalized=T)[1:6,1:6]
```

    ##                   PSR120    PSR089   PSR114    PSR149   PSR128   PSR085
    ## ENSG00000223972 6.631483 6.4516979 5.855893 5.7327512 5.806163 5.954638
    ## ENSG00000227232 8.886280 8.6444790 8.818059 8.8783571 8.602315 8.609485
    ## ENSG00000243485 1.063856 0.7062895 0.000000 1.7382567 1.313383 0.833736
    ## ENSG00000237613 0.000000 1.5334673 1.998407 0.8308821 0.000000 0.833736
    ## ENSG00000268020 0.000000 0.0000000 0.000000 0.0000000 0.000000 0.000000
    ## ENSG00000240361 0.000000 0.0000000 0.000000 0.0000000 0.000000 0.000000

    ##                     PSR120      PSR089     PSR114      PSR149     PSR128
    ## ENSG00000223972  98.145989  86.5295295  56.916106  52.1777631  54.953735
    ## ENSG00000227232 472.191260 399.1727201 450.336210 469.5998675 387.646619
    ## ENSG00000243485   1.090511   0.6316024   0.000000   2.3363177   1.485236
    ## ENSG00000237613   0.000000   1.8948072   2.995585   0.7787726   0.000000
    ## ENSG00000268020   0.000000   0.0000000   0.000000   0.0000000   0.000000
    ## ENSG00000240361   0.000000   0.0000000   0.000000   0.0000000   0.000000
    ##                      PSR085
    ## ENSG00000223972  61.0189953
    ## ENSG00000227232 389.5828163
    ## ENSG00000243485   0.7822948
    ## ENSG00000237613   0.7822948
    ## ENSG00000268020   0.0000000
    ## ENSG00000240361   0.0000000

Now, we’ll compare normal to variance transformed data:

``` r
par(mfrow=c(1,3), mar=c(4,4,.1,.1))
meanSdPlot(counts(dds))
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
```

<img src="Data_QC_plots_files/figure-gfm/sd-transform-1.png" width="50%" /><img src="Data_QC_plots_files/figure-gfm/sd-transform-2.png" width="50%" /><img src="Data_QC_plots_files/figure-gfm/sd-transform-3.png" width="50%" />

*Interpretation:*

Below we plot the standard deviation of the count data across samples,
against the mean. What we see is that transforming using the variance
stabilising method evens the standard deviation across samples. While
this is not appropriate when we are lddooking at the differential
expression, it can be useful for assess the quality of data across the
samples.

### 3. Data quality assessment by sample clustering and visualization

### 3.1 Heatmap of count matrix

Firstly, format the data as generate a heatmap.

``` r
#BiocManager::install('pheatmap')
#library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]
df <- as.data.frame(colData(dds)[,c("Condition", 'time.point', 'Plate', 'sex')])
```

Generate various orderings of the df to see the impact of ordering by
condition, timepoint and clustering

``` r
df_ordered_cond_tp <- as.data.frame(colData(dds)) %>% arrange(Condition, time.point)
df_ordered_cond_tp <- df_ordered_cond_tp[,c('Condition', 'time.point')]

df_ordered_tp_cond <- as.data.frame(colData(dds)) %>% arrange(time.point, Condition)
df_ordered_tp_cond <- df_ordered_tp_cond[,c('Condition', 'time.point')]

df_ordered_cond <- as.data.frame(colData(dds)) %>% arrange(Condition) 
df_ordered_cond <- df_ordered_cond[,c('Condition', 'time.point')]

df_ordered_tp <- as.data.frame(colData(dds)) %>% arrange(time.point) 
df_ordered_tp <- df_ordered_tp[,c('Condition', 'time.point')]
```

Comparison of condition and time.point ordering:

``` r
pheatmap(assay(vsd)[select,row.names(df_ordered_cond)], cluster_rows=F, show_rownames=F, show_colnames = F,
         cluster_cols=F, annotation_col=df_ordered_cond, main = 'Top 50 Read Counts Ordered by Condition')
pheatmap(assay(vsd)[select,row.names(df_ordered_tp)], cluster_rows=F, show_rownames=F, show_colnames = F,
         cluster_cols=F, annotation_col=df_ordered_tp, main = 'Top 50 Read Counts Ordered by Time Point')
```

<img src="Data_QC_plots_files/figure-gfm/heatmap-1.png" width="50%" /><img src="Data_QC_plots_files/figure-gfm/heatmap-2.png" width="50%" />

``` r
pheatmap(assay(vsd)[select,row.names(df_ordered_cond_tp)], cluster_rows=F, show_rownames=F, show_colnames = F,
         cluster_cols=F, annotation_col=df_ordered_cond,  main = 'Top 50 Read Counts Ordered by Condition, Time Point')
pheatmap(assay(vsd)[select,row.names(df_ordered_tp_cond)], cluster_rows=F, show_rownames=F, show_colnames = F,
         cluster_cols=F, annotation_col=df_ordered_tp, main = 'Top 50 Read Counts Ordered by Time Point, Condition')
```

<img src="Data_QC_plots_files/figure-gfm/heatmap-2-1.png" width="50%" /><img src="Data_QC_plots_files/figure-gfm/heatmap-2-2.png" width="50%" />

**Heatmap of sample-to-sample distances**

``` r
sampleDists <- dist(t(assay(vsd)))
```

``` r
#library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$sample, vsd$ID, vsd$Condition, vsd$time.point, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hmap_sampledist <- 
  pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, show_colnames = T, show_rownames = T, treeheight_col = 0)
```

![](Data_QC_plots_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
#ggsave('heatmap_sampledist.png', hmap_sampledist, height=25, width=11, dpi=1200)
```

Next, generate a PCA plot to identify groups within the samples

``` r
pc1 <- plotPCA(vsd, intgroup=c("sex"),ntop=50, pcsToUse = 1:2)
```

    ## using ntop=50 top features by variance

``` r
pc1
```

![](Data_QC_plots_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
#ggsave('pca_sex.png', dpi=1200)
```

Correct for sex effects:

``` r
# Create a sex corrected set
mat <- assay(vsd)
mm <- model.matrix(~ Condition, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$sex, design=mm)
assay(vsd) <- mat
plotPCA(vsd,  intgroup=c("sex"),ntop=200, pcsToUse = 1:2)
```

    ## using ntop=200 top features by variance

![](Data_QC_plots_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->
