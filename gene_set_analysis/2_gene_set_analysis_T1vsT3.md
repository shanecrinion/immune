Gene Set Analysis (Timepoint 1 vs Timepoint 3)
================
Shane Crinion
2025-10-20

### 1. Import Data

Import DESeq2 object:

``` r
library(DESeq2)
# ├── import using dds 
dds <- readRDS('~/Documents/dev/dds_allParticipants_covar_included.rds')

# import gene names
gene_info <- readRDS('~/Documents/dev/dev2/remake/geneIDs.rds')
```

Extract the results table for the Patient vs Control analysis:

``` r
results <- results(dds, alpha=0.05, test='Wald', name='time.point_T3_vs_T1')
results$ensembl <- row.names(results)
```

We have to extract the Ensembl IDs and convert them to Entrez IDs, gene
symbol and descriptions. We do this using biomaRt.

Load the libraries and Ensembl dataset.

``` r
library(biomaRt)
library("AnnotationDbi")
library("org.Hs.eg.db")
```

    ## 

``` r
# Select the Ensembl dataset for human
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
```

We then extract the mapping information using the Ensembl mart.

``` r
# Retrieve the mapping between Ensembl IDs, Entrez IDs, gene symbols, and gene names
mapping <- getBM(
    attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "description"),
    filters = "ensembl_gene_id",
    values = row.names(results),
    mart = mart)
```

Merge this information to the results.

``` r
results <- merge(as.data.frame(results), mapping, by.x='ensembl', by.y='ensembl_gene_id', all.x=T)
```

### 2. KEGG and Gene Ontology pathway analysis (gene set test with GAGE and pathway visualisation with pathview)

We performed pathway analysis using the gage (Generally Applicable
Gene-set Enrichment for Pathway Analysis) package. This will identify
the significantly perturbed KEGG pathways. We used the pathview function
to visualise the gene expression pertubations in significant KEGG
pathways. GSEA uses the expression profiles of the genes.

Import the libraries and lists for KEGG pathways

``` r
# load libraries containing KEGG pathways
library(gage)
library(gageData)

# Import lists of genes in each KEGG pathway
data(kegg.sets.hs) 
# Import index of signalling and metabolic pathways
data(sigmet.idx.hs) 
# Import list of genes in Gene Ontology lists
data(go.subs.hs)

# Filter gene sets to only signalling and metabolic pathways
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
```

Extract fold changes for all genes.

``` r
foldchanges = results$log2FoldChange
names(foldchanges) = row.names(results)
```

Perform GAGE (Generally Applicable Gene-set Enrichment) to infer gene
sets that are significantly perturbed relative to all genes considered.

``` r
# Get the results
KEGGres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
GOres = gage(foldchanges, gsets=go.subs.hs, same.dir=TRUE)
```

Extract the upregulated and down regulated pathways

``` r
KEGGres.upregulated <- as.data.frame(KEGGres$greater)
KEGGres.downregulated <- as.data.frame(KEGGres$less)
```

View the results

``` r
head(KEGGres.upregulated[order(-KEGGres.upregulated$q.val),])
```

    ##                                           p.geomean stat.mean     p.val
    ## hsa04740 Olfactory transduction           0.2296305 0.7417934 0.2296305
    ## hsa00512 Mucin type O-Glycan biosynthesis 0.2375139 0.7225464 0.2375139
    ## hsa04614 Renin-angiotensin system         0.3127328 0.4934384 0.3127328
    ## hsa04744 Phototransduction                0.3154045 0.4836298 0.3154045
    ## hsa03060 Protein export                   0.3273122 0.4509276 0.3273122
    ## hsa00100 Steroid biosynthesis             0.3991872 0.2574644 0.3991872
    ##                                               q.val set.size      exp1
    ## hsa04740 Olfactory transduction           0.9996615       85 0.2296305
    ## hsa00512 Mucin type O-Glycan biosynthesis 0.9996615       22 0.2375139
    ## hsa04614 Renin-angiotensin system         0.9996615       17 0.3127328
    ## hsa04744 Phototransduction                0.9996615       26 0.3154045
    ## hsa03060 Protein export                   0.9996615       20 0.3273122
    ## hsa00100 Steroid biosynthesis             0.9996615       18 0.3991872

``` r
head(KEGGres.downregulated[order(-KEGGres.downregulated$q.val),])
```

    ##                                           p.geomean stat.mean     p.val
    ## hsa04740 Olfactory transduction           0.7703695 0.7417934 0.7703695
    ## hsa00512 Mucin type O-Glycan biosynthesis 0.7624861 0.7225464 0.7624861
    ## hsa04744 Phototransduction                0.6845955 0.4836298 0.6845955
    ## hsa04614 Renin-angiotensin system         0.6872672 0.4934384 0.6872672
    ## hsa03060 Protein export                   0.6726878 0.4509276 0.6726878
    ## hsa00450 Selenocompound metabolism        0.5984895 0.2527316 0.5984895
    ##                                               q.val set.size      exp1
    ## hsa04740 Olfactory transduction           0.7703695       85 0.7703695
    ## hsa00512 Mucin type O-Glycan biosynthesis 0.7672516       22 0.7624861
    ## hsa04744 Phototransduction                0.6959121       26 0.6845955
    ## hsa04614 Renin-angiotensin system         0.6959121       17 0.6872672
    ## hsa03060 Protein export                   0.6898264       20 0.6726878
    ## hsa00450 Selenocompound metabolism        0.6200696       14 0.5984895

These analyses detected no significant pathways associated with
time-point differences.
