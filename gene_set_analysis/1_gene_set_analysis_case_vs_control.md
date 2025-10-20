Gene Set Analysis (Cases vs Controls)
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
results <- results(dds, alpha=0.05, test='Wald', name='Condition_Patient_vs_Control')
results$ensembl <- row.names(results)
```

We have to extract the Ensembl IDs and convert them to Entrez IDs, gene
symbol and descriptions. We do this using biomaRt.

Load the libraries and Ensembl dataset.

``` r
library(biomaRt)
library("AnnotationDbi")
library("org.Hs.eg.db")
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

    ##                                                  p.geomean stat.mean     p.val
    ## hsa00511 Other glycan degradation                0.9392496 -1.609945 0.9392496
    ## hsa04150 mTOR signaling pathway                  0.8853015 -1.210652 0.8853015
    ## hsa04540 Gap junction                            0.8914406 -1.240114 0.8914406
    ## hsa04120 Ubiquitin mediated proteolysis          0.8935649 -1.249065 0.8935649
    ## hsa00270 Cysteine and methionine metabolism      0.8530750 -1.059128 0.8530750
    ## hsa01040 Biosynthesis of unsaturated fatty acids 0.8334798 -0.984011 0.8334798
    ##                                                      q.val set.size      exp1
    ## hsa00511 Other glycan degradation                0.9392496       15 0.9392496
    ## hsa04150 mTOR signaling pathway                  0.8991497       48 0.8853015
    ## hsa04540 Gap junction                            0.8991497       79 0.8914406
    ## hsa04120 Ubiquitin mediated proteolysis          0.8991497      123 0.8935649
    ## hsa00270 Cysteine and methionine metabolism      0.8748094       32 0.8530750
    ## hsa01040 Biosynthesis of unsaturated fatty acids 0.8632552       16 0.8334798

``` r
head(KEGGres.downregulated[order(-KEGGres.downregulated$q.val),])
```

    ##                                                      p.geomean stat.mean
    ## hsa04020 Calcium signaling pathway                   0.9461460  1.613036
    ## hsa00230 Purine metabolism                           0.9281511  1.466743
    ## hsa00290 Valine, leucine and isoleucine biosynthesis 0.9299003  1.537567
    ## hsa04740 Olfactory transduction                      0.9336038  1.510517
    ## hsa03410 Base excision repair                        0.9035763  1.319145
    ## hsa03040 Spliceosome                                 0.8950511  1.257334
    ##                                                          p.val     q.val
    ## hsa04020 Calcium signaling pathway                   0.9461460 0.9461460
    ## hsa00230 Purine metabolism                           0.9281511 0.9394388
    ## hsa00290 Valine, leucine and isoleucine biosynthesis 0.9299003 0.9394388
    ## hsa04740 Olfactory transduction                      0.9336038 0.9394388
    ## hsa03410 Base excision repair                        0.9035763 0.9265974
    ## hsa03040 Spliceosome                                 0.8950511 0.9237387
    ##                                                      set.size      exp1
    ## hsa04020 Calcium signaling pathway                        169 0.9461460
    ## hsa00230 Purine metabolism                                145 0.9281511
    ## hsa00290 Valine, leucine and isoleucine biosynthesis       11 0.9299003
    ## hsa04740 Olfactory transduction                            85 0.9336038
    ## hsa03410 Base excision repair                              31 0.9035763
    ## hsa03040 Spliceosome                                      117 0.8950511

These analyses detected no significant pathways associated with
case-control differences.
