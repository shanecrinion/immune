#!/usr/bin/env Rscript 

# Shane Crinion / shanecrinion@gmail.com / s.crinion1@universityofgalway.ie
# 04-07-2024

#### 1. Read and Format Data ####

# 1. Read in the count data
setwd('~/Documents/dev/dev2/')

# Read in the counts generated using FeatureCounts
counts <- read.table('counts.Homo_sapiens.GRCh37.75.txt', header=T)

# Extract the sample numbers from the column names
colnames(counts) <-  gsub("bam_files_fixed\\.(PSR\\d+)_aligned_reads_sorted\\.bam", 
                          "\\1", colnames(counts))

# Remove samples where QC failed for accompanying time point
#counts <- counts[,!colnames(counts) %in% c('PSR013', 'PSR189', 'PSR007', 'PSR086')]

# Read in the timepoint and condition info from the sample submission file
status <- read.csv('iRELATE mRNA and cDNA T3.xlsx - NGS sample submission format.csv')

# Extract the relevant columns
status <- status[c('sample', 'time.point', 'X', 'ID', 'X.1')]
colnames(status)[3] <- 'Condition'
colnames(status)[5] <- 'Plate'

# Import additional phenotypic information
#install.packages("readxl")
library("readxl")
pheno <- read_excel('remake/iRELATE_full_wf.xlsx')

# Extract list covar of interest of interest to extract from pheno file
covar_to_test <- c('tobacco', 'SubjectID', 'age', 'sex', 'bmi')

# Import olz info to combine with other covariates
olz <- read.csv('olz_patient.csv')

# Merge the pheno and olz data
status <- merge(status, olz, by.x='sample', by.y='SubjectID', all.x=T, all.y=F)

# # Merge the pheno and status data
status <- merge(status,pheno[,covar_to_test], 
                 by.x='sample', by.y='SubjectID', 
                 all.x=T, all.y=F)
# 
# # Import cellular deconvolution information
# cell_type_proportions <- read.csv('CIBERSORTx_results.csv')
# 
# # list the cell types of interst
# cell_types <- c('B.cells.naive', 'T.cells.CD8', 'T.cells.CD4.naive', 
#                 'T.cells.CD4.memory.resting','NK.cells.resting', 
#                 'Mast.cells.resting', 'Neutrophils')
# 
# # extract from the ref matrix
# cell_type_proportions <- cell_type_proportions[,c('Mixture', cell_types)]
# 
# # merge cell type fractions with other covariates
# status <- merge(status,cell_type_proportions, 
#                 by.x='ID', by.y='Mixture', 
#                 all.x=T, all.y=F)

# Extract T1 only 
status <- subset(status, time.point=='T1')

# subset samples
status <- subset(status, ID %in% colnames(counts))

# Assign unique sample ID to row names
row.names(counts) <- counts$Geneid

# Filter to only include T1
counts.dds <- counts[,status$ID]

# Format the data 
# Assign row names
row.names(status) <- status$ID

# Format time.point
status$time.point <- as.factor(status$time.point)
status$time.point <- relevel(status$time.point, ref='T1')

# Convert and format covariates using in design formula

# Case/Control Status
status$Condition <- as.factor(status$Condition)
status$Condition <- relevel(status$Condition, ref='Control')

# ID
status$ID <- as.factor(status$ID)

# Sample
status$sample <- as.factor(status$sample)

# Plate
status$Plate <- as.factor(status$Plate)

# Sex
status$sex <- as.factor(status$sex)

# Tobacco
status$tobacco <- as.factor(status$tobacco)

# Create a scaled and centred age to account for large SD
status$scaleAge <- as.numeric(scale(status$age, center = T))

# Assign value as 0 for Olanzepine levels in Controls
status[status$Condition=='Control' & is.na(status$OLZ.equivalents),'OLZ.equivalents'] <- 0

#### 2. Create DESeq2 Object ####

# Identify columns that have NAs
na_present <- sapply(status[row.names(status) %in% names(counts.dds),], function(x) any(is.na(x)))

# Print result
print(na_present)

# Identify and remove samples with NA detected in OLZ.equivalents and bmi
row.names(status[is.na(status$bmi),])
row.names(status[is.na(status$OLZ.equivalents),]) 
row.names(status[is.na(status$Neutrophils),])

# Add to list to remove
na_remove <- c(row.names(status[is.na(status$bmi),]), 
                row.names(status[is.na(status$OLZ.equivalents),]))
 
#na_remove <- row.names(status[is.na(status$bmi),])
#na_remove # 6 for Olz, 1 for BMI

# create the subset for this design
design_subset <- status[row.names(status) %in% names(counts.dds) &
         !row.names(status) %in% na_remove,]


library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts.dds[row.names(design_subset)], 
                              colData=design_subset,
                              design=~ sex + bmi + age +
                                tobacco + OLZ.equivalents + Condition)

# Perform size factor estimation and normalization
dds <- DESeq(dds) 
resultsNames(dds)

# extract results
results <- results(dds, alpha=0.05, test='Wald', name='Condition_Patient_vs_Control')

# count DEGs
sum(results$padj<0.05,na.rm=T) # w 467  
sum(results$padj<0.05 & abs(results$log2FoldChange) > 0.3,na.rm=T) # 304

# extract filtered results
filtered_results <- results[!is.na(results$padj) & results$padj < 0.05, ]
write.csv(x = as.data.frame(filtered_results), file='DEGs.padj05.csv',quote = F)

# annotate gene IDs
# Convert Ensembl names to Gene names 
# BiocManager::install('biomaRt')
library(biomaRt)

# Obtain correct Mart
ensembl <- useMart("ensembl")

# Extract human dataset
ensembl_dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

# Get Ensembl IDs
ensembl_ids <- row.names(dds)

# Get Gene IDs for list of ENSEMBL IDs
gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = ensembl_ids,
                   mart = ensembl_dataset)

# append gene names to results
counts_w_gene_names <- merge(counts(dds), gene_info, 
                 by.x='row.names', by.y='ensembl_gene_id', all.x=T, all.y=F)

# extract DEGs
results_w_gene_name.filtered <- merge(as.data.frame(filtered_results), gene_info, 
by.x='row.names', by.y='ensembl_gene_id', all.x=T, all.y=F)


write.table(x=counts_w_gene_names, file='counts.27-07-2024.txt', quote = F, row.names = T, sep='\t')
write.csv(x = results_w_gene_name, file='DEGs.padj05.w_geneNames.csv',quote = F)
saveRDS(counts(dds), file = 'dds.counts.rds')
saveRDS(dds, file='dds.rds')
