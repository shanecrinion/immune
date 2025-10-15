# #!/usr/bin/env Rscript 

# Shane Crinion / shanecrinion@gmail.com / s.crinion1@universityofgalway.ie
# 19-06-2024

setwd('~/Documents/dev/dev2/')

# Read in the counts generated using FeatureCounts
counts <- read.table('counts.Homo_sapiens.GRCh37.75.txt', header=T)

# Extract the sample numbers from the column names
colnames(counts) <-  gsub("bam_files_fixed\\.(PSR\\d+)_aligned_reads_sorted\\.bam", 
                          "\\1", colnames(counts))

# Remove samples where QC failed for accompanying time point
counts <- counts[,!colnames(counts) %in% c('PSR013', 'PSR189', 'PSR007', 'PSR086')]

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

# Extract columns of interest and add to dds object
pheno_of_interest <- as.data.frame(pheno[,c('SubjectID','age','sex')])

# Extract T1 only 
#status <- subset(status, time.point=='T1')

# Merge the additional phenotypic information with the dds object
status <- merge(status,pheno_of_interest, 
                      by.x='sample', by.y='SubjectID', 
                      all.x=T, all.y=F)

table(status[status$time.point=='T1','Condition'])
table(status[status$time.point=='T3','Condition'])

range(status$age,na.rm = T)
mean(status$age,na.rm = T)

# Create a modified version of the status file including only the sequenced samples
status <- subset(status, status$ID %in% colnames(counts))

# Assign unique sample ID to row names
row.names(counts) <- counts$Geneid
# Filter to only include T1
counts.dds <- counts[,status$ID]

# Create a modified version of the counts file to remove additional columns for 'dds' creation
counts.dds <- counts[,-c(1:6)]
write.table(counts, 'remake/read_counts.txt', quote = F,sep = '\t')

# Format the data 
row.names(status) <- status$ID
status$time.point <- as.factor(status$time.point)
status$time.point <- relevel(status$time.point, ref='T1')
status$sample <- as.factor(status$sample)
status$ID <- as.factor(status$ID)
status$Condition <- as.factor(status$Condition)
status$Plate <- as.factor(status$Plate)
status$sex <- as.factor(status$sex)


write.table(status, 'rnaseq_metadata.txt', quote = F, sep='\t')


# Create DESeqDataSet object with two conditions, two timepoints and an interaction term
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts.dds, 
                              colData=status[names(counts.dds),], 
                              design=~ sex + condition)

# How many genes was dispersion estimated for?
nrow(dds) # 63677 - all genes from total RNA

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


# Filter results to only include those with > 10 reads in > 23 samples
#smallestGroupSize <- 23
#keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
#dds <- dds[keep,]
#nrow(dds) # 22933

# Perform size factor estimation and normalization
dds <- DESeq(dds) 

saveRDS(dds, 'remake/dds_unfiltered_T1.rds')
saveRDS(gene_info, 'remake/geneIDs.rds')
