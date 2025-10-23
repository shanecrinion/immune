# #!/usr/bin/env Rscript 

# Shane Crinion / shanecrinion@gmail.com / s.crinion1@universityofgalway.ie
# 19-06-2024

# 1. Read in the count data)

# Read in the counts generated using FeatureCounts
counts <- read.table('counts.Homo_sapiens.GRCh37.75.txt', header=T)

# Extract the sample numbers from the column names
colnames(counts) <-  gsub("bam_files_fixed\\.(PSR\\d+)_aligned_reads_sorted\\.bam", 
                          "\\1", colnames(counts))

# Extract ensembl id to relabel
ensembl_ids <- counts$Geneid

expression <- subset(counts, select=-c(Geneid, Chr, Start, End, Strand, Length))

row.names(expression) <- ensembl_ids

# Remove samples where QC failed for accompanying time point
#counts <- counts[,!colnames(counts) %in% c('PSR013', 'PSR189', 'PSR007', 'PSR086')]

#library("readxl")
# subjectID, age, sex, tobacco usage and bmi
pheno.1 <- read.csv('iRELATE_full_wf - iRELATE_full_wf.csv')
pheno.1 <- pheno.1[,c('SubjectID', 'tobacco', 'age', 'sex', 'bmi')]

# timepoint and condition case/control (from the sample submission file)
pheno.2 <- read.csv('iRELATE mRNA and cDNA T3.xlsx - NGS sample submission format.csv')
pheno.2 <- pheno.2[c('sample', 'time.point', 'X', 'ID', 'X.1')]
names(pheno.2) <- c('SubjectID', 'time.point', 'Condition', 'ID', 'Plate')
colnames(pheno.2)[3] <- 'Condition'
colnames(pheno.2)[5] <- 'Plate'

# olenzapine dosage 
pheno.3 <- read.csv('olz_patient.csv')

# merge the simpler datasets first
pheno = merge(pheno.1, pheno.3, by='SubjectID', all.x=T)
pheno = merge(pheno, pheno.2, by='SubjectID', all.x=T)

# Filter to only include samples that were sequenced and passed QC. Assign sample IDs as row.names
pheno <- subset(pheno, ID %in% rownames(t(expression)))
row.names(pheno) <- pheno$ID

# Assign value as 0 for Olanzepine levels in Controls
pheno[pheno$Condition=='Control' & is.na(pheno$OLZ.equivalents),'OLZ.equivalents'] <- 0

# Convert to factors and level as necessary
pheno = pheno %>% mutate(
  time.point = factor(time.point, levels=c('T1', 'T3')),
  Condition = factor(Condition, levels=c('Control', 'Patient')),
  tobacco = factor(tobacco, levels = c(0,1)),
  sex = as.factor(sex),
  OLZ.equivalents = as.numeric(scale(OLZ.equivalents, center=T)),
  age = as.numeric(scale(age, center = T)),
  bmi = as.numeric(scale(bmi, center=T))) %>%
  dplyr::select(!status) # no longer needed

# We will remove those samples as they will cause errors in the deseq object design
pheno = pheno[complete.cases(pheno), ]

# Now we create the DESeq object using the count matrix, meta-data and the design formula.
### Check that sample names match in both files
all(colnames(expression[,rownames(pheno)]) %in% rownames(pheno))
all(colnames(expression[,rownames(pheno)]) == rownames(pheno))

# Now create the full model DESeq DataSet object and generate the results.
library(DESeq2)
dds_all <- DESeqDataSetFromMatrix(countData = expression[,rownames(pheno)], colData=pheno, design=~ sex + bmi + age + tobacco + OLZ.equivalents + time.point + Condition + time.point:Condition)

# Perform size factor estimation and normalization
dds_all<- DESeq(dds_all) 

# Save dds object
saveRDS(dds_all, file='dds_covarIncluded.rds')

# Do the below to extract the differential expression results for the desired comparison (change 'name' to Cases vs Control if you wish to investigate those results)
# dds_all_results <- results(dds_all, alpha=0.05, test='Wald', name='time.point_T3_vs_T1')
