# Author: Aldas Žarnauskas
# Title: Tumour Purity Estimation
# Date: 20/007/2023
# Last time updated: 28/008/2023
# Objective: Estimate tumour purity in the RCC samples
# Description: In this code I uploaded counts and coldata of RCC datasets. 
# I estimated tumour purity individually for each batch. Subsequently, I created
# new coldata file containing scaled tumour purity scores. 

#Clear environment
rm(list=ls(all.names = TRUE))


### ------------------------------- library --------------------------------- ##
################################################################################

library(DESeq2)
library(biomaRt)
library(dplyr)
library(org.Hs.eg.db)
library(CePa)
library(AnnotationDbi)
library(estimate)
library(tibble)

### ----------------------- Data Import + Wrangling ------------------------ ###
################################################################################

#Read in counts
counts <- read.csv("C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/RCC_counts_v2.csv",
                   header = T, sep = ',')
#check if counts have na values
any(is.na(counts))
counts[is.na(counts)] <- 0
any(is.na(counts))


#In your case this will be coldata_EOC.csv  
coldata <- read.csv(
  "C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/colData_kidney_cancer.csv")



#------------------ Tumour Purity Estimation by Batches -----------------------#
################################################################################

unique(coldata$batch) # "GSE188486" "GSE69197"  "GSE217386" "GSE143630" "GSE153262"
# all batches were generated with Illumina sequencing

GSE188486_colData <- coldata %>% 
  {.[.$batch == 'GSE188486',]} %>% 
  {rownames(.) <- NULL; .}
GSE188486_counts <- counts %>% 
  {.[, c(GSE188486_colData$run_accession)]} %>% as.matrix()


GSE69197_colData <- coldata %>% 
  {.[.$batch == 'GSE69197',]} %>% 
  {rownames(.) <- NULL; .}
GSE69197_counts <- counts %>% 
  {.[, c(GSE69197_colData$run_accession)]} %>% as.matrix()


GSE217386_colData <- coldata %>% 
  {.[.$batch == 'GSE217386',]} %>% 
  {rownames(.) <- NULL; .}
GSE217386_counts <- counts %>% 
  {.[, c(GSE217386_colData$run_accession)]} %>% as.matrix()


GSE143630_colData <- coldata %>% 
  {.[.$batch == 'GSE143630',]} %>% 
  {rownames(.) <- NULL; .}
GSE143630_counts <- counts %>% 
  {.[, c(GSE143630_colData$run_accession)]} %>% as.matrix()


GSE153262_colData <- coldata %>% 
  {.[.$batch == 'GSE153262',]} %>% 
  {rownames(.) <- NULL; .}
GSE153262_counts <- counts %>% 
  {.[, c(GSE153262_colData$run_accession)]} %>% as.matrix()


GSE150474_colData <- coldata %>% 
  {.[.$batch == 'GSE150474',]} %>% 
  {rownames(.) <- NULL; .}
GSE150474_counts <- counts %>% 
  {.[, c(GSE150474_colData$run_accession)]} %>% as.matrix()


# Set the dataset_name, all_coldata, all_counts as the name of the dataset

# Choose the title of the dataset used to name the file of the estimateScore() function 
dataset_name <- "GSE150474_scores_ctl.gct"

#Remove cell line column and use as rownames instead
all_coldata <- GSE150474_colData %>% 
  column_to_rownames(., var = 'run_accession')

#Editing counts data
#Create a copy of counts data
all_counts <- GSE150474_counts 

# quality -check, ensure columns and rows match in the same order
all(rownames(all_coldata) %in% colnames(all_counts))
all(rownames(all_coldata) == colnames(all_counts))



#-------------------- Set Counts Rownames to EntrezIDs ------------------------#
################################################################################
# Use not normalised gene expression raw counts as the input because estimate
# package estimates the proportion of the stromal, immune, and tumour purity 
# using single sample Gene Set Enrichment Analsysi.

# uncomment if you haven't yet loaded human ensembl dataset
# ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

annotations <- getBM(
  attributes = c("ensembl_gene_id_version", 'entrezgene_id'),
  filters = "ensembl_gene_id_version",
  values = rownames(all_counts),
  mart = ensembl)

sum(is.na(annotations$entrezgene_id))
# 35973 ensemblID versions did map to entrezID
sum(!is.na(annotations$entrezgene_id))
# 30286 ensemblID versions did not map to entrezID

#------------------------ Estimating Tumor Purity -----------------------------#
################################################################################

# common_genes Contains info on 10412 common genes from estimateR package
data("common_genes")

# filter annotation table to contain only common genes
annotations_common <- annotations[annotations$entrezgene_id %in% common_genes$EntrezID,] %>% 
  {.[!duplicated(.$entrezgene_id),]; .} # filter out duplicated genes

nrow(annotations_common)
#10160 common genes

# swap ensemblIDs as rownames in the counts table to the entrezIDs
counts_filtered <- as.data.frame(all_counts[
  rownames(all_counts) %in% annotations_common$ensembl_gene_id_version,]) %>% 
  {.$ensembl_gene_id_version <- rownames(.); .} %>% # add a column for ensembl_gene_id_version
  left_join(., annotations_common, by = 'ensembl_gene_id_version') %>% 
  .[!duplicated(.$entrezgene_id),] %>% 
  {rownames(.) <- as.character(.$entrezgene_id); .} %>% 
  dplyr::select(-c('ensembl_gene_id_version', 'entrezgene_id'))


#Write as text file
write.table(counts_filtered, "counts_filtered_withctrl.txt", quote = F, row.names = T, sep = "\t")

#Run filterCommonGenes
#Unifies different number of genes per platform against 10,412 common genes.

#Specify counts file output
out.file <- "counts_ctl.gct"
#Run function
filterCommonGenes(input.f = "counts_filtered_withctrl.txt", output.f = out.file, id = "EntrezID")

#Estimate score
#Specify output file
out.file <- dataset_name
#Calculate scores
estimateScore(input.ds = "counts_ctl.gct",
                      output.ds = out.file,
                      platform = "illumina") #Specifies illumina RNAseq platform



# add tumour purity scores to the colData_kidney_cancer.csv
GSE188486_scores <- read.gct(
  "c:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/GSE188486_scores_ctl.gct")
GSE69197_scores <- read.gct(
  "c:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/GSE69197_scores_ctl.gct")
GSE217386_scores <- read.gct(
  "c:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/GSE217386_scores_ctl.gct")
GSE143630_scores <- read.gct(
  "c:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/GSE143630_scores_ctl.gct")
GSE153262_scores <- read.gct(
  "c:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/GSE153262_scores_ctl.gct")
GSE150474_scores <- read.gct(
  "c:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/GSE150474_scores_ctl.gct")


kidney_cancer_metadata.df <- coldata %>% 
  {.$purity <- c(GSE188486_scores[3,], GSE69197_scores[3,], 
                 GSE217386_scores[3,], GSE143630_scores[3,], 
                 GSE153262_scores[3,], GSE150474_scores[3,]) %>% unlist(); .} %>% 
  {.$purityscaled <- scale(.$purity, center = T); .} 



write.csv(kidney_cancer_metadata.df, 
          'C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/colData_kidney_cancer_tumourPurity.csv',
          row.names = F)
