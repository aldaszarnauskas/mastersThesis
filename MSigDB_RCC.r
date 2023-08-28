# Author: Aldas Žarnauskas
# Title: Obtain RCC cancer Metadata
# Date: 20/007/2023
# Last time updated: 20/007/2023
# Objective: obtain RCC metadata

#Clear environment
rm(list=ls(all.names = TRUE))

### ----------------------------- Libraries --------------------------------- ##
################################################################################
library(GEOquery)
library(dplyr)
library(HelpersMG)


### -------------------- Obtaining Metadata --------------------------------- ##
################################################################################

# Define the GEO accession number
gse_list <- c("GSE188486", "GSE217386", "GSE69197", "GSE143630", "GSE153262", 
              'GSE150474')

metadata <- list()

for (gse in gse_list){
  gse_data <- getGEO(gse, GSEMatrix = TRUE, getGPL = FALSE)
  gse_data <- gse_data[[1]] %>% as.data.frame()
  metadata[[gse]] <- gse_data
}

names(metadata) 
# "GSE188486" "GSE217386" "GSE152938" "GSE69197"  "GSE143630" "GSE153262"


### ------------------------- Data Import+Cleaning -------------------------- ##
################################################################################

# PRJNA779050 ------------------------------------------------------------------
gse_data_GSE188486 <- metadata[['GSE188486']] 
names(gse_data_GSE188486) 
# 2-geo_accession, 44-subtype
# Filter columns and rows of interest
gse_data_GSE188486_filtered <- gse_data_GSE188486[
  gse_data_GSE188486$library_strategy == "RNA-Seq",] %>% 
  .[, c(2, 44)] 
# Clean the dataframe
gse_data_GSE188486_cleaned <- gse_data_GSE188486_filtered %>% 
  {colnames(.) <- c('sample', 'subtype'); .} %>% 
  {.$subtype <- ifelse(.$subtype == 'Chromophobe RCC', 'chRCC', 
                       ifelse(
                         .$subtype == 'Clear cell RCC', 'ccRCC', 'pRCC')); .} %>% 
  {.$batch <- 'GSE188486'; .} %>% 
  {.$paired_single <- 'paired-end'; .}


# PRJNA898782 ------------------------------------------------------------------
gse_data_GSE217386 <- metadata[['GSE217386']]
names(gse_data_GSE217386) 
# 2-geo_accession, 8-subtype, 13-grade, 14-stage, 17-metastatic
# Filter columns and rows of interest
gse_data_GSE217386_filtered <- gse_data_GSE217386[
  gse_data_GSE217386$library_strategy == "RNA-Seq",] %>% # filter out non-RNA-seq data
  .[, c(2, 8, 13, 14, 17, 24)] # Select columns of interest
gse_data_GSE217386_cleaned <- gse_data_GSE217386_filtered %>% 
  {colnames(.) <- c('sample', 'subtype', 'grade', 'stage', 'metastatic', 'patient'); .} %>% 
  {.$subtype <- ifelse(.$subtype == 'clear cell renal cell carcinoma', 'ccRCC', 'normal'); .} %>% 
  {.$grade <- .$grade %>% lapply(function(x) sub("tumor grade: ", "", x)); .} %>%
  {.$grade <- .$grade %>% as.character(); .} %>% 
  {.$stage <- ifelse(.$stage == 'tumor stage: Stage I', 'Stage1', 
                     ifelse(.$stage == 'tumor stage: Stage II', 'Stage2', 'Stage3')); .} %>% 
  {.$metastatic <- .$metastatic %>% lapply(function(x) sub("m: ", "", x)); .} %>%
  {.$metastatic <- .$metastatic %>% as.character(); .} %>%
  {.$batch <- 'GSE217386'; .} %>% 
  {.$paired_single <- 'paired-end'; .} %>% 
  mutate(patient = str_sub(patient, 1, 2)) %>% 
  {.$patient <- setNames(paste('patient', sep = '', replicate(length(unique(.$patient)), 
                                   paste(sample(0:9, 4), collapse = ""), 
                                   simplify = FALSE)), unique(.$patient));.}
  

# PRJNA284822 ------------------------------------------------------------------
gse_data_GSE69197 <- metadata[['GSE69197']]
names(gse_data_GSE69197) 
# 2-geo_accession, 8-subtype
# Filter columns and rows of interest
gse_data_GSE69197_filtered <- gse_data_GSE69197[gse_data_GSE69197$library_strategy == "RNA-Seq",] %>% 
  .[, c(2, 8, 17)] 

gse_data_GSE69197_cleaned <- gse_data_GSE69197_filtered %>% 
  {colnames(.) <- c('sample', 'subtype', 'patient'); .} %>% 
  {.$subtype <- ifelse(.$subtype == 'RNA sequencing Uninvolved Kidney', 'normal', 'ccRCC'); .}%>% 
  {.$batch <- 'GSE69197'; .} %>% 
  {.$paired_single <- 'paired-end'; .} %>% 
  mutate(patient = str_sub(patient, -1, -1)) %>% 
  {.$patient <- setNames(
    paste('patient', sep = '', replicate(length(unique(.$patient)),
                                         paste(sample(0:9, 4), collapse = ""), 
                                         simplify = FALSE)), unique(.$patient));.}


# PRJNA601152 ------------------------------------------------------------------
gse_data_GSE143630 <- metadata[['GSE143630']]
names(gse_data_GSE143630) 
# 2-geo_accession, 8-subtype,13-stage, 12-metastatic
# Filter columns and rows of interest
gse_data_GSE143630_filtered <- gse_data_GSE143630[gse_data_GSE143630$library_strategy == "RNA-Seq",] %>% 
  .[, c(2, 8, 12, 13)] 

gse_data_GSE143630_cleaned <- gse_data_GSE143630_filtered %>% 
  {colnames(.) <- c('sample', 'subtype', 'metastatic', 'stage'); .} %>% 
  {.$subtype <- ifelse(.$subtype == 'Renal Cell Carcinoma', 'ccRCC', 'normal'); .} %>% 
  {.$stage <- ifelse(.$stage == 'tumor stage: 1', 'Stage1', 
                     ifelse(.$stage == 'tumor stage: 2', 'Stage2', 'Stage3')); .} %>% 
  {.$metastatic <- ifelse(.$metastatic == 'metastatic: No', 'M0', 
                     ifelse(.$metastatic == 'metastatic: Yes', 'M1', 3)); .} %>% 
  {.$batch <- 'GSE143630'; .} %>% 
  {.$paired_single <- 'single-end'; .}


# PRJNA641885 ------------------------------------------------------------------
gse_data_GSE153262 <- metadata[['GSE153262']]
names(gse_data_GSE153262) 
# 2-geo_accession, 8-subtype,13-stage, 14-treatment 
# Filter columns and rows of interest
gse_data_GSE153262_filtered <- gse_data_GSE153262[
  gse_data_GSE153262$library_strategy == "RNA-Seq",] %>% 
  .[, c(2, 8, 13, 14)]

gse_data_GSE153262_cleaned <- gse_data_GSE153262_filtered %>%
  {colnames(.) <- c('sample', 'subtype', 'stage', 'treatment'); .} %>%
  {.$subtype <- ifelse(
    .$subtype == 'renal cell carcinoma', 'RCC', 'normal'); .} %>% 
  {.$stage <- ifelse(.$stage == 'pathologic t stage: 1B', 'Stage1', 
                     ifelse(.$stage == 'pathologic t stage: 2A', 'Stage2', 
                            ifelse(.$stage == 'pathologic t stage: 2', 'Stage2', 
                                   ifelse(.$stage == 'pathologic t stage: 3A', 'Stage3', 'Stage4')))); .} %>%
  {.$treatment <- ifelse(.$treatment == 'treatment: nephrectomy only', 'Nephrectomy', 
                         'SBRT_nephrectomy'); .} %>% 
  {.$batch <- 'GSE153262'; .} %>% 
  {.$paired_single <- 'single-end'; .}



# PRJNA632545 ------------------------------------------------------------------
gse_data_GSE150474 <- metadata[['GSE150474']]
names(gse_data_GSE150474) 
# 2-geo_accession, 8-subtype
# Filter columns and rows of interest
gse_data_GSE150474_filtered <- gse_data_GSE150474[
  gse_data_GSE150474$library_strategy == "RNA-Seq"
  & gse_data_GSE150474$source_name_ch1 == 'normal tissue',] %>% 
  .[, c(2, 8)] 

gse_data_GSE150474_cleaned <- gse_data_GSE150474_filtered %>% 
  {colnames(.) <- c('sample', 'subtype'); .} %>% 
  {.$subtype <- ifelse(.$subtype == 'normal tissue', 'normal', 'normal'); .} %>% 
  {.$batch <- 'GSE150474'; .} %>% 
  {.$paired_single <- 'paired-end'; .}
  


# Merging metadata & adding run_accession --------------------------------------

kidney_cancer_metadata.df <- dplyr::bind_rows(gse_data_GSE153262_cleaned, 
                                              gse_data_GSE143630_cleaned, 
                                              gse_data_GSE69197_cleaned, 
                                              gse_data_GSE217386_cleaned, 
                                              gse_data_GSE188486_cleaned,
                                              gse_data_GSE150474_cleaned) %>% 
  {.[is.na(.)] <- NA; .}


# some samples from normal tissue are assigned a stage in the metadata. Those
# samples are collected from adjacent cancer tissues. The stage of the adjacent
# cancer tissue is used to label those normal samples. Although, to avoid 
# confusion, the stage will be removed if samples comes from a normal tissue.

kidney_cancer_metadata.df <- kidney_cancer_metadata.df %>%
  mutate(stage = ifelse(subtype == 'normal', 'normal', stage),
         metastatic = ifelse(subtype == 'normal', 'normal', metastatic),
         treatment = ifelse(subtype == 'normal', 'normal', treatment),
         early_late_stage = ifelse(stage == 'Stage1', 'early', ifelse(
           stage == 'Stage2', 'early', ifelse(stage == 'Stage3', 'late', ifelse(
             stage == 'Stage4', 'late', ifelse(
               subtype == 'normal', 'normal', NA
             ))
           )
         )))

# download files containing run_accession (e.g. SRR16892490) & sample alias 
# (e.g. GSM5683704) so that you could merge them with the final colData
# PRJNA779050 PRJNA284822 PRJNA898782 PRJNA601152 PRJNA641885
url <- c('https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA779050&result=read_run&fields=run_accession,sample_alias&format=tsv&download=true&limit=0',
         'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA284822&result=read_run&fields=run_accession,sample_alias&format=tsv&download=true&limit=0',
         'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA898782&result=read_run&fields=run_accession,sample_alias&format=tsv&download=true&limit=0',
         'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA601152&result=read_run&fields=run_accession,sample_alias&format=tsv&download=true&limit=0',
         'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA641885&result=read_run&fields=run_accession,sample_alias&format=tsv&download=true&limit=0',
         'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA632545&result=read_run&fields=run_accession,sample_alias&format=tsv&download=true&limit=0')

run_alias <- c(
  'PRJNA779050_run_alias.tsv',
  'PRJNA284822_run_alias.tsv',
  'PRJNA898782_run_alias.tsv',
  'PRJNA601152_run_alias.tsv',
  'PRJNA641885_run_alias.tsv',
  'PRJNA632545_run_alias.tsv'
  )

download.file(url, run_alias)

PRJNA779050_alias <- read.csv('C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/PRJNA779050_run_alias.tsv',
                              sep = '\t', header = T)
PRJNA284822_alias <- read.csv('C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/PRJNA284822_run_alias.tsv',
                              sep = '\t', header = T)
PRJNA898782_alias <- read.csv('C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/PRJNA898782_run_alias.tsv',
                              sep = '\t', header = T)
PRJNA601152_alias <- read.csv('C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/PRJNA601152_run_alias.tsv',
                              sep = '\t', header = T)
PRJNA641885_alias <- read.csv('C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/PRJNA641885_run_alias.tsv',
                              sep = '\t', header = T)
PRJNA632545_alias <- read.csv('C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/PRJNA632545_run_alias.tsv',
                                  sep = '\t', header = T)

kidney_cancer_metadata.df <- rbind(
  PRJNA779050_alias[,c(1,3)],
  PRJNA284822_alias[,c(1,3)],
  PRJNA898782_alias[,c(1,3)],
  PRJNA601152_alias[,c(1,3)],
  PRJNA641885_alias[,c(1,3)],
  PRJNA632545_alias[,c(1,3)]) %>% 
  as.data.frame() %>% 
  {colnames(.) <- c('run_accession', 'sample'); .} %>% 
  inner_join(., kidney_cancer_metadata.df, by = 'sample') %>% 
  mutate()


write.csv(kidney_cancer_metadata.df, 
          'C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/colData_kidney_cancer.csv',
          row.names = F)

