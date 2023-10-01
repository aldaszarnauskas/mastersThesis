# Author: Aldas Žarnauskas
# Title: Differential Expression Analysis
# Date: 21/007/2023
# Last time updated: 21/007/2023
# Objective: perform DEA on RCC datasets

#Clear environment
rm(list=ls()) 

### ------------------------------- library --------------------------------- ##
################################################################################
library(dplyr)
library(DESeq2)
library(tibble)
library(tidyr)


### ----------------------- Data Import + Wrangling ------------------------- ##
################################################################################


# load RAAS genes 
gene_set.df <- read.csv('c:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/ras_raas_genes.csv')

# read in counts
counts <- read.csv("C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/RCC_counts.csv",
                   header = T, sep = ',')

# check if counts have na values
any(is.na(counts))
counts[is.na(counts)] <- 0
any(is.na(counts))

 
coldata <- read.csv(
  "C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/colData_kidney_cancer.csv")


### ---------------- Differential Expression Analysis ---------------------- ###
################################################################################


differential_expression <- function(counts, coldata, comparison, batch = F,
                                    purityscaled = F, add_title = '',
                                    matched = F){
  
  
  # Constructs a formula based on set parameters
  if ('normal' %in% unique(coldata[[comparison]])){
    
    ref_val <- 'normal'
    comparisons <- unique(coldata[[comparison]]) %>% 
      .[. != ref_val] %>%
      c(ref_val, .)
    
  } else {
    
    comparisons <- unique(coldata[[comparison]]) %>% sort()
    
  }
  
  
  if (all(batch, purityscaled, matched)){
    
    deseq_formula <- as.formula(paste("~ batch + purityscaled + patient + ", 
                                      comparison))
    
  } else if (all(batch, purityscaled)){
    
    deseq_formula <- as.formula(paste("~ batch + purityscaled + ", comparison))
    
  } else if (all(batch, matched)){
    
    deseq_formula <- as.formula(paste("~ batch + patient + ", comparison))
    
  } else if (batch){
    
    deseq_formula <- as.formula(paste("~ batch + ", comparison))
    
  } else if (all(purityscaled, matched)){
    
    deseq_formula <- as.formula(paste("~ purityscaled + patient + ", comparison))
    
  } else if (purityscaled){
    
    deseq_formula <- as.formula(paste("~ purityscaled + ", comparison))
    
  } else if (matched){
    
    deseq_formula <- as.formula(paste("~ patient + ", comparison))
    
  } else {
    
    deseq_formula <- as.formula(paste("~ ", comparison))
    
  }

  
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = coldata,
                                design = deseq_formula)
  
  
  # specify control level
  dds[[comparison]] <- relevel(dds[[comparison]], ref = comparisons[1])
  
  # run DESeq
  dds <- DESeq(dds)
  
  # count the number of subtypes
  nr_comparisons <- length(comparisons)
  
  # eliminate rows with low counts
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep,]
  
  
  # store the results of pairwise comaparisons between all subtypes
  result_tabs <- list()
  
  
  for (i in 1:(nr_comparisons-1)){
    
    a_comparison <- comparisons[i] # reference subtype
    
    
    # relevel a reference subtype after the first interation
    if (i > 1){
      dds[[comparison]] <- relevel(dds[[comparison]], ref = a_comparison)
      dds <- nbinomWaldTest(dds)  
    }
    
    
    # this variable will be used to obtain coefficients for the lfcshrink()
    result_names <- resultsNames(dds)
    
    
    # obtain coefficients of subtypes
    coefs <- list()
    for (comp in 2:nr_comparisons){
      comp <- comparisons[comp]
      coefs[[comp]] <- grep(comp, result_names)
    }
    if (nr_comparisons < 3){
      
      coefficient <- coefs[[comparisons[2]]]
      res_shrunk <- lfcShrink(dds, coef = coefficient)
      category <- paste(add_title, comparisons[2], '_vs_', a_comparison, sep = '')
      result_tabs[[category]] <- res_shrunk
      
    } else{
      
      # performs multiple pairwise comparisons + lfcshrinkage
      for (j in (i+1):nr_comparisons){
        coefficient <- coefs[[comparisons[j]]]
        res <- results(dds, contrast = c(comparison, comparisons[j], a_comparison))
        res_shrunk <- lfcShrink(dds, coef = coefficient, res = res)
        category <- paste(add_title, comparisons[j], '_vs_', a_comparison, sep = '')
        result_tabs[[category]] <- res_shrunk 
        
      }
    }
  }
  
  return(result_tabs)
  
}

### ------------------------------ Peforming DE ---------------------------- ###
################################################################################

### Normal vs Cancer
normalVScancer_colData <- coldata %>% 
  filter(batch == 'GSE188486' | batch == 'GSE217386' | batch == 'GSE69197') %>% 
  {rownames(.) <- NULL; .} %>% 
  mutate(subtype = ifelse(subtype == 'normal', 'normal', 'cancer'))
normalVScancer_counts <- counts %>% 
  {.[, c(normalVScancer_colData$run_accession)]} %>% as.matrix()
normalVScancer_colData <- normalVScancer_colData %>% 
  column_to_rownames(., var = 'run_accession') %>% 
  {.$run_accession <- rownames(.);.}
# Quality check, ensure columns and rows match in the same order
all(rownames(normalVScancer_colData) %in% colnames(normalVScancer_counts))
all(rownames(normalVScancer_colData) == colnames(normalVScancer_counts))

DE_normalVStumour <- differential_expression(
  coldata =  normalVScancer_colData,
  counts = normalVScancer_counts, comparison = 'subtype', purityscaled = F,
  batch = T)
DE_normalVStumourPE <- differential_expression(
  coldata =  normalVScancer_colData,
  counts = normalVScancer_counts, comparison = 'subtype', purityscaled = T,
  batch = T,
  add_title = 'PE')


### Subtypes
subtypes_colData <- coldata %>% 
  filter((batch == 'GSE188486' | batch == 'GSE217386' | .$batch == 'GSE69197') &
           !(patient %in% c('patient6801', 'patient3596'))) %>% 
  {rownames(.) <- NULL; .}
subtypes_counts <- counts %>% 
  {.[, c(subtypes_colData$run_accession)]} %>% as.matrix()
subtypes_colData <- subtypes_colData %>% 
  column_to_rownames(., var = 'run_accession') %>% 
  {.$run_accession <- rownames(.);.}
# Quality check, ensure columns and rows match in the same order
all(rownames(subtypes_colData) %in% colnames(subtypes_counts))
all(rownames(subtypes_colData) == colnames(subtypes_counts))

DE_subtypes <- differential_expression(
  coldata =  subtypes_colData,
  counts = subtypes_counts, comparison = 'subtype',
  add_title = 'threedatasets',
  batch = T)
DE_subtypesPE <- differential_expression(
  coldata =  subtypes_colData,
  counts = subtypes_counts, comparison = 'subtype', purityscaled = T,
  add_title = 'threedatasetsPE',
  batch = T)


### Stages
stages_colData <- coldata %>% 
  {.[.$batch == 'GSE217386' | .$batch == 'GSE143630',]} %>% 
  {rownames(.) <- NULL; .}
stages_counts <- counts %>% 
  {.[, c(stages_colData$run_accession)]} %>% as.matrix()
stages_colData <- stages_colData %>% 
  column_to_rownames(., var = 'run_accession') %>% 
  {.$run_accession <- rownames(.);.}
# Quality check, ensure columns and rows match in the same order
all(rownames(stages_colData) %in% colnames(stages_counts))
all(rownames(stages_colData) == colnames(stages_counts))

DE_stages <- differential_expression(coldata =  stages_colData, 
                         counts = stages_counts, comparison = 'stage',
                         add_title = 'GSE217386&GSE143630',
                         batch = T)
DE_stagesPE <- differential_expression(coldata =  stages_colData, 
                                     counts = stages_counts, comparison = 'stage',
                                     purityscaled = T,
                                     add_title = 'GSE217386&GSE143630PE',
                                     batch = T)


### Metastatic progression 
metastatic_colData <- coldata %>% 
  {.[.$batch == 'GSE143630',]} %>% 
  {rownames(.) <- NULL; .}
metastatic_counts <- counts %>% 
  {.[, c(metastatic_colData$run_accession)]} %>% as.matrix()
metastatic_colData <- metastatic_colData %>% 
  column_to_rownames(., var = 'run_accession') %>% 
  {.$run_accession <- rownames(.);.}
# Quality check, ensure columns and rows match in the same order
all(rownames(metastatic_colData) %in% colnames(metastatic_counts))
all(rownames(metastatic_colData) == colnames(metastatic_counts))

DE_metastatic <- differential_expression(
  coldata =  metastatic_colData,
  counts = metastatic_counts, comparison = 'metastatic')
DE_metastaticPE <- differential_expression(
  coldata =  metastatic_colData,
  counts = metastatic_counts, comparison = 'metastatic', purityscaled = T,
  add_title = 'PE')


### Treatment
treatment_colData <- coldata %>% 
  {.[.$batch == 'GSE153262',]} %>% 
  {rownames(.) <- NULL; .}
treatment_counts <- counts %>% 
  {.[, c(treatment_colData$run_accession)]} %>% as.matrix()
treatment_colData <- treatment_colData %>% 
  column_to_rownames(., var = 'run_accession') %>% 
  {.$run_accession <- rownames(.);.}
# Quality check, ensure columns and rows match in the same order
all(rownames(treatment_colData) %in% colnames(treatment_counts))
all(rownames(treatment_colData) == colnames(treatment_counts))

DE_treatment <- differential_expression(
  coldata =  treatment_colData, 
  counts = treatment_counts, comparison = 'treatment',
  batch = F)
DE_treatmentPE <- differential_expression(
  coldata =  treatment_colData, 
  counts = treatment_counts, comparison = 'treatment',
  batch = F, purityscaled = T,
  add_title = 'PE')


### Normal vs ccRCC GSE217386 & GSE69197 datasets
nVScc2dat_colData <- coldata %>% 
  {.[.$batch == 'GSE217386' | .$batch == 'GSE69197',]} %>% 
  {rownames(.) <- NULL; .}
nVScc2dat_counts <- counts %>% 
  {.[, c(nVScc2dat_colData$run_accession)]} %>% as.matrix()
nVScc2dat_colData <- nVScc2dat_colData %>% 
  column_to_rownames(., var = 'run_accession') %>% 
  {.$run_accession <- rownames(.);.}
# Quality check, ensure columns and rows match in the same order
all(rownames(nVScc2dat_colData) %in% colnames(nVScc2dat_counts))
all(rownames(nVScc2dat_colData) == colnames(nVScc2dat_counts))

DE_nVScc2dat <- differential_expression(
  coldata =  nVScc2dat_colData, 
  counts = nVScc2dat_counts, comparison = 'subtype',
  batch = T,
  add_title = 'nVScc2dat')
DE_nVScc2datPE <- differential_expression(
  coldata =  nVScc2dat_colData, 
  counts = nVScc2dat_counts, comparison = 'subtype',
  batch = T,
  add_title = 'nVScc2datPE', purityscaled = T)


### Subtypes GSE188486 dataset
subGSE188486_colData <- coldata %>% 
  {.[.$batch == 'GSE188486',]} %>% 
  {rownames(.) <- NULL; .}
subGSE188486_counts <- counts %>% 
  {.[, c(subGSE188486_colData$run_accession)]} %>% as.matrix()
subGSE188486_colData <- subGSE188486_colData %>% 
  column_to_rownames(., var = 'run_accession') %>% 
  {.$run_accession <- rownames(.);.}
# Quality check, ensure columns and rows match in the same order
all(rownames(subGSE188486_colData) %in% colnames(subGSE188486_counts))
all(rownames(subGSE188486_colData) == colnames(subGSE188486_counts))

DE_subGSE188486 <- differential_expression(
  coldata =  subGSE188486_colData, 
  counts = subGSE188486_counts, comparison = 'subtype',
  batch = F,
  add_title = 'subGSE188486')
DE_subGSE188486PE <- differential_expression(
  coldata =  subGSE188486_colData, 
  counts = subGSE188486_counts, comparison = 'subtype',
  batch = F,
  add_title = 'subGSE188486PE', purityscaled = T)


### Early vs Late stages GSE217386 dataset
eVSlGSE217386_colData <- coldata %>% 
  filter(batch == 'GSE217386', !is.na(stage),
         stage != 'normal') %>% 
  {rownames(.) <- NULL; .}
eVSlGSE217386_counts <- counts %>% 
  {.[, c(eVSlGSE217386_colData$run_accession)]} %>% as.matrix()
eVSlGSE217386_colData <- eVSlGSE217386_colData %>% 
  column_to_rownames(., var = 'run_accession') %>% 
  {.$run_accession <- rownames(.);.}
# Quality check, ensure columns and rows match in the same order
all(rownames(eVSlGSE217386_colData) %in% colnames(eVSlGSE217386_counts))
all(rownames(eVSlGSE217386_colData) == colnames(eVSlGSE217386_counts))

DE_eVSlGSE217386 <- differential_expression(
  coldata =  eVSlGSE217386_colData,
  counts = eVSlGSE217386_counts, comparison = 'early_late_stage',
  batch = F,
  add_title = 'eVSlGSE217386')
DE_eVSlGSE217386PE <- differential_expression(
  coldata =  eVSlGSE217386_colData,
  counts = eVSlGSE217386_counts, comparison = 'early_late_stage',
  batch = F,
  add_title = 'eVSlGSE217386PE', purityscaled = T)


### Normal vs ccRCC GSE217386 dataset
nVSccGSE217386_colData <- coldata %>% 
  {.[.$batch == 'GSE217386',]} %>% 
  {rownames(.) <- NULL; .}
nVSccGSE217386_counts <- counts %>% 
  {.[, c(nVSccGSE217386_colData$run_accession)]} %>% as.matrix()
nVSccGSE217386_colData <- nVSccGSE217386_colData %>% 
  column_to_rownames(., var = 'run_accession') %>% 
  {.$run_accession <- rownames(.);.}
# Quality check, ensure columns and rows match in the same order
all(rownames(nVSccGSE217386_colData) %in% colnames(nVSccGSE217386_counts))
all(rownames(nVSccGSE217386_colData) == colnames(nVSccGSE217386_counts))

DE_nVSccGSE217386 <- differential_expression(nVSccGSE217386_counts, 
                                             nVSccGSE217386_colData, 
                                             comparison = 'subtype',
                                             batch = F,
                                             purityscaled = F, 
                                             add_title = 'GSE217386',
                                             matched = F)
DE_nVSccGSE217386.matched <- differential_expression(nVSccGSE217386_counts, 
                                             nVSccGSE217386_colData, 
                                             comparison = 'subtype',
                                             batch = F,
                                             purityscaled = F, 
                                             add_title = 'GSE217386.matched',
                                             matched = T)
DE_nVSccGSE217386PE <- differential_expression(nVSccGSE217386_counts, 
                                             nVSccGSE217386_colData, 
                                             comparison = 'subtype',
                                             batch = F,
                                             purityscaled = T, 
                                             add_title = 'GSE217386PE',
                                             matched = F)
DE_nVSccGSE217386.matchedPE <- differential_expression(nVSccGSE217386_counts, 
                                               nVSccGSE217386_colData, 
                                               comparison = 'subtype',
                                               batch = F,
                                               purityscaled = T, 
                                               add_title = 'GSE217386.matchedPE',
                                               matched = T)


### Normal vs ccRCC GSE217386 dataset
nVSccGSE69197_colData <- coldata %>% 
  {.[.$batch == 'GSE69197',]} %>% 
  {rownames(.) <- NULL; .}
nVSccGSE69197_counts <- counts %>% 
  {.[, c(nVSccGSE69197_colData$run_accession)]} %>% as.matrix()
nVSccGSE69197_colData <- nVSccGSE69197_colData %>% 
  column_to_rownames(., var = 'run_accession') %>% 
  {.$run_accession <- rownames(.);.}
# Quality check, ensure columns and rows match in the same order
all(rownames(nVSccGSE69197_colData) %in% colnames(nVSccGSE69197_counts))
all(rownames(nVSccGSE69197_colData) == colnames(nVSccGSE69197_counts))

DE_nVSccGSE69197 <- differential_expression(nVSccGSE69197_counts, 
                                             nVSccGSE69197_colData, 
                                             comparison = 'subtype',
                                             batch = F,
                                             purityscaled = F, 
                                             add_title = 'GSE69197',
                                             matched = F)
DE_nVSccGSE69197.matched <- differential_expression(nVSccGSE69197_counts, 
                                            nVSccGSE69197_colData, 
                                            comparison = 'subtype',
                                            batch = F,
                                            purityscaled = F, 
                                            add_title = 'GSE69197.matched',
                                            matched = T)
DE_nVSccGSE69197PE <- differential_expression(nVSccGSE69197_counts, 
                                            nVSccGSE69197_colData, 
                                            comparison = 'subtype',
                                            batch = F, 
                                            add_title = 'GSE69197PE',
                                            matched = F, purityscaled = T)
DE_nVSccGSE69197.matchedPE <- differential_expression(nVSccGSE69197_counts, 
                                                    nVSccGSE69197_colData, 
                                                    comparison = 'subtype',
                                                    batch = F,
                                                    add_title = 'GSE69197.matchedPE',
                                                    matched = T, purityscaled = T)


### --------------------------------- Saving DE ---------------------------- ###
################################################################################

# Saving DE expression results in R file
all_DE_res <- c(
  DE_normalVStumour,
  DE_normalVStumourPE,
  DE_subtypes,
  DE_subtypesPE,
  DE_stages,
  DE_stagesPE,
  DE_metastatic,
  DE_metastaticPE,
  DE_treatment,
  DE_treatmentPE,
  DE_nVScc2dat,
  DE_nVScc2datPE,
  DE_subGSE188486,
  DE_subGSE188486PE,
  DE_eVSlGSE217386,
  DE_eVSlGSE217386PE,
  DE_nVSccGSE217386,
  DE_nVSccGSE217386.matched,
  DE_nVSccGSE217386PE,
  DE_nVSccGSE217386.matchedPE,
  DE_nVSccGSE69197,
  DE_nVSccGSE69197.matched,
  DE_nVSccGSE69197PE,
  DE_nVSccGSE69197.matchedPE
)

write_rds(all_DE_res, 'c:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/Figures_DE/all_DE_results.rds')
test <- read_rds('c:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/Figures_DE/all_DE_results.rds', refhook = NULL)


# Storing ACE2 lf2C and p.adj values in a csv file
ace2_expr <- data.frame()
gene_expression <- list()
for (tab in names(all_DE_res)){
  results <- all_DE_res[[tab]] %>% as.data.frame() %>% rownames_to_column(var = 'ensembl_gene_id_version')
  joined <- results %>% inner_join(gene_set.df, by = 'ensembl_gene_id_version')
  gene_expression[[tab]] <- joined
  ace2_padj <- joined[joined$external_gene_name == 'ACE2',] %>% dplyr::select(6)
  ace2_lf2c <- joined[joined$external_gene_name == 'ACE2',] %>% dplyr::select(3)
  cat(tab,', padj of ACE2:',as.character(ace2_lf2c), '\n')
  ace2_expr <- rbind(ace2_expr,data.frame(
    name = tab,
    padj = as.character(ace2_padj),
    lf2c = as.character(ace2_lf2c)
    ))
}

write.csv(ace2_expr, 'c:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/Figures_DE/ace2_expression_per_comparison.csv')
