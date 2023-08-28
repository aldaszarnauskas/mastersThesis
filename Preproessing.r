# Author: Aldas Žarnauskas
# Title: Normalisation, PCA & Normalised Heatplot
# Date: 02/008/2023
# Last time updated: 26/008/2023
# Objective: Visualise and normalise raw counts of rcc samples
# Overview: In this script there are defined four functions. One function performs
# required normalisations on raw counts and saves it in a csv file. The second
# function plots PCA plots of the normalised data. The third function plots the
# Heatplot of z-scores of normalised counts. And the fourth function plots a 
# boxplot of ACE2 normalised expression. Additional function was created to
# combine all functions mentioned above. In addition, the latter function 
# generates two sets of plots, one set represents data without tumour purity
# estimation and the other set represents data with tumour purity estimation.

#Clear environment
rm(list=ls()) 

### ----------------------------- Libraries --------------------------------- ##
################################################################################
library(DESeq2)
library(dplyr)
library(pheatmap)
library(tibble)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(grid)
library(gridExtra)
library(patchwork)
library(ggplotify)
library(RColorBrewer)

### ----------------------------- Data Import ------------------------------- ##
################################################################################

ras_raas_genes <- read.csv('C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/ras_raas_genes.csv', 
                           header = T)
ace2_ensembl <- ras_raas_genes[
  ras_raas_genes$external_gene_name == 'ACE2',]$ensembl_gene_id_version


counts <- read.csv("C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/RCC_counts_v2.csv",
                   header = T, sep = ',') %>% 
  .[, -c(which(colnames(.) == 'SRR10883885.1'))]

# Replace NA values in counts with 0, if necessary
any(is.na(counts))
counts[is.na(counts)] <- 0
any(is.na(counts))

# Upload metadata
coldata <- read.csv(
  "C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/colData_kidney_cancer_tumourPurity.csv")



### ------------------------ Normalisation+PCA+Heaplot ---------------------- ##
################################################################################

# batchRem_normalised takes in counts as input and returns variance stabilised
# transformed counts. In addition, it also performs batch normalisation, if 
# necessary.
batchRem_normalised <- function(counts, coldata, comparison, batch = F,
                                purityscaled = T,
                                matched = F){
  
  # If normal is among the features, it sets as the reference level
  if ('normal' %in% unique(coldata[[comparison]])){
    
    ref_val <- 'normal'
    comparisons <- unique(coldata[[comparison]]) %>% 
      .[. != ref_val] %>%
      c(ref_val, .)
    
  } else {
    
    comparisons <- unique(coldata[[comparison]]) %>% sort()
    
  }
  
  # Based on preferences, the following defines DESeq formula
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
    
    deseq_formula <- as.formula(paste("~ purityscaled + matched + ", comparison))
    
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
  
  
  # Removes batch effects if batch=TRUE
  if (batch){
    
    vsd <- vst(dds, blind=F)
    vsd_prior_batch_removal <- vsd
    mat <- assay(vsd)
    mm <- model.matrix(as.formula(paste('~ ', comparison)), colData(vsd))
    mat <- limma::removeBatchEffect(mat, batch=vsd$batch, design=mm)
    assay(vsd) <- mat
    all_counts <- mat
    
    result <- list(vsd_prior_batch_removal = vsd_prior_batch_removal,
                   vsd = vsd,
                   all_counts = all_counts)
      
  } else {
    
    vsd <- vst(dds, blind=F)
    
    result <- list(vsd = vsd,
                   all_counts = assay(vsd))
    
  }
  
  return(result)
  
}

# pca1 plots variance stabilised couts if batch effects were not removed
pca1 <- function(vsd_prior_batch_removal, title, comparison, caption_name = ''){
  
  
  p1 <- plotPCA(vsd_prior_batch_removal, intgroup = comparison) +
    theme_classic() +
    geom_point(aes(color = eval(as.symbol(comparison))),
               size = 2) +  
    theme(legend.position = 'right',
          legend.title = element_text(
            colour = "black", face = "bold.italic", family = "Helvetica"),
          legend.text = element_text(
            face = "italic", colour = "black", family = "Helvetica"),
          axis.title = element_text(
            family = "Helvetica", size = (10), colour = "black"),
          axis.text = element_text(
            family = "Courier", colour = "black", size = (10)),
          plot.caption= element_text(size=7, face = 'italic',
                                     color="Black")) +
    labs(caption = paste(caption_name, '
                         ', title),
         color = comparison) + 
    scale_fill_brewer(palette = "Set1", aesthetics = c("colour", "fill"))
  
  return(p1)
  
}

# pca2 plots variance stabilised couts if batch effects were removed
pca2 <- function(vsd_prior_batch_removal, vsd, title, comparison, 
                 caption = caption){
  
  # PCA prior batch removal
  p1 <- plotPCA(vsd_prior_batch_removal, intgroup = comparison) +
    theme_classic() +
    geom_point(aes(color = eval(as.symbol(comparison))), size = 2) +  
    theme(legend.position = 'none',
          plot.title = element_text(hjust = 0.5,
                                    family = "Helvetica", 
                                    face = "bold", size = (15), ),
          axis.title = element_text(
            family = "Helvetica", size = (10), colour = "black"),
          axis.text = element_text(
            family = "Courier", colour = "black", size = (10))) + 
    ggtitle('Prior batch removal') +
    scale_fill_brewer(palette = "Set1", aesthetics = c("colour", "fill")) 
  
  
  p1_legend <- p1 + theme(legend.position = 'top',
                          legend.title = element_text(
                            colour = "black", face = "bold.italic", 
                            family = "Helvetica"),
                          legend.text = element_text(
                            face = "italic", colour = "black", 
                            family = "Helvetica")) + 
    labs(color = comparison) 
  
  
  
  
  
  # PCA after batch removal
  p2 <- plotPCA(vsd, intgroup = comparison) +
    theme_classic() +
    geom_point(aes(color = eval(as.symbol(comparison))), size = 2) +  
    theme(legend.position = 'none',
          plot.title = element_text(hjust = 0.5,
                                    family = "Helvetica", 
                                    face = "bold", size = (15), ),
          axis.title = element_text(
            family = "Helvetica", size = (10), colour = "black"),
          axis.text = element_text(
            family = "Courier", colour = "black", size = (10))) + 
    ggtitle('After batch removal') +
    scale_fill_brewer(palette = "Set1", aesthetics = c("colour", "fill"))
  
  # Extract the legend from plot1
  legend_plot <- cowplot::get_legend(p1_legend)

  
  caption <- textGrob(
    label = paste(caption, '
                  ', title),
    gp = gpar(fontface = "italic", fontsize = 7)
  )
  
  # Combine PCA1, PCA2 and the legend of PCA plots
  combined_plots <- grid.arrange(arrangeGrob(p1,p2, ncol=2,nrow = 1, widths = c(0.5,0.5)),
                                 arrangeGrob(legend_plot, ncol=1, nrow=1), 
                                 heights=c(0.9,0.1), nrow = 2,
                                 bottom = caption)
  
  return(combined_plots)
  
}

# nornalised_heatmap plots the heatmap of z-score normalised counts
nornalised_heatmap <- function(normalised_counts, coldata, comparison,
                               genes_of_interest, title,
                               cutree_rows = 1, cutree_cols = 1,
                               clustering_distance_cols = 'euclidean',
                               clustering_distance_rows = 'euclidean',
                               caption = ''){
  
  # Selecting the counts of the genes in the RAAS pathway 
  df <- normalised_counts %>% 
    as.data.frame() %>% 
    t(.) %>% scale(.) %>% t(.) %>% as.data.frame() %>% 
    rownames_to_column(., var = 'ensembl_gene_id_version') %>% 
    merge(genes_of_interest[, c(1, 2)], by = 'ensembl_gene_id_version') %>% 
    column_to_rownames(., 'external_gene_name') %>% .[, 2:ncol(.)]
  
  # Creating annotation dataframe for the heatmap
  annotation <- as.data.frame(factor(coldata[[comparison]]))
  colnames(annotation) <- comparison
  annotation$batch <- coldata$batch
  rownames(annotation) <- rownames(coldata)
  
  # Replace NA values in counts with 0, if necessary
  any(is.na(df))
  df[is.na(df)] <- 0
  any(is.na(df))
  
  max_val <- max(abs(df))
  
  
  # Define symmetric breaks
  my_breaks <- seq(-max_val, max_val, length.out = 101)
  
  heatmap <- pheatmap(
    df, 
    annotation_col = annotation,
    cluster_rows = T,
    cluster_cols = T,
    color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
    breaks = my_breaks,
    legend = T,
    clustering_distance_cols = clustering_distance_cols,
    clustering_distance_rows = clustering_distance_rows,
    cutree_rows = cutree_rows,
    cutree_cols = cutree_cols,
    fontsize_row = 8,  
    fontsize_col = 8
  )
  
  
  combined_plot <- as.ggplot(heatmap) 
  combined_plot <- combined_plot + plot_annotation(
    caption = paste('Normalised Z-Score Heatplot of ', title, '
                    ', caption, sep = '')
    )
  
  return(combined_plot)
  
}


plot_gene <- function(counts, coldata, comparison, 
                      gene_symbol, gene_ensembl, title,
                      caption = ''){
  
  # Obtains expression of gene of interes
  gene_expression.df <- counts[rownames(counts) == gene_ensembl,] %>%
    as.data.frame() %>% 
    {colnames(.) <- 'gene_expression'; .} %>% 
    rownames_to_column(var = 'run_accession') %>% 
    merge(coldata, by = 'run_accession') %>% 
    {. <- .[order(.[[comparison]]),]} %>% 
    {.$run_accession <- factor(
      .$run_accession, levels = unique(.$run_accession[order(.[[comparison]])])); .}
    
  if ("normal" %in% gene_expression.df[[comparison]]){
    
    gene_expression.df[[comparison]] <- factor(
      gene_expression.df[[comparison]], 
      levels = c("normal", 
                 sort(unique(gene_expression.df[[comparison]])) %>% 
                   .[-which(. == "normal")]
                 ))
    
  }
  
  
  boxplot <- ggplot(gene_expression.df,
         aes(x = eval(as.symbol(comparison)), y = gene_expression, 
             fill = eval(as.symbol(comparison)))) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA,
                 show.legend = FALSE) +
    geom_jitter(aes(color = eval(as.symbol(comparison))), width = 0.3, 
                alpha = 0.5) +
    scale_fill_brewer(palette = "Set1", aesthetics = c("colour", "fill")) +
    labs(x = comparison, y = paste("Normalized ", gene_symbol, " Expression", 
                                   sep = '')) +
    theme(
      axis.title = element_text(
        family = "Helvetica", size = (10), colour = "black"),
      axis.text = element_text(
        family = "Courier", colour = "black", size = (10)),
      plot.caption= element_text(size=7, face = 'italic',
                                 color="Black"),
      legend.position = 'none'
    ) + 
    theme_classic()

  combined_plot <- boxplot + plot_annotation(
    caption = paste(caption,"
                    Normalized Expression of", gene_symbol, "in", title)
  ) +
    guides(
      color = F,
      fill = F,
      shape = F)
  
  return(combined_plot)

}


# normalised_counts_plots combines the functions above
normalised_counts_plots <- function(all_colData, all_counts,
                                    all_comparison, title, title_purity,
                                    batch = T,
                                    cut_rows = 1, cut_cols = 1,
                                    clustering_distance_cols = 'euclidean',
                                    clustering_distance_rows = 'euclidean',
                                    matched = F,
                                    caption = '',
                                    filename = '',
                                    tumourEstimated = F){
  
  # Prior tumour purity estimation
  all_results <- batchRem_normalised(
    all_counts, all_colData, batch = batch, 
    purityscaled = F, comparison = all_comparison,
    matched = matched)
  
  all_vsd_prior <- all_results$vsd_prior_batch_removal
  all_vsd <- all_results$vsd
  all_normalised <- all_results$all_counts
  write.csv(all_normalised, file = paste('normalised_counts', filename,
                                         '.csv', sep = ''), 
            row.names = T)
  
  if (batch){
    
    pca_plot <- pca2(all_vsd_prior, all_vsd, title, caption = caption,
                     all_comparison)
    
  } else {
    
    pca_plot <- pca1(all_vsd, title, all_comparison, caption_name = caption)
    
  }
  norm_heat <- nornalised_heatmap(all_normalised, all_colData, comparison = all_comparison,
                                  ras_raas_genes, title,
                                  cutree_rows = cut_rows, 
                                  cutree_cols = cut_cols,
                                  clustering_distance_cols = clustering_distance_cols,
                                  clustering_distance_rows = clustering_distance_rows,
                                  caption = caption)
  
  ace2_boxplot <- plot_gene(counts = all_normalised, 
                            coldata = all_colData, 
                            comparison = all_comparison, 
                            gene_symbol = 'ACE2', 
                            gene_ensembl = ace2_ensembl, 
                            title = title,
                            caption = caption)
  
  
  if (batch){
    
    ggsave(paste('pca_plot_of_', filename, sep = '', '.png'), pca_plot,
           dpi = 1000,
           width = 10, height = 5,)
    
  } else{
    
    ggsave(paste('pca_plot_of_', filename, sep = '', '.png'), pca_plot,
           dpi = 1000,
           width = 5, height = 5,)
    
  }
  
  ggsave(paste('normalised_heatmap_of_', filename, sep = '', '.png'), norm_heat,
         dpi = 1000,
         width = 10, height = 10)
  ggsave(paste('normalised_count_boxplot_of_', filename, sep = '', '.png'), ace2_boxplot,
         dpi = 1000,
         width = 5, height = 5,)
  
  
  if (tumourEstimated){
    
    # Prior tumour purity estimation
    caption <- paste(caption, ' with ', 'Tumour purity estimated', sep = '')
    filename <- paste(filename, 'purityEstimated', sep = '')
    
    all_results_purity <- batchRem_normalised(
      all_counts, all_colData, batch = batch,
      purityscaled = T, comparison = all_comparison)
    
    
    all_vsd_prior_purity <- all_results_purity$vsd_prior_batch_removal
    all_vsd_purity <- all_results_purity$vsd
    all_normalised_purity <- all_results_purity$all_counts
    write.csv(all_normalised_purity, file = paste('normalised_counts', filename,
                                                  '.csv', sep = ''),
              row.names = T)
    
    if (batch){
      
      pca_plot <- pca2(all_vsd_prior_purity, all_vsd_purity, title,
                       caption = caption, all_comparison)
      
    } else {
      
      pca_plot <- pca1(all_vsd_purity, title, caption = caption, all_comparison)
      
    }
    
    norm_heat <- nornalised_heatmap(all_normalised_purity, all_colData, comparison = all_comparison,
                                    ras_raas_genes, title,
                                    cutree_rows = cut_rows,
                                    cutree_cols = cut_cols, caption = caption,
                                    clustering_distance_cols = clustering_distance_cols,
                                    clustering_distance_rows = clustering_distance_rows)
    ace2_boxplot <- plot_gene(counts = all_normalised_purity,
                              coldata = all_colData,
                              comparison = all_comparison,
                              gene_symbol = 'ACE2',
                              gene_ensembl = ace2_ensembl,
                              title = title, caption = caption)
    
    if (batch){
      
      ggsave(paste('pca_plot_of_', filename, sep = '', '.png'), pca_plot,
             dpi = 1000,
             width = 10, height = 5,)
      
    } else{
      
      ggsave(paste('pca_plot_of_', filename, sep = '', '.png'), pca_plot,
             dpi = 1000,
             width = 5, height = 5,)
      
    }
    
    ggsave(paste('normalised_heatmap_of_', filename, sep = '', '.png'), norm_heat,
           dpi = 1000,
           width = 10, height = 10)
    ggsave(paste('normalised_count_boxplot_of_', filename, sep = '', '.png'), ace2_boxplot,
           dpi = 1000,
           width = 5, height = 5,)
    
  }
  
}


### ------------------------ Datasets Selection & Plotting ------------------ ##
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

normalised_counts_plots(all_colData <- normalVScancer_colData,
                        all_counts <- normalVScancer_counts,
                        all_comparison <- 'subtype',
                        title <- 'Normal vs RCC Tumour',
                        caption = 'Datasets: GSE188486, GSE217386, GSE69197',
                        filename = 'normalVScancer_colData',
                        tumourEstimated = T
                        )


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

normalised_counts_plots(all_colData <- subtypes_colData,
                        all_counts <- subtypes_counts,
                        all_comparison <- 'subtype',
                        title <- 'RCC Subtypes',
                        # clustering_distance_cols = 'correlation',
                        # clustering_distance_rows = 'correlation',
                        # cut_rows = 2,
                        # cut_cols = 4,
                        caption = 'Datasets: GSE188486, GSE217386, GSE69197
                        Patients 6801 & 3596 were removed as outliers',
                        filename = 'subtypesGSE188486&GSE217386&GSE69197',
                        tumourEstimated = T)



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

normalised_counts_plots(all_colData <- nVScc2dat_colData,
                        all_counts <- nVScc2dat_counts,
                        all_comparison <- 'subtype',
                        title <- 'Normal vs ccRCC',
                        batch = T,
                        matched = F,
                        caption = 'Datasets: GSE217386, GSE69197',
                        filename = 'nVSccRCCGSE217386&GSE69197',
                        tumourEstimated = T)


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

normalised_counts_plots(all_colData <- subGSE188486_colData,
                        all_counts <- subGSE188486_counts,
                        all_comparison <- 'subtype',
                        title <- 'RCC Subtypes',
                        batch = F,
                        caption = 'Datasets: GSE188486',
                        filename = 'subtypesGSE188486',
                        tumourEstimated = T)



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

normalised_counts_plots(all_colData <- stages_colData,
                        all_counts <- stages_counts,
                        all_comparison <- 'stage',
                        title <- 'RCC Stages',
                        batch = T,
                        caption = 'Datasets: GSE217386, GSE143630',
                        filename = 'stagesGSE217386&GSE143630',
                        tumourEstimated = T)


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

normalised_counts_plots(all_colData <- treatment_colData,
                        all_counts <- treatment_counts,
                        all_comparison <- 'treatment',
                        title <- 'RCC Treatment',
                        batch = F,
                        caption = 'Datasets: GSE153262',
                        filename = 'treatmentGSE153262',
                        tumourEstimated = T
                        )

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

normalised_counts_plots(all_colData <- metastatic_colData,
                        all_counts <- metastatic_counts,
                        all_comparison <- 'metastatic',
                        title <- 'RCC Metastatic Progression',
                        batch = F,
                        caption = 'Datasets: GSE143630',
                        filename = 'metastaticGSE143630',
                        tumourEstimated = T
)



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

normalised_counts_plots(all_colData <- eVSlGSE217386_colData,
                        all_counts <- eVSlGSE217386_counts,
                        all_comparison <- 'early_late_stage',
                        title <- 'Early vs Late',
                        batch = F,
                        caption = 'Datasets: GSE217386',
                        filename = 'eVSlGSE217386',
                        tumourEstimated = T
)



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

normalised_counts_plots(all_colData <- nVSccGSE217386_colData,
                        all_counts <- nVSccGSE217386_counts,
                        all_comparison <- 'subtype',
                        title <- 'Normal vs ccRCC',
                        batch = F,
                        matched = F,
                        caption = 'Datasets: GSE217386',
                        filename = 'nVSccGSE217386',
                        tumourEstimated = T
)




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

normalised_counts_plots(all_colData <- nVSccGSE69197_colData,
                        all_counts <- nVSccGSE69197_counts,
                        all_comparison <- 'subtype',
                        title <- 'Normal vs ccRCC GSE69197 dataset',
                        title_purity <- 'Normal vs ccRCC GSE69197 dataset, Tumour Purity Excluded',
                        batch = F,
                        matched = F,
                        caption = 'Datasets: GSE69197',
                        filename = 'nVSccGSE69197',
                        tumourEstimated = T
)
 
