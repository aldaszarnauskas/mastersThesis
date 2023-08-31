# Author: Aldas Žarnauskas
# Title: Visualisation Pipelines
# Date: 10/007/2023
# Last time updated: 26/008/2023
# Objective: Visualise DE results

# Clear environment
rm(list=ls(all.names = TRUE)) 

### ----------------------------- Libraries --------------------------------- ##
################################################################################
library(pheatmap)
library(tibble)
library(ggplot2)
library(ggrepel)
library(tidyverse)


### ----------------------------- Data Import ------------------------------- ##
################################################################################

all_DE_res <- readRDS('C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/Figures_DE/all_DE_results.rds')
ras_raas_genes <- read.csv('C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/ras_raas_genes.csv', header = T)


### ----------------------------- Volcano Plot ------------------------------ ##
################################################################################

volcano_plot <- function(DE, Type, caption = ''){
  
  # Set rownames as ensembl_gene_id_version column
  results <- DE[[Type]] %>% as.data.frame() %>% 
    {rownames_to_column(., var = 'ensembl_gene_id_version')}
  
  # Set significance thresholds
  log2FC_threshold <- 1
  pvalue_threshold <- -log10(0.05)
  
  
  # Define genes of interest
  genes_of_interest <- c('ACE2', 'ACE', 'REN', 'AGT', 'AGTR1', 'AGTR2', 
                         'CYP11B2', 'MAS1')
  ras_raas_genes_top <- subset(
    ras_raas_genes, external_gene_name %in% genes_of_interest)
  ras_raas_volcano <- filter(
    results, ensembl_gene_id_version %in% 
      ras_raas_genes_top$ensembl_gene_id_version) %>% 
    base::merge(ras_raas_genes_top, by = 'ensembl_gene_id_version')
  
  ras_genes <- results %>%
    filter(ensembl_gene_id_version %in% ras_raas_volcano$ensembl_gene_id_version)
  
  p <- ggplot(results, aes(x = log2FoldChange, y = -log10(padj))) +
    theme_bw()  +
    labs(x = "Log2(fold change)", y = "-log10(adjusted P-value)")  +
    geom_point(aes(color = ifelse(
      ensembl_gene_id_version %in% ras_genes$ensembl_gene_id_version, 
      "Genes of interest",
      ifelse(
        log2FoldChange >= log2FC_threshold & -log10(padj) >= pvalue_threshold, 
        "Upregulated", 
        ifelse(log2FoldChange <= -log2FC_threshold & -log10(padj) >= pvalue_threshold, 
               "Downregulated",
               ifelse(abs(log2FoldChange) <= log2FC_threshold & -log10(padj) >= pvalue_threshold, 
                      'Not significant', 
                      ifelse(
                        abs(log2FoldChange) >= log2FC_threshold & -log10(padj) <= pvalue_threshold, 
                        "Not significant", "No change"))))))) +
    scale_color_manual(values = c("Genes of interest" = "#984EA3", 
                                  "Upregulated" = "#E41A1C", 
                                  "Downregulated" = "#377EB8",
                                  "Not significant" = '#4DAF4A',
                                  'No change' = 'grey')) +
    geom_point(data = ras_genes,
               shape = 21,
               size = 2, 
               fill = "#984EA3", 
               colour = "black") +
    theme(
      legend.position = "right",
      legend.title = element_text(
        colour = "black", face = "bold.italic", family = "Helvetica"),
      legend.text = element_text(
        face = "italic", colour = "black", family = "Helvetica"),
      axis.title = element_text(
        family = "Helvetica", size = (10), colour = "black"),
      axis.text = element_text(
        family = "Courier", colour = "black", size = (10)),
      plot.caption= element_text(size=7, face = 'italic',
                                 color="Black"),
      panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    ) +
    plot_annotation(
      caption = paste(caption, '
                      Volcano Plot of ', Type)) +
    guides(color = guide_legend(title = "Expression Change")) +
    scale_x_continuous(breaks = c(seq(-10, 10, 2)),     
                       limits = c(-10, 10)) + 
    
    geom_hline(yintercept = pvalue_threshold, linetype = "dashed", 
               color = "black") +
    geom_vline(xintercept = c(log2FC_threshold, -log2FC_threshold), 
               linetype = "dashed", color = "black") +
    geom_label_repel(data = ras_raas_volcano,
                     aes(label = external_gene_name),
                     force = 2,
                     nudge_y = 1) 
  
  return(p)
  
}

# Generate Volcano Plots for all DE in al_DE_res
DE_names <- names(all_DE_res)

for (i in DE_names){
  volcano <- volcano_plot(all_DE_res, i, caption = i)
  ggsave(paste(i, '.png', sep = ''),
         dpi = 1000,
         width = 10, height = 10)
}

