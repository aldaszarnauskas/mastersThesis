# Author: Aldas Žarnauskas
# Title: Gene Set Enrichment Analysis with clusterProfiler
# Date: 14/006/2023
# Last time updated: 17/007/2023
# Objective: perform gene set enrichment analysis on a custom gene set list

# Clear environment
rm(list=ls(all.names = TRUE))

### ----------------------------- Libraries --------------------------------- ##
################################################################################
library(clusterProfiler)
library(msigdbr)
library(biomaRt)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggstance)
library(forcats)


### ----------------------------- Data Import ------------------------------- ##
################################################################################

# load differential expression data
all_DE_res <- readRDS('C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/DE_plots/all_DE_results.rds')


### ---------------------------------- GSEA --------------------------------- ##
################################################################################
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_set_enrichment <- function(DE){
  
  # Select LF2C column
  gene.list <- threedatasetsPEchRCC_vs_normal$log2FoldChange
  
  # Retrieve corresponding external_gene_name and ensembl_gene_id of the gene.list
  gene_symbols <- getBM(attributes = c("external_gene_name", 'ensembl_gene_id', 'ensembl_gene_id_version'),
                        filters = "ensembl_gene_id_version",
                        values = names(gene.list),
                        mart = ensembl)
  
  # Remove duplicates
  gene_symbols <- gene_symbols[!(gene_symbols$external_gene_name == ''),]
  
  
  # Merge gene.list and gene_symbols by ensembl_gene_id_version
  gene_list.df <- gene.list %>% as.data.frame() %>% 
    {colnames(.) <- 'log2FoldChange'; .} %>% 
    {rownames_to_column(., var = 'ensembl_gene_id_version')} %>% 
    {merge(., gene_symbols, by = 'ensembl_gene_id_version')}
  
  
  
  # Remove duplicates from gene_list.df
  if (any(duplicated(gene_list.df$external_gene_name))){
    
    dup_gene_symbols <- gene_list.df %>%
      dplyr::filter(duplicated(external_gene_name)) %>%
      dplyr::pull(external_gene_name)
    
    
    gene_list.df %>%
      dplyr::filter(external_gene_name %in% dup_gene_symbols) %>%
      dplyr::arrange(external_gene_name)
    
    filt_gene.df <- gene_list.df %>%
      dplyr::arrange(dplyr::desc(log2FoldChange)) %>%
      # Filter out the duplicated rows using `dplyr::distinct()`
      dplyr::distinct(external_gene_name, .keep_all = TRUE)
    
    gene_list.df <- filt_gene.df %>% 
      {. <- .$log2FoldChange} %>% 
      {names(.) <- filt_gene.df$external_gene_name; .}
    
  }
  
  # ensure gene_list is a named numeric vector
  if (!is.numeric(gene_list.df) || is.null(names(gene_list.df))) {
    stop("gene_list must be a named numeric vector.")
  }
  
  # check for NA values
  if (any(is.na(gene_list.df))) {
    # remove NA values
    gene_list.df <- gene_list.df[!is.na(gene_list.df)] %>%
      # remove NaN values
      .[!is.nan(.)]
  }
  
  # Construct GSEA Input -------------------------------------------------------
  
  # Obtain MSigDB pathways
  msigdb <- msigdbr(species = "Homo sapiens")
  
  # obtain msigdb
  gda_custom <- data.frame(
    geneId = msigdb$human_ensembl_gene,
    geneSymbol = msigdb$human_gene_symbol,
    pathwayId = msigdb$gs_id,
    pathwayName = msigdb$gs_name)
  
  # Select pathways of interest
  pathway_types <- c("KEGG_")
  
  # filter pathways
  filtered_gda_custom <- gda_custom[
    grepl(paste(pathway_types, collapse = "|"), gda_custom$pathwayName), ]
  
  # Construct the term objects as GSEA input
  term2gene_custom=filtered_gda_custom[, c("pathwayId", "geneSymbol")]
  term2name_custom=filtered_gda_custom[, c("pathwayId", "pathwayName")]
  
  
  # Perform GSEA ---------------------------------------------------------------
  gse <-  GSEA(
    gene_list.df, TERM2GENE=term2gene_custom, TERM2NAME=term2name_custom, 
    pvalueCutoff = 1) 
  
  return(gse)
  
}

### ----------------------------- Running GSEA ------------------------------ ##
################################################################################

comparisons_of_interest <- c(
  "threedatasetsPEchRCC_vs_normal", 
  "threedatasetsPEccRCC_vs_normal",
  "threedatasetsPEpRCC_vs_normal",
  "threedatasetsPEccRCC_vs_chRCC",
  "threedatasetsPEpRCC_vs_chRCC",
  "threedatasetsPEpRCC_vs_ccRCC")

threedatasetsPEchRCC_vs_normal <- all_DE_res[["threedatasetsPEchRCC_vs_normal"]]
threedatasetsPEccRCC_vs_normal <- all_DE_res[["threedatasetsPEccRCC_vs_normal"]]
threedatasetsPEpRCC_vs_normal <- all_DE_res[["threedatasetsPEpRCC_vs_normal"]]
threedatasetsPEccRCC_vs_chRCC <- all_DE_res[["threedatasetsPEccRCC_vs_chRCC"]]
threedatasetsPEpRCC_vs_chRCC <- all_DE_res[["threedatasetsPEpRCC_vs_chRCC"]]
threedatasetsPEpRCC_vs_ccRCC <- all_DE_res[["threedatasetsPEpRCC_vs_ccRCC"]]


GSEA_threedatasetsPEchRCC_vs_normal <- gene_set_enrichment(threedatasetsPEchRCC_vs_normal)
GSEA_threedatasetsPEccRCC_vs_normal <- gene_set_enrichment(threedatasetsPEccRCC_vs_normal)
GSEA_threedatasetsPEpRCC_vs_normal <- gene_set_enrichment(threedatasetsPEpRCC_vs_normal)
GSEA_threedatasetsPEccRCC_vs_chRCC <- gene_set_enrichment(threedatasetsPEccRCC_vs_chRCC)
GSEA_threedatasetsPEpRCC_vs_chRCC <- gene_set_enrichment(threedatasetsPEpRCC_vs_chRCC)
GSEA_threedatasetsPEpRCC_vs_ccRCC <- gene_set_enrichment(threedatasetsPEpRCC_vs_ccRCC)


### ----------------------------- GSEA Visualisation ---=-------------------- ##
################################################################################

gsea_plots <- function(gsea_results, pathNR){
  
  # order gsea results based on abs(NES)
  gse.df <- gsea_results %>% as.data.frame() %>% 
    {. <- .[order(abs(.$NES), decreasing =T), ]; .} # ordered on absolute NES
  
  # Select top pathways
  gse.top <- gse.df[1:pathNR, ]
  
  # Select IDs and desciptions of the top pathways
  gse.topIDs <- gse.top$ID
  gse.topDescription <- gse.top$Description
  
  # order gsea results object by abs(NES)
  y <- mutate(gsea_results, ordering = abs(NES)) %>%
    arrange(desc(ordering)) 
  
  # filter the y object to contain only the top pathways
  y_filtered <- y %>%
    filter(Description %in% c(gse.topDescription))
  
  
  # Plot the inverted barplot
  ggplot(y_filtered, aes(
    NES, fct_reorder(Description, NES),
    fill = ifelse(NES > 0 & p.adjust < 0.05, "Significant Positively Enriched",
                  ifelse(
                    NES < 0 & p.adjust < 0.05, 
                    "Significant Negatively Enriched", "Not Significant"))),
    showCategory=(pathNR)) + 
    geom_barh(stat = 'identity', width = 0.6) +
    scale_fill_manual(
      values = c(
        "Significant Positively Enriched" = "#E41A1C", 
        "Significant Negatively Enriched" = "#377EB8", 
        "Not Significant" = "grey")) +
    ylab(NULL) + 
    guides(fill = guide_legend(title = "Enrichment & P-adjusted")) +
    theme_classic() +
    theme(legend.position = "bottom",
          axis.ticks.y = element_blank()) +
    labs(caption = paste('
                         '))
  
  
  ggsave(paste(substitute(gsea_results), 'pathway_enrichment', '.png', sep = ''),
         dpi = 1000,
         width = 10, height = 10)
  
}

gsea_plots(GSEA_threedatasetsPEchRCC_vs_normal, 20)
gsea_plots(GSEA_threedatasetsPEccRCC_vs_normal, 10)
gsea_plots(GSEA_threedatasetsPEpRCC_vs_normal, 10)
gsea_plots(GSEA_threedatasetsPEccRCC_vs_chRCC, 10)
gsea_plots(GSEA_threedatasetsPEpRCC_vs_chRCC, 10)
gsea_plots(GSEA_threedatasetsPEpRCC_vs_ccRCC, 10)
