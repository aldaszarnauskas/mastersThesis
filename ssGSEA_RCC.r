# Author: Aldas Žarnauskas
# Title: Single Sample Gene Set Enrichment Analysis (ssGSEA)
# Date: 22/008/2023
# Last time updated: 22/008/2023
# Objective: Perform single sample gene set enrichment analysis (ssGSEA)

#Clear environment
rm(list=ls(all.names = TRUE))

###-------------------------------- Library ---------------------------------###
################################################################################
library(GSVA)
library(msigdbr)
library(tibble)
library(Biobase)
library(GSEABase)
library(dplyr)
library(org.Hs.eg.db)
library(pheatmap)
library(biomaRt)
library(ggplotify)
library(patchwork)
library(ggplot2)
library(RColorBrewer)

###-------------------------------- Data Import -----------------------------###
################################################################################
# Tissue data
counts <- read.csv("C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/Norm_PCA_BOX_Heat/normalised_countssubtypesGSE188486&GSE217386&GSE69197purityEstimated.csv") %>% 
  column_to_rownames(var = 'X')
coldata <- read.csv("C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/colData_kidney_cancer_tumourPurity.csv") %>% 
  filter(run_accession %in% colnames(counts)) %>% 
  column_to_rownames(var = 'run_accession')
coldata.annotated <- AnnotatedDataFrame(coldata) 

counts.ES <- ExpressionSet(counts %>% as.matrix(),
                           phenoData = coldata.annotated)


# Connect to ensembl BioMart database
ensembl <- useEnsembl("ensembl")

# Select a dataset, here homo sapiens, and update mart object using useDataset
ensembl <- useDataset('hsapiens_gene_ensembl', mart = ensembl)


###------------------------------------ MSigDB ------------------------------###
################################################################################
# Get a list of gene sets for human species
msigdb <- msigdbr(species = "Homo sapiens")

# Define hallmark pathways
hallmark_pathways.df <- msigdb[grepl('HALLMARK_', msigdb$gs_name),] %>% 
  merge(., (
    
    getBM(
      attributes = c("ensembl_gene_id",'ensembl_gene_id_version'),
      filters = "ensembl_gene_id",
      values = .$ensembl_gene,
      mart = ensembl)
    
    ),
    by.x = "ensembl_gene", by.y = "ensembl_gene_id"
    )

hallmark_pathways.list <- split(
  x = hallmark_pathways.df$ensembl_gene_id, 
  f = hallmark_pathways.df$gs_name
)


# RAAS pathway 
raas.list <- read.csv('C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/ras_raas_genes.csv') 
raas.list <- list(
  "Renin_Angiotensin_Aldosterone_Pathway" = raas.list$ensembl_gene_id_version)

# HALLMARK_ANGIOGENESIS pathway 
HALLMARK_ANGIOGENESIS.list <- list(
    "HALLMARK_ANGIOGENESIS" = hallmark_pathways.list$HALLMARK_ANGIOGENESIS)

# HALLMARK_G2M_CHECKPOINT pathway 
HALLMARK_G2M_CHECKPOINT.list <- list(
  "HALLMARK_G2M_CHECKPOINT" = hallmark_pathways.list$HALLMARK_G2M_CHECKPOINT)

# HALLMARK_INFLAMMATORY_RESPONSE pathway 
HALLMARK_INFLAMMATORY_RESPONSE.list <- list(
  "HALLMARK_INFLAMMATORY_RESPONSE" = hallmark_pathways.list$HALLMARK_INFLAMMATORY_RESPONSE)

# HALLMARK_HYPOXIA pathway 
HALLMARK_HYPOXIA.list <- list(
  "HALLMARK_HYPOXIA" = hallmark_pathways.list$HALLMARK_HYPOXIA)


###------------------------------------ ssGSEA ------------------------------###
################################################################################

gsva_halmark <- gsva(counts.ES, hallmark_pathways.list)
gsva_ANGIOGENESIS <- gsva(counts.ES, HALLMARK_ANGIOGENESIS.list)
gsva_G2M_CHECKPOINT <- gsva(counts.ES, HALLMARK_G2M_CHECKPOINT.list)
gsva_INFLAMMATORY_RESPONSE <- gsva(counts.ES, HALLMARK_INFLAMMATORY_RESPONSE.list)
gsva_HYPOXIA <- gsva(counts.ES, HALLMARK_HYPOXIA.list)
gsva_raas <- gsva(counts.ES, raas.list)

###------------------------------ Visualisation -----------------------------###
################################################################################

annotation <- as.data.frame(factor(coldata[['subtype']])) %>% 
  {colnames(.) <- 'subtype'; .} %>% 
  {rownames(.) <- rownames(coldata); .} %>% 
  {.$batch <- coldata$batch; .}

# Hallmark heatmap -------------------------------------------------------------
max_val <- max(abs(exprs(gsva_halmark)))
my_breaks <- seq(-max_val, max_val, length.out = 101)


hallmark.heatmap <- pheatmap(exprs(gsva_halmark), 
         cluster_rows = F,
         annotation_col = annotation,
         color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
         breaks = my_breaks,
         fontsize_row = 8,  
         fontsize_col = 8
         )

combined_plot <- as.ggplot(hallmark.heatmap) 
combined_plot <- combined_plot + plot_annotation(
  caption = 'Dataset: GSE188486&GSE217386&GSE69197
  Purity estimated
  Subtypes
  Normalised enrichment score
  '
)

ggsave(paste('hallmark_heatplot', '.png', sep = ''),
       dpi = 1000,
       width = 10, height = 10)

# Heatmaps of selected pathways ------------------------------------------------

pheatmap_plot <- function(expression, title){
  
  max_val <- max(abs(exprs(expression)))
  my_breaks <- seq(-max_val, max_val, length.out = 101)
  
  heatmap <- pheatmap(exprs(expression), 
                           cluster_rows = F,
                           annotation_col = annotation,
                           color = colorRampPalette(rev(
                           brewer.pal(11, "RdYlBu")))(100),
                           breaks = my_breaks,
                           fontsize_row = 8,  
                           fontsize_col = 8,
                           cellheight=10, cellwidth = 10,
                           show_rownames = F,
                           main = title,
                           legend = T,
                           annotation_legend = F
  )
  
}

# RAAS
raas.heatmap <- pheatmap_plot(gsva_raas, 'Renin_Angiotensin_Aldosterone_Pathway')

# HALLMARK_ANGIOGENESIS
ANGIOGENESIS.heatmap <- pheatmap_plot(gsva_ANGIOGENESIS, 'HALLMARK_ANGIOGENESIS')

# HALLMARK_G2M_CHECKPOINT
G2M_CHECKPOINT.heatmap <- pheatmap_plot(gsva_G2M_CHECKPOINT, 'HALLMARK_G2M_CHECKPOINT')

# HALLMARK_INFLAMMATORY_RESPONSE
INFLAMMATORY_RESPONSE.heatmap <- pheatmap_plot(gsva_INFLAMMATORY_RESPONSE, 'HALLMARK_INFLAMMATORY_RESPONSE')

# HALLMARK_HYPOXIA
HYPOXIA.heatmap <- pheatmap_plot(gsva_HYPOXIA, 'HALLMARK_HYPOXIA')

# Extract the legend of the heatmaps above
legend.heatmap <- pheatmap(exprs(gsva_HYPOXIA), 
                            cluster_rows = F,
                            cluster_cols = F,
                            annotation_col = annotation,
                            color = colorRampPalette(c('white'))(100),
                            fontsize_row = 8,  
                            fontsize_col = 8,
                            show_rownames = F,
                            show_colnames = F,
                            cellheight=0, cellwidth = 0,
                            legend = F,
                            annotation_legend = T,
                            annotation_names_col = F
)


combined <- (as.ggplot(HYPOXIA.heatmap) / as.ggplot(G2M_CHECKPOINT.heatmap) / 
  as.ggplot(INFLAMMATORY_RESPONSE.heatmap)) | (as.ggplot(ANGIOGENESIS.heatmap) /
     as.ggplot(raas.heatmap) / as.ggplot(legend.heatmap)) 
combined <- combined + plot_annotation(
  caption = 'Dataset: GSE188486&GSE217386&GSE69197
  Purity estimated
  Subtypes
  Normalised enrichment score
  '
)

ggsave(paste('selected_pathways_heat', '.png', sep = ''),
       dpi = 1000,
       width = 18, height = 9)
