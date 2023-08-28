# Author: Aldas Žarnauskas
# Title: Obtaining RAS & RAAS Gene Sets
# Date: 09/007/2023
# Last time updated: 09/007/2023
# Objective: obtain a list of genes associated with RAS & RAAS


### ----------------------------- Libraries --------------------------------- ##
################################################################################
library(msigdbr)
library(biomaRt)

### ----------------------- Obtaining Gene Sets ----------------------------- ##
################################################################################
# Connect to ensembl
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Get a list of gene sets for human species
msigdb <- msigdbr(species = "Homo sapiens")

# Select RAS & RAAS
ras_pathway <- msigdb[msigdb$gs_name == "KEGG_RENIN_ANGIOTENSIN_SYSTEM", ]
raas_pathway <- msigdb[
  msigdb$gs_name == "WP_RENINANGIOTENSINALDOSTERONE_SYSTEM_RAAS", ]



ras_gene_symbols <- ras_pathway[,4] %>% unlist() %>% unname()
raas_gene_symbols <- raas_pathway[,4] %>% unlist() %>% unname()
ras_pathway <- ras_gene_symbols %>% sort()
ras_raas_pathways <- unique(c(ras_gene_symbols, raas_gene_symbols)) %>% sort()



# Get annotations
annotations_ras_pathway <- getBM(
  attributes = c("ensembl_gene_id_version", "external_gene_name", 
                 'description', 'entrezgene_id'),
  filters = "external_gene_name",
  values = ras_pathway,
  mart = ensembl)

# Get annotations
annotations_ras_raas_pathways <- getBM(
  attributes = c("ensembl_gene_id_version", "external_gene_name", 
                 'description', 'entrezgene_id'),
  filters = "external_gene_name",
  values = ras_raas_pathways,
  mart = ensembl)


# write gene sets info to a csv file
write.csv(annotations_ras_raas_pathways, 
          'C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/ras_raas_genes.csv', row.names = F)

write.csv(annotations_ras_pathway, 
          'C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/ras_genes.csv', row.names = F)
