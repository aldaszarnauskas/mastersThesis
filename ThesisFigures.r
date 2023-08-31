# Author: Aldas Žarnauskas
# Title: Thesis Figures
# Date: 28/008/2023
# Last time updated: 28/008/2023
# Objective: 


# Clear environment
rm(list=ls(all.names = TRUE))

#---------------------------------- Libraries ---------------------------------#
################################################################################

library(tibble)
library(dplyr)
library(ggplot2)
library(ggstance)
library(forcats)
library(ComplexHeatmap)
library(biomaRt)
library(RColorBrewer)


#--------------------------- Data Import --------------------------------------#
################################################################################

#Tissue data
normalised_counts <- read.csv("C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/Norm_PCA_BOX_Heat/normalised_countssubtypesGSE188486&GSE217386&GSE69197purityEstimated.csv")
coldata <- read.csv("C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/colData_kidney_cancer_tumourPurity.csv")

# Replace NA values in counts with 0, if necessary
any(is.na(normalised_counts))
counts[is.na(normalised_counts)] <- 0
any(is.na(normalised_counts))

# #connect to ensembl BioMart database
# ensembl <- useEnsembl("ensembl")
# 
# #select a dataset, here homo sapiens, and update mart object using useDataset
# ensembl <- useDataset('hsapiens_gene_ensembl', mart = ensembl)


# Loading RAAS Genes -----------------------------------------------------------
genes_of_interest <- read.csv('C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/ras_raas_genes.csv')
genes_of_interest.df <- genes_of_interest[, c(1, 2)] 
# ------------------------------------------------------------------------------
df <- normalised_counts %>% 
  column_to_rownames(., var = "X") %>% 
  as.data.frame() %>% 
  t(.) %>% scale(.) %>% t(.) %>% as.data.frame() %>% 
  rownames_to_column(., var = 'ensembl_gene_id_version') %>% 
  merge(genes_of_interest[, c(1, 2)], by = 'ensembl_gene_id_version') %>% 
  column_to_rownames(., 'external_gene_name') %>% .[, 2:ncol(.)]


all_DE_res <- readRDS("C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/DE_plots/all_DE_results.rds")

threedatasetsPEchRCC_vs_normal <- all_DE_res[["threedatasetsPEchRCC_vs_normal"]]
threedatasetsPEccRCC_vs_normal <- all_DE_res[["threedatasetsPEccRCC_vs_normal"]]
threedatasetsPEpRCC_vs_normal <- all_DE_res[["threedatasetsPEpRCC_vs_normal"]]
threedatasetsPEccRCC_vs_chRCC <- all_DE_res[["threedatasetsPEccRCC_vs_chRCC"]]
threedatasetsPEpRCC_vs_chRCC <- all_DE_res[["threedatasetsPEpRCC_vs_chRCC"]]
threedatasetsPEpRCC_vs_ccRCC <- all_DE_res[["threedatasetsPEpRCC_vs_ccRCC"]]




threedatasetsPEchRCC_vs_normal.df <- threedatasetsPEchRCC_vs_normal %>% 
  as.data.frame() %>% 
  {rownames_to_column(., var = 'ensembl_gene_id_version')}
threedatasetsPEchRCC_vs_normal_RAS <- filter(
  threedatasetsPEchRCC_vs_normal.df, ensembl_gene_id_version %in% 
    genes_of_interest.df$ensembl_gene_id_version) %>% 
  base::merge(genes_of_interest.df, by = 'ensembl_gene_id_version') %>% 
  {column_to_rownames(., var = "external_gene_name")} %>% 
  .[3] %>% as.matrix()


threedatasetsPEccRCC_vs_normal.df <- threedatasetsPEccRCC_vs_normal %>% 
  as.data.frame() %>% 
  {rownames_to_column(., var = 'ensembl_gene_id_version')}
threedatasetsPEccRCC_vs_normal_RAS <- filter(
  threedatasetsPEccRCC_vs_normal.df, ensembl_gene_id_version %in% 
    genes_of_interest.df$ensembl_gene_id_version) %>% 
  base::merge(genes_of_interest.df, by = 'ensembl_gene_id_version') %>% 
  {column_to_rownames(., var = "external_gene_name")} %>% 
  .[3] %>% as.matrix()

threedatasetsPEpRCC_vs_normal.df <- threedatasetsPEpRCC_vs_normal %>% 
  as.data.frame() %>% 
  {rownames_to_column(., var = 'ensembl_gene_id_version')}
threedatasetsPEpRCC_vs_normal_RAS <- filter(
  threedatasetsPEpRCC_vs_normal.df, ensembl_gene_id_version %in% 
    genes_of_interest.df$ensembl_gene_id_version) %>% 
  base::merge(genes_of_interest.df, by = 'ensembl_gene_id_version') %>% 
  {column_to_rownames(., var = "external_gene_name")} %>% 
  .[3] %>% as.matrix()

threedatasetsPEccRCC_vs_chRCC.df <- threedatasetsPEccRCC_vs_chRCC %>% 
  as.data.frame() %>% 
  {rownames_to_column(., var = 'ensembl_gene_id_version')}
threedatasetsPEccRCC_vs_chRCC_RAS <- filter(
  threedatasetsPEccRCC_vs_chRCC.df, ensembl_gene_id_version %in% 
    genes_of_interest.df$ensembl_gene_id_version) %>% 
  base::merge(genes_of_interest.df, by = 'ensembl_gene_id_version') %>% 
  {column_to_rownames(., var = "external_gene_name")} %>% 
  .[3] %>% as.matrix()

threedatasetsPEpRCC_vs_chRCC.df <- threedatasetsPEpRCC_vs_chRCC %>% 
  as.data.frame() %>% 
  {rownames_to_column(., var = 'ensembl_gene_id_version')}
threedatasetsPEpRCC_vs_chRCC_RAS <- filter(
  threedatasetsPEpRCC_vs_chRCC.df, ensembl_gene_id_version %in% 
    genes_of_interest.df$ensembl_gene_id_version) %>% 
  base::merge(genes_of_interest.df, by = 'ensembl_gene_id_version') %>% 
  {column_to_rownames(., var = "external_gene_name")} %>% 
  .[3] %>% as.matrix()

threedatasetsPEpRCC_vs_ccRCC.df <- threedatasetsPEpRCC_vs_ccRCC %>% 
  as.data.frame() %>% 
  {rownames_to_column(., var = 'ensembl_gene_id_version')}
threedatasetsPEpRCC_vs_ccRCC_RAS <- filter(
  threedatasetsPEpRCC_vs_ccRCC.df, ensembl_gene_id_version %in% 
    genes_of_interest.df$ensembl_gene_id_version) %>% 
  base::merge(genes_of_interest.df, by = 'ensembl_gene_id_version') %>% 
  {column_to_rownames(., var = "external_gene_name")} %>% 
  .[3] %>% as.matrix()







#----------------------- ACE2 Expression in Cancers ---------------------------#
################################################################################

ace2expr <- data.frame(
  "Small intestine" = c(237.9, "High"),
  "Duodenum" = c(179.2, "High"),
  "Kidney" = c(85.2, "High"),
  "Gallbladder" = c(64.1, "High"),
  "Testis" = c(48.5, "High"),
  "Heart muscle" = c(41.3, "Not detected"),
  "Colon" = c(9.6, "High"),
  "Adipose tissue" = c(8.6, "Not detected"),
  "Rectum" = c(5.6, "High"),
  "Liver" = c(5.6, "Not detected"),
  "Choroid plexus" = c(5.4, "Not detected"),
  "Thyroid gland" = c(4.3, "Not detected"),
  "Breast" = c(4.2, "Not detected"),
  "Seminal vesicle" = c(4.0, "Medium"),
  "Stomach" = c(2.6, "Not detected"),
  "Pancreas" = c(2.4, "Medium"),
  "Ovary" = c(1.8, "Not detected"),
  "Esophagus" = c(1.7, "Not detected"),
  "Parathyroid gland" = c(1.5, "Not detected"),
  "Fallopian tube" = c(1.5, "Medium"),
  "Salivary gland" = c(1.4, "Not detected"),
  "Appendix" = c(1.3, "Low"),
  "Placenta" = c(1.2, "High"),
  "Vagina" = c(1.1, "Not detected"),
  "Urinary bladder" = c(0.8, "Not detected"),
  "Lung" = c(0.7, "Not detected"),
  "Cervix" = c(0.6, "Not detected"),
  "Skeletal muscle" = c(0.6, "Not detected"),
  "Bone marrow" = c(0.6, "Not detected"),
  "Adrenal gland" = c(0.5, "Not detected"),
  "Tongue" = c(0.5, "Not detected"),
  "Epididymis" = c(0.5, "Medium"),
  "Endometrium" = c(0.5, "Not detected"),
  "Smooth muscle" = c(0.5, "Not detected"),
  "Thymus" = c(0.5, "Not detected"),
  "Prostate" = c(0.4, "Not detected"),
  "Lymph node" = c(0.4, "Not detected"),
  "Retina" = c(0.3, "Not detected"),
  "Skin" = c(0.3, "Not detected"),
  "Midbrain" = c(0.2, "Not detected"),
  "Tonsil" = c(0.2, "Not detected"),
  "Cerebral cortex" = c(0.1, "Not detected"),
  "Basal ganglia" = c(0.1, "Not detected"),
  "Hypothalamus" = c(0.1, "Not detected"),
  "Amygdala" = c(0.1, "Not detected"),
  "Hippocampal formation" = c(0.1, "Not detected"),
  "Spinal cord" = c(0.1, "Not detected"),
  "Pituitary gland" = c(0.1, "Not detected"),
  "Spleen" = c(0.1, "Not detected"),
  "Cerebellum" = c(0.0, "Not detected")) %>% t() %>% 
  {colnames(.) <- c("RNA_expression_nTPM", "Protein_expression_score"); .} %>% 
  as.data.frame() %>% 
  {rownames_to_column(., var = 'Organ')} %>%
  mutate(RNA_expression_nTPM = as.numeric(RNA_expression_nTPM)) %>%
  arrange(desc(RNA_expression_nTPM))
  


ggplot(ace2expr, aes(
  RNA_expression_nTPM, fct_reorder(Organ, RNA_expression_nTPM),
  fill = Protein_expression_score)) + 
  geom_barh(stat = 'identity', width = 0.6) +
  scale_x_log10()+
  ylab(NULL) + 
  
  theme_classic() +
  xlab("RNA Expression (log10 nTPM)") +
  #geom_text(xaes(label = Organ), hjust = -0.1, size = 3, color = "black") +
  # geom_text(aes(x = ifelse(RNA_expression_nTPM > 1, 0.4, 2.6), 
  #               label = Organ), size = 3,
  #           data = ace2expr) +
  # 
  guides(fill = guide_legend(title = "Protein expression")) +
  scale_x_continuous(
    trans = "log10",
    breaks = c(0.1, 1, 300)) +
  scale_fill_brewer(palette = "Set1", aesthetics = c("colour", "fill")) +
  theme(legend.position = c(0.65,0.078),
        axis.ticks.y = element_blank(),
        #axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  labs(caption = paste('
                         '))




ggsave(paste('RNA_Protein_expr', '.png', sep = ''),
       dpi = 1000,
       width = 4.5, height = 10)



#---------------------------- Complex Heatmap ---------------------------------#
################################################################################
common_genes <- intersect(
  intersect(
    intersect(
      intersect(
        intersect(
          intersect(
            rownames(df), 
            rownames(threedatasetsPEchRCC_vs_normal_RAS)
          ),
          rownames(threedatasetsPEccRCC_vs_normal_RAS)
        ),
        rownames(threedatasetsPEpRCC_vs_normal_RAS)
      ),
      rownames(threedatasetsPEccRCC_vs_chRCC_RAS)
    ),
    rownames(threedatasetsPEpRCC_vs_chRCC_RAS)
  ),
  rownames(threedatasetsPEpRCC_vs_ccRCC_RAS)
)

color_scale <- colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100)

ht1 = Heatmap(
  df[common_genes,] %>% as.matrix(), 
  col = color_scale,
  name = "z-score")

ht2 = Heatmap(threedatasetsPEchRCC_vs_normal_RAS[common_genes,],
              col = color_scale,
              name = "chVSnormall2FC")

ht3 = Heatmap(threedatasetsPEccRCC_vs_normal_RAS[common_genes,],
              col = color_scale,
              name = "ccVSnormall2FC")

ht4 = Heatmap(threedatasetsPEpRCC_vs_normal_RAS[common_genes,],
              col = color_scale,
              name = "pVSnormall2FC")

ht5 = Heatmap(threedatasetsPEccRCC_vs_chRCC_RAS[common_genes,],
              col = color_scale,
              name = "ccVSchl2FC")

ht6 = Heatmap(threedatasetsPEpRCC_vs_chRCC_RAS[common_genes,],
              col = color_scale,
              name = "pVSchl2FC")

ht7 = Heatmap(threedatasetsPEpRCC_vs_ccRCC_RAS[common_genes,],
              col = color_scale,
              name = "pVSccl2FC")


ht1 + ht2 + ht3 + ht4 + ht5 + ht6 + ht7



#-------------------------------- Methylation RCC -----------------------------#
################################################################################

meth <- read.csv("C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/Probe-level-methylation-data-2023-08-29 (1).csv")
CpG <- read.csv("C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/Aggregation-methylation-data-2023-08-29.csv")

meth_RCC <- filter(
  meth, 
  tissue == "KIRC" |
  tissue == "KIRP" |
  tissue == "KICH"
  ) %>% 
  #{.$tissue <- ifelse(.$type == "Normal", "normal", .$tissue); .}
  {.$tissue <- ifelse(
    .$type == "Normal" & .$tissue == 'KIRC', "cc_normal", 
    ifelse(.$type == "Normal" & .$tissue == 'KIRP', "p_normal", 
           ifelse(.$type == "Normal" & .$tissue == 'KICH', "ch_normal",
                  .$tissue))); .}
  



CpG_RCC <- filter(
  CpG, 
  tissue == "KIRC" |
  tissue == "KIRP" |
  tissue == "KICH"
) %>% 
  #{.$tissue <- ifelse(.$type == "Normal", "normal", .$tissue); .}
  {.$tissue <- ifelse(
    .$type == "Normal" & .$tissue == 'KIRC', "cc_normal", 
    ifelse(.$type == "Normal" & .$tissue == 'KIRP', "p_normal", 
           ifelse(.$type == "Normal" & .$tissue == 'KICH', "ch_normal",
                  .$tissue))); .}

ggplot(CpG_RCC, 
       aes(x = eval(as.symbol("tissue")), 
           y = eval(as.symbol("value")), fill = eval(as.symbol("tissue")))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = eval(as.symbol("tissue"))), width = 0.3, 
              alpha = 0.5) +
  labs(x = "tissue", y = "ACE2 Beta value", caption = "CpG_RCC") +
  scale_fill_brewer(palette = "Set1", aesthetics = c("colour", "fill"), 
                    guide = "none") +
  theme(
    axis.title = element_text(
      family = "Helvetica", size = (10), colour = "black"),
    axis.text = element_text(
      family = "Courier", colour = "black", size = (10)),
    plot.caption= element_text(size=7, face = 'italic',
                               color="Black"),
    legend.position = 'none') +
  #guides(color = "none", size = "none", fill = "none") +
  theme_classic() #+
  ggsignif::geom_signif(
    comparisons = comparison.list,
    annotations = annotations,
    y_position = c(max_val, (max_val+(max_val*.3)), (max_val+(max_val*.6)),
                   (max_val+(max_val*.9)), (max_val+(max_val*1.2)),
                   (max_val+(max_val*1.5)))) 
  


ggsave(paste('CpG_methylation', sep = '', '.png'),
       dpi = 1000,
       width = 5, height = 5)



ggplot(CpG_RCC, aes(x = tissue, y = value, fill = type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = tissue), width = 0.3, alpha = 0.5) +
  labs(x = "Tissue", y = "ACE2 Beta value", caption = "CpG_RCC") +
  #scale_fill_brewer(palette = "Set1", aesthetics = c("colour", "fill"), guide = "none") +
  theme(
    axis.title = element_text(family = "Helvetica", size = 10, colour = "black"),
    axis.text = element_text(family = "Courier", colour = "black", size = 10),
    plot.caption = element_text(size = 7, face = 'italic', color = "black"),
    legend.position = 'right'
  ) 
