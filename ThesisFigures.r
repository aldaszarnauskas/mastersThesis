library(tibble)
library(dplyr)
library(ggplot2)
library(ggstance)
library(forcats)

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
