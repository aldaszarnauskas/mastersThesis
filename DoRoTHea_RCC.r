# Author: Aldas Žarnauskas
# Title: Transcription Factor Analysis (TFA)
# Date: 26/006/2023
# Last time updated: 22/008/2023
# Objective: perform TFEA & visualise the results

#Clear environment
rm(list=ls(all.names = TRUE))

###-------------------------------- Library ---------------------------------###
################################################################################
library(biomaRt)
library(dplyr)
library(dorothea)
library(bcellViper)
library(viper)
library(ggplot2)
library(tibble)
library(pheatmap)
library(ggsignif)
library(dplyr)
library(rstatix)
library(car)


###----------------------- Read in and format data --------------------------###
################################################################################

# Tissue data
counts <- read.csv("C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/Norm_PCA_BOX_Heat/normalised_countssubtypesGSE188486&GSE217386&GSE69197purityEstimated.csv")
coldata <- read.csv("C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/kidney_cancer/colData_kidney_cancer_tumourPurity.csv")

# Connect to ensembl BioMart database
ensembl <- useEnsembl("ensembl")

# Select a dataset, here homo sapiens, and update mart object using useDataset
ensembl <- useDataset('hsapiens_gene_ensembl', mart = ensembl)


# Loading RAAS Genes -----------------------------------------------------------
gene_set.df <- read.csv('C:/Users/aldas/OneDrive/Bioinformatics/master_thesis/scripts_data/ras_raas_genes.csv')

# ------------------------------------------------------------------------------
# In this variable, store the transcript IDs of the genes of interest
ensembl_IDs<- counts$X

# Get annotations
annotations <- getBM(attributes = c("ensembl_gene_id_version", "external_gene_name"),
                     filters = "ensembl_gene_id_version",
                     values = ensembl_IDs, 
                     mart = ensembl)


# Join with counts table
dorothea_counts_final <- counts %>%
  mutate(ensembl_gene_id_version = X) %>%
  full_join(., annotations, by = 'ensembl_gene_id_version') %>%
  as.data.frame() %>%
  na.omit() %>% 
  .[.$external_gene_name != '',] %>% 
  .[!duplicated(.$external_gene_name), ]  %>%
  {rownames(.) <- .$external_gene_name; .} %>% 
  select(-c(X, ensembl_gene_id_version, external_gene_name))



###------------------------------ DoRoTHea ----------------------------------###
################################################################################

# Accessing (human) dorothea regulons
data(dorothea_hs, package = "dorothea")

# STAT3 is not a TF of ACE2 in DoRoTHea, although, the article proves differently
ace2_tfs <- data.frame(tf = c('STAT3', "HNF1α", "HNF1β", "FOXA2"), confidence = "NA", 
                    target = 'ACE2', mor = 1)
dorothea_hs <- rbind(dorothea_hs, ace2_tfs)

# Running viper with dorothea regulons
regulonsAB <- dorothea_hs %>%
  filter(confidence %in% c("A", "B") | target=="ACE2")
regulonsAtoE <- dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C", "D", "E", "NA"))


# Extracting TF activities
tf_activities <- run_viper(dorothea_counts_final, regulonsAB, tidy = FALSE, 
                           options =  list(method = "scale", minsize = 4, 
                                           eset.filter = FALSE, cores = 1, 
                                           verbose = FALSE, nes = TRUE))


###------------------------------ Statistics --------------------------------###
################################################################################


dorothea_stat_test <- function(tf_activities, TFname, comparison){
  
  
  # Filter to retain the TF of interest and comparison column
  tf_activities <- tf_activities %>% 
    select(TFname, subtype)
  
  # Run an anova
  res.aov <- aov(as.formula(paste0(colnames(tf_activities)[1], "~ subtype")), data = tf_activities)
  
  # Get residuals
  aov.residuals <- residuals(object = res.aov)
  
  # If only two groups are present, test for normality
  if (length(unique(tf_activities$subtype)) == 2){
    
    # Do shapiro test (tests for normality)
    stest <- shapiro.test(x = aov.residuals)
    
    # Extract pvalue
    p_norm_single <- stest$p.value
    
    # If passes the shapiro test (p > 0.05, normally distributed) carry out t.test
    if (p_norm_single > 0.05){
      
      # T-test
      t_test <- t.test(as.formula(paste0(colnames(tf_activities)[1], "~ subtype")), data = tf_activities)
      
      # Save p-value in the pvalue variable
      pvalue <- t_test$p.value
      
      # Save type of test in a variable
      typeoftest <- "t-test (parametric, 2 samples)"
      
      # Multiple comp
      multiplecomp <- "none"
      
    }
    
    # If not normally distributed (p>0.05), do a mann whitney u test (wilcox test)
    else {
      
      # W-test
      m_test <- wilcox.test(as.formula(paste0(colnames(tf_activities)[1], "~ subtype")), data = tf_activities)
      
      # Save pvalue in 
      pvalue <- m_test$p.value
      
      # Save type of test
      typeoftest <- "Mann Whitney U-test (non-parametric, 2 samples)"
      
      # Multiple comp
      multiplecomp <- "none"
    }
    
  }
  
  # If only one group/subtype has been selected, return error message asking user to select another group
  else if (length(unique(tf_activities$subtype)) < 2){
    
    # Extract pvalue
    pvalue <- NA
    
    typeoftest <- "No statistical test could be carried out, please select more groups"
    
    multiplecomp <- "none"
  }
  
  # If multiple groups have been selected, need to first test for normality
  else {
    
    # Do shapiro test (tests for normality)
    stest <- shapiro.test(x = aov.residuals)
    
    # Extract pvalue
    p_norm_multiple <- stest$p.value
    
    # Levene test (for homogeneity of variance)
    ltest <- leveneTest(as.formula(paste0(colnames(tf_activities)[1], "~ subtype")), data = tf_activities)
    
    # Associated  pvalue
    p_levene_multiple <-ltest$`Pr(>F)`[1]
    
    # If passes the shapiro test (p > 0.05) and levenes test, carry out anova and tukey
    if (p_norm_multiple > 0.05 & p_levene_multiple >0.05){
      
      # Extract anova p value
      pvalue <- summary(res.aov)[[1]][[5]][[1]]
      
      # If p <0.05, do a post hoc test, if not, return none for multiplecomp
      if(pvalue < 0.05){
        
        # Perform Tukey hsd
        tukey <- TukeyHSD(res.aov)
        
        # Extract comparison names and pvalues for table
        tuk_comp <- rownames(tukey$subtype)
        tuk_pval <- tukey$subtype[,4]
        
        # Join together for multiplecomp
        multiplecomp <- c(tuk_pval)
        
        typeoftest <- "ANOVA with TukeyHSD post hoc"
        
      }
      
      # If p > 0.05, don't carry out a multiple test
      else{
        typeoftest <- "ANOVA only"
        multiplecomp <- "none"
      }
      
    }
    
    # If only passes shapiro and not levene, carry out welch's anova followed by games howell posthoc
    else if (p_norm_multiple > 0.05 & p_levene_multiple < 0.05){ 
      
      # Perform Welch's ANOVA
      w_anova <- oneway.test(as.formula(paste0(colnames(tf_activities)[1], "~ subtype")), data  = tf_activities, var.equal = FALSE)
      
      # Extract pvalue
      pvalue <- w_anova$p.value
      
      # If p<0.05, carry out a multiple test, if not stop here
      
      if(pvalue <0.05){
        
        # Carry out games howell multiple comp
        g_h_comp <- games_howell_test(as.formula(paste0(colnames(tf_activities)[1], "~ subtype")), data  = tf_activities)
        
        # Extract comparisons and pvalues
        
        comp_d <- paste(g_h_comp$group1, "-", g_h_comp$group2, sep = '')
        pvals_d <- g_h_comp$p.adj
        names(pvals_d) <- comp_d
        
        multiplecomp <- c(pvals_d)
        
        typeoftest <- "Welch's ANOVA with post hoc Games-Howell Test"
        
      }
      
      else{
        typeoftest<- "Welch's ANOVA only"
        multiplecomp <- "none"
      }
      
      
      
    }
    
    # If doesnt pass shapiro, carry out KW and dunns test
    else {
      
      # KW test
      kw_test <- kruskal.test(as.formula(paste0(colnames(tf_activities)[1], "~ subtype")), data = tf_activities)
      
      # Extract pvalue
      pvalue <- kw_test$p.value
      
      # If p<0.05, carry out a Dunn's test, if not, stop here
      if(pvalue<0.05){
        
        # Dunns test
        d_test <- dunn_test(as.formula(paste0(colnames(tf_activities)[1], "~ subtype")), data = tf_activities, p.adjust.method = "BH")
        
        # Multiple comparisons
        comp_d <- paste(d_test$group1, "-", d_test$group2, sep = '')
        pvals_d <- d_test$p.adj
        names(pvals_d) <- comp_d
        
        multiplecomp <- c(pvals_d)
        
        typeoftest <- "Kruskal Wallis test with post hoc Dunn's Test (non-parametric)"
        
      }
      
      else{
        typeoftest <- "Kruskal Wallis test only"
        multiplecomp <- 'none'
      }
      
    }
    
    
    
  }
  
  return(list(typeoftest, pvalue, multiplecomp))  
  
}

# Store transcription factors associated with RAS pathway genes
ras_tf <- c()
for (gene in (gene_set.df$external_gene_name %>% unlist())){
  tfs <- regulonsAB[regulonsAB$target == gene,]$tf %>% unlist()
  ras_tf <- c(ras_tf, tfs)
}
ras_tf <- ras_tf %>% unique() %>% sort() %>% intersect(rownames(tf_activities))

# Prepare tf_activities in a suitable format
tf_activities_RCC <- tf_activities[ras_tf,] %>% # Select only TFs associated with RAS
  t(.) %>% as.data.frame() %>%
  rownames_to_column(., var = 'run_accession') %>%  
  merge(coldata[, c('run_accession', 'subtype')], by = 'run_accession')


# Perform multiple pairwise comparison for each RAAS TF
tf_res <- list()
for (tf in ras_tf){
  if (tf %in% colnames(tf_activities_RCC)){
    res <- dorothea_stat_test(tf_activities_RCC, tf)
    tf_res[[tf]] <- res
  }
}

###---------------------------- Visualisation -------------------------------###
################################################################################

tf.activities.plot <- function(tf_activities, TF, comparison,
                               caption = ''){
  
  p.val = tf_res[[TF]][[3]] %>% unname()
  
  comparison.list <- list()
  comp <- names(tf_res[['NR3C1']][[3]])
  for (c in comp){
    c <- strsplit(c, '-')
    comparison.list <- append(comparison.list, c)
  }
  
  
  stars <- function(p) {
    case_when(
      p <= 0.001 ~ "***",
      p <= 0.01  ~ "**",
      p <= 0.05  ~ "*",
      TRUE       ~ "ns"
    )
  }
  
  max_val = max(tf_activities[[TF]])
  
  annotations <- stars(p.val)
  
  ggplot(tf_activities, 
         aes(x = eval(as.symbol(comparison)), 
             y = eval(as.symbol(TF)), fill = eval(as.symbol(comparison)))) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(aes(color = eval(as.symbol(comparison))), width = 0.3, 
                alpha = 0.5) +
    labs(x = comparison, y = "Normalised Enrichment Score", caption = caption) +
    scale_fill_brewer(palette = "Set1", aesthetics = c("colour", "fill"), 
                      guide = "none") +
    ggsignif::geom_signif(
      comparisons = comparison.list,
      annotations = annotations,
      y_position = c(max_val, (max_val+(max_val*.3)), (max_val+(max_val*.6)),
                     (max_val+(max_val*.9)), (max_val+(max_val*1.2)),
                     (max_val+(max_val*1.5)))) +
    theme(
      axis.title = element_text(
        family = "Helvetica", size = (10), colour = "black"),
      axis.text = element_text(
        family = "Courier", colour = "black", size = (10)),
      plot.caption= element_text(size=7, face = 'italic',
                                 color="Black"),
      legend.position = 'none') +
    guides(color = FALSE, size = FALSE, fill = F) +
    theme_classic()
}


tf.regulons <- function(transcription_factor){
  
  tfs <- regulonsAtoE[regulonsAtoE$tf==transcription_factor,] %>%
    filter(target %in% gene_set.df$external_gene_name) %>% 
    mutate(TF.confidence = paste(
      target, '(', confidence, ', ', mor, ')', sep = '')) %>% 
    .$TF.confidence %>% paste(collapse = ', ')
  
  return(tfs)
  
}



tf_activities_RCC[["subtype"]] <- factor(
  tf_activities_RCC[["subtype"]], 
  levels = c("normal", 
             sort(unique(tf_activities_RCC[["subtype"]])) %>% 
               .[-which(. == "normal")]
  ))
  


for (tf in names(tf_res)){
  if (tf_res[[tf]][2] <= .05){
    comparison <- "subtype"
    tf.plot <- tf.activities.plot(tf_activities_RCC, comparison = 'subtype',
                                  TF = tf,
                       caption = paste('Dataset: GSE188486,
                                   ',tf, ' regulates: ', tf.regulons(tf), "
                                       TF Activity ", "of ", tf, ' in ', comparison,
                                       sep = ''))
    ggsave(paste(tf, '_boxplot', sep = '', '.png'),
           dpi = 1000,
           width = 5, height = 5)
    
  }
}  

