# mastersThesis

mastersThesis repository contains scripts used by Aldas Žarnauskas to conduct his master's project. The project used raw RNA-seq data followed by various downstream analyses, including differential gene expression analysis (DGEA), gene set enrichment analysis (GSEA), single sample gene set enrichment analysis (ssGSEA) and transcription factor enrichment analysis (TFEA).

Objectives of scripts:
CountGenerationBASH.txt - take raw RNA-seq counts and generate a count table

MSigDB_RCC.r - obtain RNA-seq metadata

mSigDB_ras_raas.r - obtain renin-angiotensin system genes information

tumourpurityESTIMATE_RCC.r - estimate tumour purity in the RNA-seq counts.

Preproessing.r - preprocess RNA-seq count tables (e.g. batch removal, normalisation, initial visualisation)

DEAnalysis_RCC.R - perform DGEA 

DE_visualisation_RCC.r - visualise DGEA results

DoRoTHea_RCC.r - perform TFEA

GSEA_clusterProfiler_RCC.r - perform GSEA

ssGSEA_RCC.r - perform ssGSEA



