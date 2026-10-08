# R code for RNA-seq analysis

This repository contains the R code used for RNA-seq differential expression analysis and figure generation for the associated manuscript.

Two independent RNA-seq datasets (**WT** and **BTZ**) were analysed to assess transcriptional responses to **KTX treatment (0.5 and 1.5 µM)**. The BTZ dataset contains both **adherent and suspension** cells, whereas the WT dataset contains suspension cells only. Comparisons between the WT and BTZ datasets were therefore restricted to **suspension cells**.

The R script includes:

- Differential expression analysis using **DESeq2**
- Volcano plots of KTX treatment effects
- KEGG pathway enrichment analysis
- Identification of differentially expressed genes shared between the WT and BTZ datasets
- Analysis and heatmap visualisation of ER-stress/UPR-associated genes across KTX doses and cell states

The main R packages used include `DESeq2`, `tximport`, `EnhancedVolcano`, `clusterProfiler`, `org.Hs.eg.db`, `ggplot2`, `dplyr`, `tidyr`, and `pheatmap`.
