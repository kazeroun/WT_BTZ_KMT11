
################################################################################
###### WT dataset
################################################################################
setwd("/project/otmc/mkazerou/K36/RNA_analysis/NEW/BiometasFiles_from_Mehmet/")
WORKDIR="/project/otmc/mkazerou/K36/RNA_analysis/NEW/BiometasFiles_from_Mehmet/"
#setwd(paste0(WORKDIR,"/quant"))
print(getwd())

library(DESeq2)
library(tximport)
library(EnhancedVolcano)
library(ggplot2)

# Read raw gene counts
counts <- read.table("Biometas_k36_RawCountsgene.txt",
                     header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)

rownames(counts) <- counts$Symbol
counts$Symbol <- NULL

counts <- counts[!grepl("^MT-", rownames(counts)), ]
counts <- counts[rowSums(counts) > 0, ]

keep <- rowSums(counts >= 10) >= 3
counts <- counts[keep, ]

head(rownames(counts))


sampleTable=read.csv("/project/otmc/mkazerou/K36/RNA_analysis/NEW/BiometasFiles_from_Mehmet/sample_info.csv")

sampleTable$sample <- c("S_05uM.1","S_05uM.2","S_05uM.3",
                        "S_15uM.1","S_15uM.2","S_15uM.3",
                        "S_DMSO.1","S_DMSO.2","S_DMSO.3")

sampleTable

colData <- sampleTable
rownames(colData) <- colData$sample

colData$group <- factor(colData$group, levels = c("DMSO","KTX_0.5","KTX_1.5"))

counts <- counts[, rownames(colData)]
stopifnot(all(colnames(counts) == rownames(colData)))

dds <- DESeqDataSetFromMatrix(counts, colData, design = ~ group)
dds <- DESeq(dds)
resultsNames(dds)

# 1.5 µM vs DMSO
res_15_vs_dmso <- results(dds, name = "group_KTX_1.5_vs_DMSO")
summary(res_15_vs_dmso)

res_df_15 <- as.data.frame(res_15_vs_dmso)
res_df_15$gene <- rownames(res_df_15)
res_df_15 <- res_df_15[!is.na(res_df_15$padj), ]

min_padj <- min(res_df_15$padj[res_df_15$padj > 0])
res_df_15$padj_plot <- ifelse(res_df_15$padj == 0,
                              min_padj * 0.1,
                              res_df_15$padj)

library(dplyr)

EnhancedVolcano(
  res_df_15,
  lab = res_df_15$gene,
  #  selectLab = top20,
  x = "log2FoldChange",
  y = "padj_plot",
  pCutoff = 0.05,
  FCcutoff = 1.0,
  pointSize = 4.0,
  labSize = 2.0,
  title = "Suspension: KTX 1.5 µM vs DMSO",
  subtitle = "",
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = "grey",
  max.overlaps = 20
)

# Subset significant differentially expressed genes
res_df_sig <- res_df_15[
  !is.na(res_df_15$padj) &
    abs(res_df_15$log2FoldChange) > 1 &
    res_df_15$padj < 0.05,
]

dim(res_df_sig)

write.csv(
  as.data.frame(res_df_sig),
  file = "WT_KTX_1.5_vs_DMSO_DGE_log2FC1_padj0.05.csv"
)

##### PATHWAY
library(clusterProfiler)
library(org.Hs.eg.db)

padj_cutoff <- 0.05
lfc_cutoff  <- 1

up_genes_15 <- rownames(res_df_15[
  res_df_15$padj < padj_cutoff & res_df_15$log2FoldChange >= lfc_cutoff,
])

down_genes_15 <- rownames(res_df_15[
  res_df_15$padj < padj_cutoff & res_df_15$log2FoldChange <= -lfc_cutoff,
])

length(up_genes_15)
length(down_genes_15)

## Map gene symbols → Entrez IDs
up_entrez_15 <- bitr(
  up_genes_15,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)$ENTREZID

down_entrez_15 <- bitr(
  down_genes_15,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)$ENTREZID

## KEGG pathway enrichment

ekegg_up_15 <- enrichKEGG(
  gene         = up_entrez_15,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

ekegg_up_15 <- setReadable(
  ekegg_up_15,
  OrgDb   = org.Hs.eg.db,
  keyType = "ENTREZID"
)

ekegg_down_15 <- enrichKEGG(
  gene         = down_entrez_15,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

ekegg_down_15 <- setReadable(
  ekegg_down_15,
  OrgDb   = org.Hs.eg.db,
  keyType = "ENTREZID"
)

## KEGG visualisation

dotplot(ekegg_up_15, showCategory = 15) +
  ggtitle("KEGG – UP genes (KTX 1.5 µM vs DMSO)")

dotplot(ekegg_down_15, showCategory = 15) +
  ggtitle("KEGG – DOWN genes (KTX 1.5 µM vs DMSO)")

write.csv(as.data.frame(ekegg_up_15),
          "KEGG_UP_KTX_1.5_vs_DMSO.csv",
          row.names = FALSE)

write.csv(as.data.frame(ekegg_down_15),
          "KEGG_DOWN_KTX_1.5_vs_DMSO.csv",
          row.names = FALSE)


################################################################################
###### BTZ dataset
################################################################################

setwd("/project/otmc/mkazerou/K36/RNA_analysis/NEW/")
WORKDIR="/project/otmc/mkazerou/K36/RNA_analysis/NEW/"
setwd(paste0(WORKDIR,"/quant"))
print(getwd())

library(DESeq2)
library(tximport)
library(EnhancedVolcano)

tx2genesymbol=read.table(paste0(WORKDIR,"/tx2genesymbol_new.txt"),header=T,stringsAsFactors=F)
sampleTable=read.csv("/project/otmc/mkazerou/K36/RNA_analysis/sample_info.csv")

files=list.files(pattern="quant")
names(files)=gsub("_quant.sf","",files)
all(file.exists(files)) # must be TRUE
txi <- tximport(files, type="salmon", tx2gene=tx2genesymbol) #gene level
TXI <- tximport(files, type="salmon", txOut = TRUE) #transcript level
TXI$counts["ERCC-00171",]
dds_tx <- DESeqDataSetFromTximport(txi, colData = sampleTable, design = ~ group)

assays(dds_tx)[["avgTxLength"]]=NULL
is_ercc <- grepl("^ERCC-", rownames(dds_tx))
dds_tx <- estimateSizeFactors(dds_tx, controlGenes = is_ercc)
sizeFactors(dds_tx)

dds <- DESeqDataSetFromTximport(txi, colData = sampleTable, design = ~ group)
#sizeFactors(dds)
sizeFactors(dds) <- sizeFactors(dds_tx)
dim(dds)
# Remove mitochondrial genes
dds_filtered <- dds[ !grepl("^MT-", rownames(dds)), ]
dim(dds_filtered)
# Keep only rows that have a count of 10 or more in at least 3 samples
keep <- rowSums(counts(dds_filtered) >= 10) >= 3
dds_final <- dds_filtered[keep,]
dim(dds_final)

# Diff Exp

# Convert columns to factors
dds_final$sample_type <- as.factor(dds_final$sample_type)
dds_final$group <- as.factor(dds_final$group)

print(levels(dds_final$sample_type))
print(levels(dds_final$group))

design(dds_final) <- ~ sample_type + group + sample_type:group

# Rerun the main DESeq analysis
dds_result_interaction <- DESeq(dds_final)

resultsNames(dds_result_interaction)

# KTX 1.5 vs. DMSO in Suspension cells
res_sus_KTX_1.5_vs_DMSO <- results(dds_result_interaction,
                                   contrast = list(c("group_KTX_1.5_vs_DMSO", "sample_typeSuspension.groupKTX_1.5")))
write.csv(as.data.frame(res_sus_KTX_1.5_vs_DMSO),
          file = "Suspension_KTX_1.5_vs_DMSO_DGE.csv")

library(EnhancedVolcano)

EnhancedVolcano(res_sus_KTX_1.5_vs_DMSO,
                lab = rownames(res_sus_KTX_1.5_vs_DMSO),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 2.0,
                title = "Suspension: KTX 1.5 vs. DMSO",
                subtitle = 'Differential Expression Analysis',
                caption = paste0("Fold Change Cutoff: 1.0; Adjusted P-value Cutoff: 0.05; Total genes: ", nrow(res_sus_KTX_1.5_vs_DMSO)),
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'grey',
                max.overlaps = 50)

# Subset significant differentially expressed genes
res_sus_KTX_sig <- res_sus_KTX_1.5_vs_DMSO[
  !is.na(res_sus_KTX_1.5_vs_DMSO$padj) &
    abs(res_sus_KTX_1.5_vs_DMSO$log2FoldChange) > 1 &
    res_sus_KTX_1.5_vs_DMSO$padj < 0.05,
]

dim(res_sus_KTX_sig)

write.csv(
  as.data.frame(res_sus_KTX_sig),
  file = "Suspension_KTX_1.5_vs_DMSO_DGE_log2FC1_padj0.05.csv"
)


#----------------------------------------------------
#PATHWAY ANALYSIS
#----------------------------------------------------
library(clusterProfiler)  # for ID mapping
library(org.Hs.eg.db)     # human mapping database

padj_cutoff <- 0.05
log2fc_cutoff <- 1

res1.5 <- res_sus_KTX_1.5_vs_DMSO[!is.na(res_sus_KTX_1.5_vs_DMSO$padj), ]
up1.5   <- rownames(res1.5[res1.5$padj < padj_cutoff & res1.5$log2FoldChange > log2fc_cutoff, ])
down1.5 <- rownames(res1.5[res1.5$padj < padj_cutoff & res1.5$log2FoldChange < -log2fc_cutoff, ])

up1.5_entrez   <- bitr(up1.5, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
down1.5_entrez <- bitr(down1.5, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# KEGG enrichment
ekegg_up_1.5 <- enrichKEGG(up1.5_entrez$ENTREZID, organism='hsa', pAdjustMethod="BH", pvalueCutoff=0.05)
ekegg_down_1.5 <- enrichKEGG(down1.5_entrez$ENTREZID, organism='hsa', pAdjustMethod="BH", pvalueCutoff=0.05)

# Plot (with safety checks)
plot_if_results <- function(enrich_obj, title) {
  if (nrow(as.data.frame(enrich_obj)) > 0) {
    print(dotplot(enrich_obj, showCategory=10, title=title))
  } else {
    message(paste("⚠️ No significant enrichment for:", title))
  }
}

# KEGG
plot_if_results(ekegg_up_1.5, "KEGG: Upregulated (Suspension KTX 1.5 vs DMSO)")
plot_if_results(ekegg_down_1.5, "KEGG: Downregulated (Suspension KTX 1.5 vs DMSO)")


write.csv(as.data.frame(ekegg_up_1.5), "KEGG_Up_Suspension_KTX_1.5_vs_DMSO.csv")
write.csv(as.data.frame(ekegg_down_1.5), "KEGG_Down_Suspension_KTX_1.5_vs_DMSO.csv")


################################################################################
###### Common genes in datasets
################################################################################

BTZ=res_sus_KTX_sig
WT=res_df_sig

genes_up_WT <- rownames(WT[WT$log2FoldChange > 0, ])
genes_down_WT <- rownames(WT[WT$log2FoldChange < 0, ])

genes_up_BTZ <- rownames(BTZ[BTZ$log2FoldChange > 0, ])
genes_down_BTZ <- rownames(BTZ[BTZ$log2FoldChange < 0, ])

common_up_15   <- intersect(genes_up_WT, genes_up_BTZ)
common_down_15 <- intersect(genes_down_WT, genes_down_BTZ)

length(common_up_15)
length(common_down_15)

write.csv(common_up_15,
          "Common_UP_genes_KTX_1.5_Biometas_vs_BTZ.csv",
          row.names = FALSE)

write.csv(common_down_15,
          "Common_DOWN_genes_KTX_1.5_Biometas_vs_BTZ.csv",
          row.names = FALSE)


################################################################################
###### ER Stress heatmap
################################################################################

genes_er <- c(
  "ERN1","XBP1","ATF6","ATF4","DDIT3","PPP1R15A","EIF2AK3",
  "PRDX4","ERO1A","P4HB","PDIA5","PDIA6","HMOX1","SLC7A11",
  "TXNIP","GPX7","GPX8",
  "DERL3","SEL1L","SYVN1","HERPUD1","HERPUD2",
  "DERL1","DERL2","OS9","ERLEC1","UBE2J1",
  "AUP1","UBE2G2","FAF2","HSPA5","HSP90B1"
)

genes_er[genes_er %in% rownames(dds_final)]

setwd("/project/otmc/mkazerou/K36/RNA_analysis/NEW/")
WORKDIR="/project/otmc/mkazerou/K36/RNA_analysis/NEW/"
setwd(paste0(WORKDIR,"/quant"))
print(getwd())

library(DESeq2)
library(tximport)
library(EnhancedVolcano)

tx2genesymbol=read.table(paste0(WORKDIR,"/tx2genesymbol_new.txt"),header=T,stringsAsFactors=F)
sampleTable=read.csv("/project/otmc/mkazerou/K36/RNA_analysis/sample_info.csv")

files=list.files(pattern="quant")
names(files)=gsub("_quant.sf","",files)
all(file.exists(files)) # must be TRUE
txi <- tximport(files, type="salmon", tx2gene=tx2genesymbol) #gene level
TXI <- tximport(files, type="salmon", txOut = TRUE) #transcript level
TXI$counts["ERCC-00171",]
dds_tx <- DESeqDataSetFromTximport(txi, colData = sampleTable, design = ~ group)

assays(dds_tx)[["avgTxLength"]]=NULL
is_ercc <- grepl("^ERCC-", rownames(dds_tx))
dds_tx <- estimateSizeFactors(dds_tx, controlGenes = is_ercc)
sizeFactors(dds_tx)

dds <- DESeqDataSetFromTximport(txi, colData = sampleTable, design = ~ group)
#sizeFactors(dds)
sizeFactors(dds) <- sizeFactors(dds_tx)
dim(dds)
# Remove mitochondrial genes
dds_filtered <- dds[ !grepl("^MT-", rownames(dds)), ]
dim(dds_filtered)
# Keep only rows that have a count of 10 or more in at least 3 samples
keep <- rowSums(counts(dds_filtered) >= 10) >= 3
dds_final <- dds_filtered[keep,]
dim(dds_final)

# Diff Exp
dds_final$sample_type <- as.factor(dds_final$sample_type)
dds_final$group <- as.factor(dds_final$group)

print(levels(dds_final$sample_type))
print(levels(dds_final$group))

design(dds_final) <- ~ sample_type + group + sample_type:group

# Rerun the main DESeq analysis
dds_result_interaction <- DESeq(dds_final)

resultsNames(dds_result_interaction)

#----------------------------------------------------

# KTX 0.5 vs. DMSO
res_KTX_0.5_vs_DMSO <- results(dds_result_interaction, 
                               contrast = c("group", "KTX_0.5", "DMSO"))

# KTX 1.5 vs. DMSO
res_KTX_1.5_vs_DMSO <- results(dds_result_interaction, 
                               contrast = c("group", "KTX_1.5", "DMSO"))

# KTX vs. DMSO in Adherent cells
res_adh_KTX_0.5_vs_DMSO <- results(dds_result_interaction, 
                                   name = "group_KTX_0.5_vs_DMSO")

res_adh_KTX_1.5_vs_DMSO <- results(dds_result_interaction, 
                                   name = "group_KTX_1.5_vs_DMSO")

# KTX 0.5 vs. DMSO in Suspension cells
res_sus_KTX_0.5_vs_DMSO <- results(dds_result_interaction,
                                   contrast = list(c("group_KTX_0.5_vs_DMSO", "sample_typeSuspension.groupKTX_0.5")))

# KTX 1.5 vs. DMSO in Suspension cells
res_sus_KTX_1.5_vs_DMSO <- results(dds_result_interaction,
                                   contrast = list(c("group_KTX_1.5_vs_DMSO", "sample_typeSuspension.groupKTX_1.5")))

# KTX 1.5 vs. KTX 0.5 in Adherent cells
res_adh_KTX_1.5_vs_0.5 <- results(dds_result_interaction,
                                  contrast = list("group_KTX_1.5_vs_DMSO", "group_KTX_0.5_vs_DMSO"))

res_sus_15 <- as.data.frame(
  res_sus_KTX_1.5_vs_DMSO[rownames(res_sus_KTX_1.5_vs_DMSO) %in% genes_er, ]
)

res_adh_15 <- as.data.frame(
  res_KTX_1.5_vs_DMSO[rownames(res_KTX_1.5_vs_DMSO) %in% genes_er, ]
)

res_sus_15$condition <- "Suspension_1.5"
res_adh_15$condition <- "Adherent_1.5"

res_sus_15$gene <- rownames(res_sus_15)
res_adh_15$gene <- rownames(res_adh_15)

res_sus_05 <- as.data.frame(
  res_sus_KTX_0.5_vs_DMSO[rownames(res_sus_KTX_0.5_vs_DMSO) %in% genes_er, ]
)

res_adh_05 <- as.data.frame(
  res_adh_KTX_0.5_vs_DMSO[rownames(res_adh_KTX_0.5_vs_DMSO) %in% genes_er, ]
)

res_sus_05$condition <- "Suspension_0.5"
res_adh_05$condition <- "Adherent_0.5"

res_sus_05$gene <- rownames(res_sus_05)
res_adh_05$gene <- rownames(res_adh_05)


res_all <- rbind(
  res_sus_15,
  res_adh_15,
  res_sus_05,
  res_adh_05
)

colnames(res_all)

library(tidyr)

heat_df <- res_all %>%
  dplyr::select(gene, condition, log2FoldChange) %>%
  tidyr::pivot_wider(
    names_from = condition,
    values_from = log2FoldChange
  )
head(heat_df)
head(heat_df$gene)

heat_df <- as.data.frame(heat_df)   # convert tibble → data.frame
rownames(heat_df) <- heat_df$gene
heat_df$gene <- NULL

heat_df <- as.data.frame(heat_df)
rownames(heat_df) <- heat_df$gene
heat_df$gene <- NULL

heat_mat <- as.matrix(heat_df)
heat_mat[!is.finite(heat_mat)] <- 0

heat_mat <- heat_mat[, c(
  "Suspension_0.5", "Suspension_1.5",
  "Adherent_0.5", "Adherent_1.5"
)]


library(pheatmap)

heat_mat <- apply(heat_mat, 2, as.numeric)

pheatmap(
  heat_mat,
#  scale = "row",
  show_rownames = TRUE,
  fontsize_row = 8,
  cellheight = 10,
  color = colorRampPalette(c("blue","white","red"))(100),
  main = "ER Stress / UPR Genes\nKTX vs DMSO (BTZ)"
)
