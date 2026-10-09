library(readr)



tts <- read_delim("ttsdown_counts_May.txt", 
                  delim = "\t", escape_double = FALSE, 
                  comment = "#", trim_ws = TRUE)

tts_gtf <- read_delim("tts_downstream.gtf", 
                      delim = "\t", escape_double = FALSE, 
                      col_names = FALSE, trim_ws = TRUE)
tts_gtf <- tts_gtf[tts_gtf$X3=="gene",]

library(GenomicRanges)
gr <- GRanges(
  seqnames = tts_gtf$X1,
  ranges = IRanges(start = tts_gtf$X4, end = tts_gtf$X5)
)

gr.red <- GenomicRanges::reduce(gr)


tts_size <- sum(width(gr.red))


tts_totals <- (colSums(tts[,7:42]) / tts_size) * 1e06


tss <- read_delim("tssup_counts_May.txt", 
                  delim = "\t", escape_double = FALSE, 
                  comment = "#", trim_ws = TRUE)

tss_gtf <- read_delim("tss_upstream.gtf", 
                      delim = "\t", escape_double = FALSE, 
                      col_names = FALSE, trim_ws = TRUE)
tss_gtf <- tss_gtf[tss_gtf$X3=="gene",]

library(GenomicRanges)
gr <- GRanges(
  seqnames = tss_gtf$X1,
  ranges = IRanges(start = tss_gtf$X4, end = tss_gtf$X5)
)

gr.red <- GenomicRanges::reduce(gr)


tss_size <- sum(width(gr.red))

tss_totals <- (colSums(tss[,7:42]) / tss_size) * 1e06


exons <- read_delim("exon_counts_May.txt", 
                    delim = "\t", escape_double = FALSE, 
                    comment = "#", trim_ws = TRUE)

exons_gtf <- read_delim("exons.gtf", 
                        delim = "\t", escape_double = FALSE, 
                        col_names = FALSE, trim_ws = TRUE)
library(GenomicRanges)
gr <- GRanges(
  seqnames = exons_gtf$X1,
  ranges = IRanges(start = exons_gtf$X4, end = exons_gtf$X5)
)

gr.red <- GenomicRanges::reduce(gr)


exons_size <- sum(width(gr.red))

exons_totals <- (colSums(exons[,7:42]) / exons_size) * 1e06


introns <- read_delim("introns_counts_May.txt", 
                      delim = "\t", escape_double = FALSE, 
                      comment = "#", trim_ws = TRUE)

introns_gtf <- read_delim("introns.gtf", 
                          delim = "\t", escape_double = FALSE, 
                          col_names = FALSE, trim_ws = TRUE)
introns_gtf <- introns_gtf[introns_gtf$X3=="gene",]

library(GenomicRanges)
gr <- GRanges(
  seqnames = introns_gtf$X1,
  ranges = IRanges(start = introns_gtf$X4, end = introns_gtf$X5)
)

gr.red <- GenomicRanges::reduce(gr)


intron_size <- sum(width(gr.red))



introns_totals <- (colSums(introns[,7:42]) / intron_size) * 1e06


igr <- read_delim("intergenic_counts_May.txt", 
                  delim = "\t", escape_double = FALSE, 
                  comment = "#", trim_ws = TRUE)

igr_gtf <- read_delim("gencode_intergenic.gtf", 
                      delim = "\t", escape_double = FALSE, 
                      col_names = FALSE, trim_ws = TRUE)
igr_gtf <- igr_gtf[igr_gtf$X3=="gene",]
library(GenomicRanges)
gr <- GRanges(
  seqnames = igr_gtf$X1,
  ranges = IRanges(start = igr_gtf$X4, end = igr_gtf$X5)
)

gr.red <- GenomicRanges::reduce(gr)


igr_size <- sum(width(gr.red))

igr_totals <- (colSums(igr[,7:42]) / igr_size) * 1e06


df <- data.frame(Exons=exons_totals,Intron=introns_totals,TSS=tss_totals,TTS=tts_totals,IGR=igr_totals)

dff <- apply(df, 1, function(x){x/sum(x)})
cnames <- colnames(dff)
cnames <- gsub(pattern = "./results/bwa/merged_library/",replacement = "",x = cnames)
cnames <- gsub(pattern = ".mLb.clN.sorted.bam",replacement = "",x = cnames)

colnames(dff) <- cnames
dff <- dff[,!grepl(pattern = "Input",x = colnames(dff))]
dff <- dff[,!grepl(pattern = "-A-",x = colnames(dff))]

library(pheatmap)
pheatmap(dff)

pdf(file = "usage.pdf",width = 7,height = 3.5)
pheatmap(dff)
dev.off()







library(tidyverse)


dff_long <- dff %>%
  as.data.frame() %>%
  rownames_to_column(var = "Region") %>%
  pivot_longer(
    cols = -Region, 
    names_to = "Sample", 
    values_to = "Value"
  ) %>%
  
  mutate(
    Condition = case_when(
      str_detect(Sample, "-T-") ~ "T",
      str_detect(Sample, "-C-") ~ "C",
      TRUE ~ NA_character_
    ),
    PointColor = case_when(
      str_detect(Sample, "-B-") ~ "B",
      str_detect(Sample, "-S-") ~ "S",
      TRUE ~ NA_character_
    )
  )


ggplot(dff_long, aes(x = Condition, y = Value)) +
  geom_boxplot(outlier.shape = NA, fill = "gray95", color = "gray40") +
  geom_jitter(aes(color = PointColor), width = 0.15, size = 3, alpha = 0.8) +
  facet_wrap(~ Region, scales = "free_y", nrow = 1) +
  theme_bw(base_size = 14) +
  labs(
    x = "Condition",
    y = "Value",
    color = "Sample Type"
  ) +
  scale_color_manual(values = c("B" = "#377eb8", "S" = "#e41a1c")) 



library(tidyverse)
library(ggpubr) 


dff_long <- dff %>%
  as.data.frame() %>%
  rownames_to_column(var = "Region") %>%
  pivot_longer(
    cols = -Region, 
    names_to = "Sample", 
    values_to = "Value"
  ) %>%
  mutate(
    Condition = case_when(
      str_detect(Sample, "-T-") ~ "T",
      str_detect(Sample, "-C-") ~ "C",
      TRUE ~ NA_character_
    ),
    PointColor = case_when(
      str_detect(Sample, "-B-") ~ "B",
      str_detect(Sample, "-S-") ~ "S",
      TRUE ~ NA_character_
    )
  )


p1 <- ggplot(dff_long, aes(x = Condition, y = Value)) +
  geom_boxplot(outlier.shape = NA, fill = "gray95", color = "gray40") +
  geom_jitter(aes(color = PointColor), width = 0.15, size = 3, alpha = 0.8) +
  
  
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.format",      
    label.x = 1.5,           
    vjust = 0             
  ) +
  
  facet_wrap(~ Region, nrow = 1) +
  
  # Style and aesthetics
  theme_bw(base_size = 14) +
  labs(
    x = "Condition",
    y = "Normalized signal intensity",
    color = "Sample Type"
  ) +
  scale_color_manual(values = c("B" = "#377eb8", "S" = "#e41a1c"))


pdf(file = "BvsS_comparison_normalized_usage.pdf",width = 9,height = 4)
print(p1)
dev. Off()
