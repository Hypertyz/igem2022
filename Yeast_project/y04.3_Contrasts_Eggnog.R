# Pathway analysis
# Go to 5.4 :  https://mkempenaar.github.io/gene_expression_analysis/chapter-5.html

library("seqinr")     # for making fasta files
library("DESeq2")     # for handeling the DESeq object
library("readxl")     # for reading eggnog output excel files
library("pheatmap")   # for generating heatmap
library(scales)       # for making a good color gradient for heatmaps
library("reshape2")   # for formatting the tables that is input to heatmap function
library("ggplot2")

#-----------Load-data-----------------------------------------------------------
setwd("C:/Users/hr/OneDrive - Danmarks Tekniske Universitet/DTU2/iGEM 2022/Transcriptomics/Yeast_data/")

source("y00_functions.R")
#path_to_data <- "C:/Users/hr/Documents/AlgaeProject_R_scripts/sample_info/"
metaData <- read_clean_metaData()
countData <- read_clean_countData()


# made in 03_DESeqObject
dds <- readRDS("R_outputs/dds.rds") #DESeq object
res <- readRDS("R_outputs/res.rds") #result table from DESeq object

# contrasts
Fa_C <- readRDS("R_outputs/Fa_C.rds")
Aa_C <- readRDS("R_outputs/Aa_C.rds")
Aa_Fa_C <- readRDS("R_outputs/Aa_Fa_C6.rds")

#------------subset into up/down regulated genes--------------------------------
p_threshold <- 0.05     #to get significant DEGs only

Fa_C_up <- subset(Fa_C, padj < p_threshold & log2FoldChange > 1.1)
nrow(Fa_C_up) # 2254

Fa_C_down <- subset(Fa_C, padj < p_threshold & log2FoldChange < 1.1)
nrow(Fa_C_down) # 2196

# Aa_C_up <- subset(Aa_C, padj < p_threshold & log2FoldChange > 0)
# nrow(Aa_C_up) # 2106
# 
# Aa_C_down <- subset(Aa_C, padj < p_threshold & log2FoldChange < 0)
# nrow(Aa_C_down) # 2047
# 
# Aa_Fa_C_up <- subset(Aa_Fa_C, padj < p_threshold & log2FoldChange > 0)
# nrow(Aa_Fa_C_up) # 2383
# 
# Aa_Fa_C_down <- subset(Aa_Fa_C, padj < p_threshold & log2FoldChange < 0)
# nrow(Aa_Fa_C_down) # 2354

# Next step: Make tables for each where you merge with the annotated eggnog table that you made for your reference genome. 
#Then, you can make your heatmap with COG categories.
# Or maybe Jake can help with some enrichment. 
# But you can also just make a pheatmap

#----Make pheatmap---
library(pheatmap)

Fa_C_up_counts <- countData[countData$feat %in% rownames(Fa_C_up), ]
Fa_C_down_counts <- countData[countData$feat %in% rownames(Fa_C_down), ]

# Aa_C_up_counts <- countData[countData$feat %in% rownames(Aa_C_up), ]
# Aa_C_down_counts <- countData[countData$feat %in% rownames(Aa_C_down), ]
# 
# Aa_Fa_C_up_counts <- countData[countData$feat %in% rownames(Aa_Fa_C_up), ]
# Aa_Fa_C_down_counts <- countData[countData$feat %in% rownames(Aa_Fa_C_down), ]

Fa_C_up_counts1 <- Fa_C_up_counts

for (i in ncol(Fa_C_up_counts1[,-1])){
  temp <- Fa_C_up_counts1[,i]
  temp <- ifelse(temp==0,0,log10(temp))
  Fa_C_up_counts1[,i] <- temp
}

Dist <- data.matrix(Fa_C_up_counts1, rownames.force = NA)
rownames(Dist) = Fa_C_up_counts$feat

pheat_save <- pheatmap(Dist, fontsize = 10)

Dist <- data.matrix(Fa_C_down_counts, rownames.force = NA)
rownames(Dist) = Fa_C_down_counts$feat
pheat_save <- pheatmap(Dist, fontsize = 10)

ggsave(
  "pheat_save_down.pdf",
  plot = pheat_save,
  path = "R_outputs/",
  height = 10,
  width = 10,
  dpi = 1000
)

#---LOOK AT OLD METHOD TO CLASSIFY GENES THAT ARE NOT WELL-ANNOTATED!!!!

