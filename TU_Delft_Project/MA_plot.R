#  Testing the effect of the lfcShrink

#-----------Libraries-and-path--------------------------------------------------
library("ggplot2") #for plotting
library("DESeq2") # for DESeq analysis
library("EnhancedVolcano")
library("patchwork") # for merging the three volcano plots

#-----------Load-data-----------------------------------------------------------
#setwd("C:/Users/hr/OneDrive - Danmarks Tekniske Universitet/DTU2/iGEM 2022/Transcriptomics/Yeast_data/")

# made in 03_DESeqObject
dds <- readRDS("R_outputs/dds.rds") #DESeq object
res <- readRDS("R_outputs/res.rds") #result table from DESeq object

# contrasts
GHB_NH4_unshrunken <- readRDS("R_outputs/GHB_NH4_unshrunken.rds")
Suc_GABA_unshrunken <- readRDS("R_outputs/Suc_GABA_unshrunken.rds")

GHB_NH4 <- readRDS("R_outputs/GHB_NH4.rds")
Suc_GABA <- readRDS("R_outputs/Suc_GABA.rds")

GHB_NH4
sum(is.na(GHB_NH4$padj))  #2630


pdf(file = "./R_outputs/MA_plots.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 6) 

par(mfrow = c(2, 2))
plotMA(GHB_NH4_unshrunken,
       colNonSig = "black",
       colSig = "blue",
       colLine = "grey40", 
       main = "GHB_NH4_unshrunken")
plotMA(GHB_NH4,
       colNonSig = "black",
       colSig = "blue",
       colLine = "grey40", 
       main = "GHB_NH4_shrunken")

plotMA(Suc_GABA_unshrunken,
       colNonSig = "black",
       colSig = "blue",
       colLine = "grey40",
       main = "Suc_GABA_unshrunken")
plotMA(Suc_GABA,
       colNonSig = "black",
       colSig = "blue",
       colLine = "grey40",
       main = "Suc_GABA_shrunken")
dev.off()


