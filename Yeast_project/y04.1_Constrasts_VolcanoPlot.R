#Explore the distribution of the mRNA reads with volcano plots

#-----------Libraries-and-path--------------------------------------------------
library("ggplot2") #for plotting
library("DESeq2") # for DESeq analysis
library("EnhancedVolcano")
library("patchwork") # for merging the three volcano plots

#-----------Load-data-----------------------------------------------------------
setwd("C:/Users/hr/OneDrive - Danmarks Tekniske Universitet/DTU2/iGEM 2022/Transcriptomics/Yeast_data/")

# made in 03_DESeqObject
dds <- readRDS("R_outputs/dds.rds") #DESeq object
res <- readRDS("R_outputs/res.rds") #result table from DESeq object

# contrasts
Fa_C <- readRDS("R_outputs/Fa_C.rds")
Aa_C <- readRDS("R_outputs/Aa_C.rds")
Aa_Fa_C <- readRDS("R_outputs/Aa_Fa_C6.rds")


#-----volcano-plots-with-contrast------------------

# Defining values for basic VolcanoPLot 
padj1 = .05
padj2 = .01
minlog2FoldChange = 1

library("EnhancedVolcano")
library(ggpubr)
theme_set(theme_pubr())

pCutoff = 0.05
FCcutoff = 1.0

#Interesting genes

# YOR388C
# YNR034W-A
# YOL052C-A
chosen_genes <- c("YOR388C","YNR034W-A","YOL052C-A")

VolcanoPlot1 <- EnhancedVolcano(Fa_C,
                                        lab = rownames(Fa_C),
                                        selectLab = chosen_genes,
                                        #selectLab = NA,
                                        x = 'log2FoldChange',
                                        y = 'pvalue',
                                        title = "Furfural",
                                        subtitle = bquote(italic(" ")),
                                        legendPosition = 'top',
                                        pCutoff = pCutoff,
                                        FCcutoff = FCcutoff,
                                        drawConnectors = TRUE,
                                        widthConnectors = 1.0,
                                        ylim = c(0,350),
                                        xlim = c(-5,8)
                                        )
                                        

VolcanoPlot2 <- EnhancedVolcano(Aa_C,
                                        lab = rownames(Aa_C),
                                        selectLab = chosen_genes,
                                        #selectLab = NA,
                                        x = 'log2FoldChange',
                                        y = 'pvalue',
                                        title = "Acetic acid",
                                        subtitle = bquote(italic(" ")),
                                        legendPosition = 'top',
                                        pCutoff = pCutoff,
                                        FCcutoff = FCcutoff,
                                        drawConnectors = TRUE,
                                        widthConnectors = 1.0,
                                        ylim = c(0,350),
                                        xlim = c(-5,8)
                                        )

VolcanoPlot3 <- EnhancedVolcano(Aa_Fa_C,
                                        lab = rownames(Aa_Fa_C),
                                        selectLab = chosen_genes,
                                        #selectLab = NA,
                                        x = 'log2FoldChange',
                                        y = 'pvalue',
                                        title = "Acitid acid and furfural",
                                        subtitle = bquote(italic(" ")),
                                        legendPosition = 'top',
                                        pCutoff = pCutoff,
                                        FCcutoff = FCcutoff,
                                        drawConnectors = TRUE,
                                        widthConnectors = 1.0,
                                        ylim = c(0,350),
                                        xlim = c(-5,8)
                                        )


VolcanoPlot1 + VolcanoPlot2 + VolcanoPlot3
Combined <- VolcanoPlot1 + VolcanoPlot2 + VolcanoPlot3

ggsave(
  "Volcano_Combined.pdf",
  plot = Combined,
  path = "R_outputs/",
  height = 7,
  width = 12,
  dpi = 1000
)


#-----------------------OLD STUFF-SAVE FOR LATER--------------------------------
# genes chosen from looking at the volcano plots
#FUN_003385-T1
#FUN_009941-T1
#FUN_006116-T1
#FUN_001827-T1
#FUN_003811-T1
#FUN_002313-T1

# par(mfrow=c(2,3))
# plotCounts(dds, gene="FUN_003385-T1", intgroup="sample_type")
# plotCounts(dds, gene="FUN_009941-T1", intgroup="sample_type")
# plotCounts(dds, gene="FUN_006116-T1", intgroup="sample_type")
# plotCounts(dds, gene="FUN_001827-T1", intgroup="sample_type")
# plotCounts(dds, gene="FUN_003811-T1", intgroup="sample_type")
# plotCounts(dds, gene="FUN_002313-T1", intgroup="sample_type")
# dev.off()

# chosen_genes <- c("FUN_003385-T1","FUN_009941-T1", "FUN_006116-T1", "FUN_001827-T1", "FUN_003811-T1", "FUN_002313-T1")
# chosen_genes_tab <- annotation[annotation$TranscriptID %in% chosen_genes,]
# View(chosen_genes_tab)


