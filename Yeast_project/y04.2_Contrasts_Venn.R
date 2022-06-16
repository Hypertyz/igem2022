# Explore contrasts

# A contrast is a linear combination of estimated log2 fold changes, 
# which can be used to test if differences between groups are equal to zero. 
# The simplest use case for contrasts is an experimental design containing a 
# factor with three levels, say A, B and C. Contrasts enable the user to generate 
# results for all 3 possible differences: log2 fold change of B vs A, of C vs A, and of C vs B. 
# The contrast argument of results function is used to extract test results of log2 fold changes of interest, for example:

#-----------Libraries-and-path--------------------------------------------------
library("ggplot2") #for plotting
library("DESeq2") # for DESeq analysis

#-----------Load-data-----------------------------------------------------------
setwd("C:/Users/hr/OneDrive - Danmarks Tekniske Universitet/DTU2/iGEM 2022/Transcriptomics/Yeast_data/")

# made in 03_DESeqObject
dds <- readRDS("R_outputs/dds.rds") #DESeq object
res <- readRDS("R_outputs/res.rds") #result table from DESeq object

# contrasts
Fa_C <- readRDS("R_outputs/Fa_C.rds")
Aa_C <- readRDS("R_outputs/Aa_C.rds")
Aa_Fa_C <- readRDS("R_outputs/Aa_Fa_C6.rds")

#-------------------------------limited Venn diagram----------------------------
# Only investigating what significant genes they have in common
pval_threshold <- 0.05
x3 <- list(
  l_Fa_C = row.names(Fa_C[which(Fa_C$padj <= pval_threshold), ]),
  l_Aa_C = row.names(Aa_C[which(Aa_C$padj <= pval_threshold), ]),
  l_Aa_Fa_C = row.names(Aa_Fa_C[which(Aa_Fa_C$padj <= pval_threshold), ])
)

source("y00_functions.R")
display_venn(x3,
             category.names = c("Fa", "Aa", "Aa_Fa"),
             # Circles
             lwd = 2,
             lty = 'blank',
             fill =c("#E69F00", "#0072B2","#009E73"),
             # Numbers
             fontfamily = "sans",
             cex = 1.6,
             # Set names
             cat.cex = 1.3,
             cat.fontface = "bold",
             cat.fontfamily = "sans",
             cat.default.pos = "outer",
             cat.dist = c(-0.03, -0.03, -0.03)
)

# Save the picture in high quality
venn.diagram(x3, filename = "R_outputs/VennD.png",
             category.names = c("YP206 vs. Axenic", "YP26 vs. Axenic", "NonAxenic vs. Axenic"),
             # Circles
             lwd = 2,
             lty = 'blank',
             fill =c("#E69F00", "#0072B2","#009E73"),
             # Numbers
             fontfamily = "sans",
             cex = 1.6,
             # Set names
             cat.cex = 1.4,
             cat.fontface = "bold",
             cat.fontfamily = "sans",
             cat.default.pos = "outer",
             cat.dist = c(-0.03, -0.03, -0.03),
             imagetype = "png"
)



#FIX THIS CODE! OLD STUFF SAVE FOR LATER!
#------Write tables with DEGs in different contrasts----------------------------
#-----WE DON'T HAVE ANY ANNOTATION YET!-----------------------------------------

# from 00_functions-----------CHANGE THIS!--------------------------------------
# source("y00_functions.R")
# top20_Aa_Fa_C <- makeTop20DEGsTable(res_Ax_NonAx, annotation)
# top20_Axenic_YP206 <- makeTop20DEGsTable(Fa_C, annotation)
# top20_Axenic_YP26 <- makeTop20DEGsTable(Aa_C, annotation)
# 
# #-----Function that write the names of genes in common
# #Maybe change this so you get ALL of the genes so that you can for example get the three common genes.
# gc(verbose = TRUE, reset = TRUE)
# # prepare data
# input <-  list(
#   top20_Aa_Fa_C = top20_Aa_Fa_C,
#   top20_Axenic_YP206     = top20_Axenic_YP206, 
#   top20_Axenic_YP26      = top20_Axenic_YP26)
# 
# # from 00_functions
# DEGs_allocation <- sharedDEGs(input)
# 
# #merge to annotation WRITE THIS AS A FUNCTION!!!!!
# 
# DEG_names <- DEGs_allocation$all_DEGs_unique
# DEG_annotation <- annotation[annotation$TranscriptID == DEG_names[1],]
# 
# for (i in DEG_names[-1]){
#   temp = annotation[annotation$TranscriptID == i,]
#   DEG_annotation=rbind(DEG_annotation,temp)
# }
# 
# #choosing the first information in InterPro section
# interpro_info <- DEG_annotation$InterPro
# interpro_info <- sapply(strsplit(interpro_info, ";"), '[',1)
# 
# DEG_annotation_simple <- cbind(DEGs_allocation,
#                                DEG_annotation$Name, 
#                                DEG_annotation$Product, 
#                                interpro_info)
# 
# write.csv2(DEG_annotation_simple,"R_outputs/top20_Shared_DEGs.csv", row.names = T)
# 
# View(DEG_annotation_simple)

# -----trash for now---
#maybe this is much smarter but idk. instead of make tables functions
# test <- annotation[annotation$TranscriptID %in% rownames(res_Ax_NonAx),]
# View(test)
