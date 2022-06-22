# Find the distribution of your data to determine statistical model. 
# DESeq uses a negative binomial model to analyse data where the mean < variance. 
#Tutorial: https://github.com/hbctraining/DGE_workshop/blob/master/lessons/01_DGE_setup_and_overview.md

#-----------Libraries-and-path--------------------------------------------------
library("ggplot2") #for plotting
library("DESeq2") # for DESeq analysis
library("vegan") #to make adonis analysis
#-----------Load-data-----------------------------------------------------------

setwd("C:/Users/hr/OneDrive - Danmarks Tekniske Universitet/DTU2/iGEM 2022/Transcriptomics/Yeast_data/")

source("y00_functions.R")
#path_to_data <- "C:/Users/hr/Documents/AlgaeProject_R_scripts/sample_info/"
metaData <- read_clean_metaData()
countData <- read_clean_countData()
head(countData)
head(metaData)
#-----------Make-DESeq-Object---------------------------------------------------
# CHANGE CHANGE CHANGE!!!!! STEAL JAKES CODE

# Make sure you have right format to make the DESeq Object:
nrow(metaData) # X
ncol(countData) # X+1

dds_dat <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~condition, tidy = TRUE)
dds <- DESeq(dds_dat)
res <- results(dds)

res <- res[order(res$padj),]
summary(res)
head(res)

metaData
#save dds and res because they are time consuming to generate in each script
#for some reason, I have to reset the wd here
#setwd("C:/Users/hr/Documents/AlgaeProject_R_scripts")    

saveRDS(dds, "R_outputs/dds.rds")
saveRDS(res, "R_outputs/res.rds")


#------make-contrasts-----------------------------------------------------------

#p_threshold <- 0.05

Fa_C <- results(dds, contrast=c("condition","Fa","C"))
#res_Ax_NonAx <- subset(res_Ax_NonAx,res_Ax_NonAx$padj < p_threshold) #142

Aa_C <- results(dds, contrast=c("condition","Aa","C"))
#res_Ax_YP206 <- subset(res_Ax_YP206,res_Ax_YP206$padj < p_threshold) #390

Aa_Fa_C <- results(dds, contrast=c("condition","Aa_Fa","C"))
#res_Ax_YP26 <- subset(res_Ax_YP26,res_Ax_YP26$padj < p_threshold) #56

# sort so that the small padj are in the top of the res tables.
Fa_C <- Fa_C[order(Fa_C$padj),]
Aa_C <- Aa_C[order(Aa_C$padj),]
Aa_Fa_C <- Aa_Fa_C[order(Aa_Fa_C$padj),]

# count all the DEGs in the three contrasts
nrow(Fa_C) 
nrow(Aa_C) 
nrow(Aa_Fa_C) 
# 6716 for all

#Check they all look OK and doesnt contain NAs
head(Fa_C)
head(Aa_C)
head(Aa_Fa_C)

summary(Fa_C)
summary(Aa_C)
summary(Aa_Fa_C)

saveRDS(Fa_C, "R_outputs/Fa_C.rds")
saveRDS(Aa_C, "R_outputs/Aa_C.rds")
saveRDS(Aa_Fa_C, "R_outputs/Aa_Fa_C.rds")




