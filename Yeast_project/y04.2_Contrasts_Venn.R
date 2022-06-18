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
FC_threshold <- 1
x3 <- list(
  l_Fa_C = row.names(Fa_C[which(Fa_C$padj <= pval_threshold & Fa_C$log2FoldChange>FC_threshold ), ]),
  l_Aa_C = row.names(Aa_C[which(Aa_C$padj <= pval_threshold & Fa_C$log2FoldChange>FC_threshold ), ]),
  l_Aa_Fa_C = row.names(Aa_Fa_C[which(Aa_Fa_C$padj <= pval_threshold & Fa_C$log2FoldChange>FC_threshold ), ])
)

x3 <- list(
  l_Fa_C = row.names(Fa_C[which(Fa_C$padj <= pval_threshold & Fa_C$log2FoldChange<FC_threshold ), ]),
  l_Aa_C = row.names(Aa_C[which(Aa_C$padj <= pval_threshold & Fa_C$log2FoldChange<FC_threshold ), ]),
  l_Aa_Fa_C = row.names(Aa_Fa_C[which(Aa_Fa_C$padj <= pval_threshold & Fa_C$log2FoldChange<FC_threshold ), ])
)


source("y00_functions.R")
display_venn(x3,
             category.names = c("Furfural", "Acidic acid", "Acidic acid + furfural"),
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



#DON'T TRUST THE CODE BELOW!!!!

#Extract the names of the genes that are upregulated and interesting: 
l_Fa_C = row.names(Fa_C[which(Fa_C$padj <= pval_threshold & Fa_C$log2FoldChange>FC_threshold ), ])
l_Aa_C = row.names(Aa_C[which(Aa_C$padj <= pval_threshold & Fa_C$log2FoldChange>FC_threshold ), ])
l_Aa_Fa_C = row.names(Aa_Fa_C[which(Aa_Fa_C$padj <= pval_threshold & Fa_C$log2FoldChange>FC_threshold ), ])

#Find the 148 DEGs that are unique to furfural
Unique_DEGs_fur <- l_Fa_C[!l_Aa_Fa_C]  #  %in% l_Fa_C& !l_Aa_C %in% l_Fa_C 
length(Unique_DEGs_fur) #only 139????

#Find the 30 DEGs that are shared with furfural and furfural + acedic acid
shared_DEGs <- l_Fa_C[l_Fa_C %in% l_Aa_Fa_C & !l_Aa_C %in% l_Aa_Fa_C ]
length(shared_DEGs)

#test <- annotation[annotation$TranscriptID %in% rownames(res_Ax_NonAx),]
# input <-  list(
#   DEG_table1 = top20_Axenic_YP206, 
#   DEG_table2 = top20_Axenic_YP26, 
#   DEG_table3 = top20_NonAxenic_YP206
# )

input <- x3
res <- Fa_C
sharedDEGs <- function(input,res){

  temp1 <- input[[1]]
  temp2 <- input[[2]]
  temp3 <- input[[3]]
  
  # Make a list of unique gene names from the table
  all_DEGs_unique <- unique(c(temp1, temp2, temp3))
  
  # Get a TRUE/FALSE vector that tels whether a certain DEG
  names <- temp1[all_DEGs_unique %in% temp1]
  temp11 <- c(all_DEGs_unique %in% temp1)
  temp22 <- c(all_DEGs_unique %in% temp2)
  temp33 <- c(all_DEGs_unique %in% temp3)

  DEGs_allocation <- cbind(names, temp11, temp22, temp33)
  
  # Get their padj values and Log2FC values in furfural
  m <- as.data.frame <- Fa_C
  test <- m[m %in% DEGs_allocation] #Doesnt work
  
  return(DEGs_allocation)
}

#DEGs_test <- sharedDEGs(input)

