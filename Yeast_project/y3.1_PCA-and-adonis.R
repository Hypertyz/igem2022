#Explore raw data with PCA plots and test if the groupings are significant with adonis

#-----------Libraries-and-path--------------------------------------------------
library("ggplot2") #for plotting
library("DESeq2") # for DESeq analysis
library("vegan") #to make adonis analysis

#-----------Load-data-----------------------------------------------------------

setwd("C:/Users/hr/OneDrive - Danmarks Tekniske Universitet/DTU2/iGEM 2022/Transcriptomics/Yeast_data/")

source("y00_functions.R")
# made in 03_DESeqObject
dds <- readRDS("R_outputs/dds.rds") #DESeq object
res <- readRDS("R_outputs/res.rds") #result table from DESeq object

metaData <- read_clean_metaData() # metadata about samples
countData <- read_clean_countData() # countdata about samples

# -------------Explore data-PCA-plot--------------------------------------------

# First we need to transform the raw count data
# vst function will perform variance stabilizing transformation

vsdata <- vst(dds, blind=FALSE)

# nice colour palette for plotting (maybe delete here, because I made it in the function)
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

PCA_top500 <- plotPCA(vsdata, intgroup="condition")

source("y00_functions.R") #to get theme_Publication()

PCA_top500_save <- PCA_top500 + theme_Publication() + 
  theme(strip.text.x = element_blank(), 
        strip.text.y = element_blank(), 
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16, hjust=1), axis.text.y = element_text(size = 16, hjust = 1),
        legend.key.size = unit(0.8, "cm"), 
        legend.text=element_text(size=15),
        legend.position = c(0.87, 0.87),legend.background = element_rect(fill = "white", color = "black")) + 
  geom_text(mapping = aes(label = name), size = 4.5, vjust=1.5) + 
  scale_colour_Publication()
PCA_top500_save

ggsave(
  "PCA_top500_save.pdf",
  plot = PCA_top500_save,
  path = "R_outputs/",
  height = 7,
  width = 12,
  dpi = 1000
)

#--------------------------------ADONIS-----------------------------------------

# Make an ANOVA test/ ADONIS test
# Tutorial: https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/adonis

# format countData for adonis
countData_adonis <-t(countData[,-1])
colnames(countData_adonis) <- countData[, 1]

ncol(countData_adonis) #11126
set.seed(123456789)
adonis_sample_type <- adonis(countData_adonis ~ condition, data = metaData, permutations = 10000) # p = 2e-04 ***
# The p value is HIGHLY significant.

