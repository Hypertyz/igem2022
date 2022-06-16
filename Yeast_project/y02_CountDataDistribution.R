# Find the distribution of your data to determine statistical model. 
# DESeq uses a negative binomial model to analyse data where the mean < variance. 
#Tutorial: https://github.com/hbctraining/DGE_workshop/blob/master/lessons/01_DGE_setup_and_overview.md

#-----------Libraries-and-path--------------------------------------------------
library("ggplot2")
library("patchwork")
setwd("C:/Users/hr/OneDrive - Danmarks Tekniske Universitet/DTU2/iGEM 2022/Transcriptomics/Yeast_data/") #DTU2/iGEM/Transcriptomics/Yeast_data

#-----------Load-count-and-metadata---------------------------------------------

source("y00_functions.R")
#path_to_data <- "C:/Users/hr/Documents/AlgaeProject_R_scripts/sample_info/"
metaData <- read_clean_metaData()
countData <- read_clean_countData()

head(countData)


#-----------Investigate count distribution plots--------------------------------

Fa1_CountDistribution <- ggplot(countData) +
                              geom_histogram(aes(x = Fa1), stat = "bin", bins = 200, color="Navy") +
                              #xlim(-5,30000)+
                              xlab("Raw expression counts") +
                              ylab("Number of genes") + theme_Publication()+
                              theme(axis.title.x = element_text(size = 15), 
                                    axis.title.y = element_text(size = 15), 
                                    axis.text.x = element_text(size=12), 
                                    axis.text.y = element_text(size=12))
Fa1_CountDistribution



# Investigate if mean < variance
mean_counts <- apply(countData[,2:13], 1, mean)
variance_counts <- apply(countData[, 2:13], 1, var)
df <- data.frame(mean_counts, variance_counts)

p_variance_mean <- ggplot(df) +
                  geom_point(aes(x=mean_counts, y=variance_counts), alpha=0.3, color= "Navy") + 
                  geom_line(aes(x=mean_counts, y=mean_counts, color="red"), size=1)+
                  scale_y_log10() +
                  scale_x_log10() + 
                  xlab("Mean of counts") +
                  ylab("Variance of counts") + theme_Publication()+
                  theme(axis.title.x = element_text(size = 15), 
                      axis.title.y = element_text(size = 15), 
                      axis.text.x = element_text(size=12), 
                      axis.text.y = element_text(size=12),
                      legend.position="none")


p_variance_mean

save <- Fa1_CountDistribution + p_variance_mean
save <- save + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 17))
save



ggsave(
  "Count_distribution1.pdf",
  plot = Fa1_CountDistribution,
  path = "R_outputs/",
  height = 6,
  width = 5,
  dpi = 1000
)

ggsave(
  "p_variance_mean.pdf",
  plot = p_variance_mean,
  path = "R_outputs/",
  height = 6,
  width = 5,
  dpi = 1000
)

ggsave(
  "Count_distribution.pdf",
  plot = save,
  path = "R_outputs/",
  height = 6,
  width = 10,
  dpi = 1000
)

# As the variance tends to be greater than the mean -> use negative binomial distribution with DESeq.
