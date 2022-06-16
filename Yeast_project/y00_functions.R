#Functions

###
### read and clean data
###

read_clean_metaData <- function(){
  #path_to_data <- "C:/Users/hr/Documents/AlgaeProject_R_scripts/sample_info/"
  #setwd(path_to_data)
  metaData <- read.csv("sample_info/metaData.csv", 
                       header = TRUE, sep=";")
  return(metaData)
}

read_clean_countData <- function(...){
  #path_to_data <- "C:/Users/hr/Documents/AlgaeProject_R_scripts/sample_info/"
  #setwd(path_to_data)
  
  #countData formated in MergeData.R
  metaData <- read.csv("sample_info/metaData.csv", 
                       header = TRUE, sep=";")
  countData <- read.csv("sample_info/MergeTables_collector_yeast.csv", 
                        header=T, sep=",")
  colnames(countData)[-1] <-metaData$Name
  return(countData)
}


###
### For Venn Diagrams------CHECK-------------------------
###

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}



#OLD STUFF THAT NEEDS TO BE REDONE TO FIT OUR DATA

###
### Make fasta file from DEGs from res table------CHECK-------------------------
###

# This function is used to get a fasta file with the genes from the DEGs that are inputted to the function

#setwd("C:/Users/hr/Documents/AlgaeProject_R_scripts")
DEGs_fasta <- function(resDEGs, DEGs_fasta){
  index = which(names(DEGs_fasta) %in% 
                  rownames(resDEGs))
  niels=DEGs_fasta[index]
  
  path = "C:/Users/hr/OneDrive - Danmarks Tekniske Universitet/DTU2/iGEM 2022/Transcriptomics/Yeast_data/R_outputs/"
  path_name = paste0(path, "orfs_", deparse(substitute(resDEGs)), ".fasta")
  
  seqinr::write.fasta(sequences = niels, 
                      names = names(niels), 
                      file.out = path_name)
}

###
### makeTop20DEGsTable-----------------------------------------------------------
###

# Make top20DEGs table
# result_table is any result table made from dds
# the function will sort the table so it is ordered by adjusted p value
# the DEGs with lowest padj will be in top. 
# the top 20 DEGs will be saved in a table called top20DEGs_name-of-table-inserted-into-function

# this is actually really not a very nice function.

#Maybe this should not be top 20 but something you could just choose to do afterwards?
makeTop20DEGsTable = function(result_table, annotation_table){
  result_table_sort <- result_table[order(result_table$padj),]
  #Extract some top genes
  top20DEG  <- head(result_table_sort, 20)
  # list of top20 gene names
  top20DEG_names <- rownames(top20DEG)
  
  # Make a table with annotation of the 20 most significant genes
  top20DEG_annotation <- annotation_table[annotation_table$TranscriptID == top20DEG_names[1],]
  
  for (i in top20DEG_names[-1]){
    temp = annotation_table[annotation_table$TranscriptID == i,]
    top20DEG_annotation=rbind(top20DEG_annotation,temp)
  }
  
  #choosing the first information in InterPro section
  interpro_info <- top20DEG_annotation$InterPro
  interpro_info <- sapply(strsplit(interpro_info, ";"), '[',1)
  
  top20DEG_annotation_simple <- cbind(top20DEG_annotation$TranscriptID,
                                      top20DEG_annotation$Name, 
                                      top20DEG_annotation$Product, 
                                      interpro_info)
  
  colnames(top20DEG_annotation_simple) <- c("TranscriptID","Name", "Product", "InterPro")
  
  path = "C:/Users/hr/Documents/AlgaeProject_R_scripts/R_outputs/"
  path_name = paste0(path, "top20DEGs_", deparse(substitute(result_table)), ".csv")
  write.csv2(top20DEG_annotation_simple, path_name, row.names = T)
  
  return(top20DEG_annotation)
}

###LOOK FOR INSPIRATION HERE: https://mkempenaar.github.io/gene_expression_analysis/chapter-5.html


#Alexis Marshal made this theme function #, base_family="Courier"
theme_Publication <- function(base_size=10) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size) #, base_family=base_family
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.2, "cm"),
            legend.spacing = unit(1, "cm"),
            legend.title = element_blank(),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#FFCC99", "#fdb462","#386cb0","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#E69F00", "#0072B2","#009E73", "#CC79A7","#F0E442", "#D55E00","#56B4E9")), ...)
  
}
