library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)

#arg <- commandArgs(trailingOnly=TRUE)
# Dev
arg <- c("BJ.rawCounts.txt", "IMR90.rawCounts.txt", "WI38.rawCounts.txt",
         "HeLa.rawCounts.txt",  "U2OS.rawCounts.txt", "VA13.rawCounts.txt")
output_file <- "enrichment.norm.png"

plotE <-function(x){

  input_name <- "input"
  TeloC_name <- "TeloC"
  file_number <- length(arg)
  bar_colors <- colorspace::sequential_hcl(10, palette = "Blues 3")
  
  img_lst <- list()
  for(ii in 1:file_number){
    
    df_tmp <- read.delim(file = arg[ii], stringsAsFactors = F)
    df_tmp["ratio"] <- df_tmp$featureReadCount / df_tmp$totalReadCount
    df_tmp["norm"] <- df_tmp$ratio / df_tmp$ratio[df_tmp$file == input_name]
    df_tmp_sub <- subset(df_tmp, df_tmp$file == TeloC_name)
      
    s_title <- stringr::str_extract(arg[[ii]],pattern = "\\w{2,5}")
    img_lst[[ii]] <- ggplot(df_tmp_sub, aes(x=featureType, y=norm)) + 
      geom_bar(stat="identity", fill=bar_colors[[ii]]) + 
      ylab("fold of reads in TeloC") + xlab("") + labs(title = s_title) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      geom_text(aes(y=norm ,label=round(norm,1)), vjust=2, color="white", size=2)
  }
  
  return(img_lst)
}

imgs <- tryCatch(plotE(arg))
export_cols <-  round(sqrt(length(imgs))) + 1
ggsave(filename = output_file,
       arrangeGrob(grobs = imgs, ncol = export_cols),
       width = 16 , height = 8, unit="in")
  