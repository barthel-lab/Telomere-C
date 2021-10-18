library(tidyverse)

# INPUT: several clipreads.metrics.txt tables
# OUTOUT: A table describing number of two-linker, one-linker, or no-linker reads
path <- str_extract(snakemake@output[[1]],pattern = "^results.align.clipreads.\\w+")
setwd(path)
prefix <- str_extract(path,"\\w+$")

# Define input tables
mtrxName <- dir() %>% grep(pattern = "(.clipreads.metrics.txt)$", value = T)

# Define subgroups
noLnkName <-  grep(x = mtrxName, pattern = "unknown-unknown", value = T)
twoLnkName <- grep(x = mtrxName, pattern = "unknown", value = T, invert = T)
oneLnkName <- grep(x = mtrxName, pattern = "unknown", value = T) %>%
  grep(pattern = "unknown-unknown", value = T, invert = T)


# 1. Load metrics values into mtrxLst list
mtrxLst <- list()

for(tb in 1:length(mtrxName)) {
  mtb <- unlist(read.table(file = mtrxName[tb], sep = "\t")) %>% 
    grep(pattern = "^[A-Z]+.+[0-9]$", value = T)
  
  mtrxLst[[tb]] <- data.frame(as.numeric(na.omit(str_extract(mtb,"\\d+$"))))
  row.names(mtrxLst[[tb]]) <- str_extract(mtb,"\\D+")
  names(mtrxLst)[[tb]] <- mtrxName[tb]
}


# 2. Merge data by linker types, then do cbind() to get table
noLnk = 0
oneLnk = 0
twoLnk = 0
bindTb = 0

noLnk <- mtrxLst[[noLnkName]]

for(n in 1:length(oneLnkName)){
  oneLnk <- oneLnk + mtrxLst[[oneLnkName[n]]]
}

for(n in 1:length(twoLnkName)){
  twoLnk <- twoLnk + mtrxLst[[twoLnkName[n]]]
}

bindTb <- cbind(twoLnk, oneLnk, noLnk)
colnames(bindTb) <- c("Two_linker","One_linker","No_linker")

# 3. Calculate percentag of ereads and bases
# Percentage of reads
percReads <- apply(bindTb["Number of examined reads",], MARGIN = 1, 
              function(X){X / sum(bindTb["Number of examined reads",])})
percReads <- round(percReads*100,2)

# Percentage of clipped reads
percCpReads <- bindTb["Number of clipped reads",] / bindTb["Number of examined reads",]
percCpReads <- round(percCpReads*100,2)

# Percentage of bases
percBases <- apply(bindTb["Number of examined bases",], MARGIN = 1, 
              function(X){X / sum(bindTb["Number of examined bases",])})
percBases <- round(percBases*100,2)

# Percentage of clipped bases
percCpBases <- bindTb["Number of clipped bases",] / bindTb["Number of examined bases",]
percCpBases <- round(percCpBases*100,2)

# Generate report
result <- data.frame("Two_linker"=0,"One_linker"=0,"No_Linker"=0)
result[1,] <- bindTb["Number of examined reads",]
result[2,] <- percReads
result[3,] <- bindTb["Number of clipped reads",]
result[4,] <- percCpReads
result[5,] <- bindTb["Number of examined bases",]
result[6,] <- percBases
result[7,] <- bindTb["Number of clipped bases",]
result[8,] <- percCpBases
row.names(result) <- c("Number of examined reads",
                       "Percent of examined reads",
                       "Number of clipped reads",
                       "Percent of clipped reads",
                       "Number of examined bases",
                       "Percent of examined bases",
                       "Number of clipped bases",
                       "Percent of clipped bases")

write.table(result,paste0(prefix,".linkerQC.txt"),sep = "\t",quote = F,col.names = NA)
