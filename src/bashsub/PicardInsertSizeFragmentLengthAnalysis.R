#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(reshape2)

## For future work I can just simply change the inDir
## RCC Cohort

inDir=args[1]

setwd(inDir)

filelist=list.files(pattern = "*_IS_metrics.txt")

cohort_IS_table=matrix(NA, nrow=length(filelist)+1, ncol=4)
colnames(cohort_IS_table) <- c("Sample", "Mean", "Mode", "Median")
cohort_IS_table <- data.frame(cohort_IS_table)

for (i in 1:length(filelist)){
  sample <- strsplit(filelist[[i]],"_IS")[[1]][1]
  
  filelines <- readLines(filelist[[i]])
  filelines <- filelines[12:(length(filelines)-1)]
  
  isdata = matrix(NA, nrow=length(filelines),ncol=2)
  for(j in 1:length(filelines)){
    isdata[j,1] <- as.numeric(strsplit(filelines[j],"\t")[[1]][1])
    isdata[j,2]  <- as.numeric(strsplit(filelines[j],"\t")[[1]][2])
  }

  colnames(isdata) <- c("length","count")
  isdata <- as.data.frame(isdata)
  if (i == 1){
    p <- ggplot() + geom_point(data=isdata, aes(x=length,y=count))
  }
  else{
    p <- p + geom_point(data=isdata, aes(x=length,y=count))
  }
  
  
  completeCount_IS<- rep(isdata$length,isdata$count)
  median_IS <- median(completeCount_IS)
  
  N_IS <- sum(isdata$count)
  average_IS <- sum(isdata$count*isdata$length)/N_IS
  mode_IS <- isdata$length[isdata$count == max(isdata$count)]
  
  
  cohort_IS_table[i,1] <- sample
  cohort_IS_table[i,2] <- average_IS
  cohort_IS_table[i,3] <- mode_IS 
  cohort_IS_table[i,4] <- median_IS 
}

cohort_IS_table[i+1,1] <- "Average"
cohort_IS_table[i+1,2] <- sum(cohort_IS_table[1:i,2])/i
cohort_IS_table[i+1,3] <- sum(cohort_IS_table[1:i,3])/i
cohort_IS_table[i+1,4] <- sum(cohort_IS_table[1:i,4])/i

write.table(cohort_IS_table, file = "cohort_IS_table.txt", sep = "\t", col.names = TRUE, quote = FALSE)

rm(list=ls())
gc()