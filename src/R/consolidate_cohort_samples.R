' consolidate_cohort_samples.R
Consolidate the cohorts output samples into a single table of samples

Usage:
    compile_cohort_samples.R -i INPUT -o OUTPUT -d DATA

Options:
    -i --input INPUT            Path to input
    -o --output OUTPUT          Output path for the table of consolidated values
    -d --data DATA              Data output to consolidate
' -> doc

if (! interactive()) {
  library(docopt)
  args <- docopt(doc, version='consolidate cohort samples v 1.0')
  print(args)
} else {
  message('Running in interactive mode. Be sure to specify args manually.')
}

library(tidyverse)

cohort_dir = args[['input']]

ifelse(!dir.exists(args[['output']]),dir.create(args[['output']]),FALSE)

setwd(cohort_dir)

file_list <- list.files(pattern = "_output.tsv")

if(args[['data']] == "MEDIPS"){
  i <- 0
  for(file in file_list){
    sample <- separate(data.frame(filename = file),col = filename, into = c("sample", "extra"), sep = "\\.")$sample
    medips_counts <- read_tsv(file, col_types='ciiid')
    
    if(i == 0){
      complete_counts <- data.frame(
        window = paste(medips_counts$bin_chr,medips_counts$bin_start,medips_counts$bin_end, sep = "."),
        bin_counts = medips_counts$bin_counts
      )
      colnames(complete_counts) <- c("window",sample)
      complete_cpm <- data.frame(
        window = paste(medips_counts$bin_chr,medips_counts$bin_start,medips_counts$bin_end, sep = "."),
        bin_counts = medips_counts$bin_cpm
      )
      colnames(complete_cpm) <- c("window",sample)
    }else{
      old_colnames <- colnames(complete_counts)
      complete_counts <- complete_counts %>% add_column(medips_counts$bin_counts)
      colnames(complete_counts) <- c(old_colnames,sample)
      complete_cpm <- complete_cpm %>% add_column(medips_counts$bin_cpm)
      colnames(complete_cpm) <- c(old_colnames,sample)
    }
    i = i + 1
    
    rm(medips_counts)
    gc()
  } 
  
  write_tsv(complete_counts, file = paste(args[['output']], "MEDIPS_Counts.tsv", sep = ""), col_names = TRUE)
  write_tsv(complete_cpm, file = paste(args[['output']], "MEDIPS_CPM.tsv", sep = ""), col_names = TRUE)
}

if(args[['data']] == "MeDEStrand"){
  i <- 0
  for(file in file_list){
    medestrand_meth <- read_tsv(file, col_types='ciiid')
    
    sample <- separate(data.frame(filename = file),col = filename, into = c("sample", "extra"), sep = "\\.")$sample
    
    if(i == 0){
      complete_meth <- data.frame(
        window = paste(medestrand_meth$bin_chr,medestrand_meth$bin_start,medestrand_meth$bin_end, sep = "."),
        bin_methyl = medestrand_meth$bin_methyl
      )
      colnames(complete_meth) <- c("window",sample)
    }else{
      old_colnames <- colnames(complete_meth)
      complete_meth <- complete_meth %>% add_column(medestrand_meth$bin_methyl)
      colnames(complete_meth) <- c(old_colnames,sample)
    }
    i = i + 1
    
    rm(medestrand_meth)
    gc()
  } 
  
  write_tsv(complete_meth, file = paste(args[['output']], "MeDEStrand_AbsMethyl.tsv", sep = ""), col_names = TRUE)
}
