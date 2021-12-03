' consolidate_cohort_samples.R
Fit bins from cfMeDIP-seq data using negative binomial regression.

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

setwd(paste(cohort_dir,"/samples",sep=""))

sample_list <- list.dirs(full.names = FALSE, recursive = FALSE)

if(args[['data']] == "MEDIPS"){
  i <- 0
  for(sample in sample_list){
    medips_file <- paste(cohort_dir,"/samples/",sample,"/merged/MEDIPS/medips_output.tsv", sep = "")
    medips_counts <- read_tsv(medips_file, col_types='ciiid')
    
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
  for(sample in sample_list){
    medestrand_file <- paste(cohort_dir,"/samples/",sample,"/merged/MeDEStrand/medestrand_output.tsv", sep = "")
    medestrand_meth <- read_tsv(medestrand_file, col_types='ciiid')
    
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
