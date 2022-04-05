' get_read_count_info.R
Obtain the count info of all samples into one table

Usage:
    get_read_count_info.R -i INPUT -o OUTPUT

Options:
    -i --input INPUT            Path to input
    -o --output OUTPUT          Output path for the table of consolidated values
' -> doc

if (! interactive()) {
  library(docopt)
  args <- docopt(doc, version='get count info v 1.0')
  print(args)
} else {
  message('Running in interactive mode. Be sure to specify args manually.')
}

library(tidyverse)

cohort_qc_flagstat_dir = args[['input']]

ifelse(!dir.exists(args[['output']]),dir.create(args[['output']]),FALSE)

setwd(cohort_qc_flagstat_dir)

aligned_files <- list.files(pattern = ".aligned.sorted.bam.flagstat.tsv")
sample_list <- separate(data.frame(filename = aligned_files),col = filename, into = c("sample", "extra"), sep = "\\.")$sample
lib_files_number <- list.files(pattern = paste(sample_list[1],"_lib",sep = ""))
read_count_info = matrix(NA, nrow=1, ncol = (5 + 2*length(lib_files_number)))
read_count_info <- data.frame(read_count_info)

lib_data_cols <- c()
for(i in 1:length(lib_files_number)){
  lib_data_cols <- c(lib_data_cols,paste("All_Reads_lib", i, sep = ""),paste("Mapped_Percent_lib", i, sep = ""))
}

colnames(read_count_info) <- c("Sample","Dup_Removed_PrimaryMap_Reads","Dup_Removed_Reads","All_Mapped_Reads","Percent_Duplicates",lib_data_cols)

for(sample in sample_list){
  flagsDupRmoved <- read.delim(paste(cohort_qc_flagstat_dir,sample,".aligned.sorted.markdup.bam.flagstat.tsv", sep = ""), header = FALSE)
  flags <- read.delim(paste(cohort_qc_flagstat_dir,sample,".aligned.sorted.bam.flagstat.tsv", sep = ""), header = FALSE)
  
  lib_files_list <- list.files(pattern = paste(sample,"_lib",sep = ""))
  lib_file_data <- c()
  for(lib_file in lib_files_list){
    lib_file_flag <- read.delim(paste(cohort_qc_flagstat_dir,lib_file, sep = ""), header = FALSE)
    lib_file_data <- c(lib_file_data,
                       as.numeric(lib_file_flag$V1[1]),
                       as.numeric(strsplit((lib_file_flag$V1[8]), split = "%")[[1]][1]))
  }
  
  sample_row <- c(sample,
                  as.numeric(flagsDupRmoved$V1[2]),
                  as.numeric(flagsDupRmoved$V1[1]),
                  as.numeric(flags$V1[1]),
                  (1-as.numeric(flagsDupRmoved$V1[1])/as.numeric(flags$V1[1]))*100,
                  lib_file_data)
  
  read_count_info <- rbind(read_count_info, sample_row)
}
read_count_info <- read_count_info[-1,]
write.table(read_count_info, file = paste(args[['output']], "flag_stats.tsv", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t", append = F)

rm(list = ls())
gc()