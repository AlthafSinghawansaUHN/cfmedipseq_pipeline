' get_duplication_rate.R
Obtain the duplication rate of all samples into one table

Usage:
    get_duplication_rate.R -i INPUT -o OUTPUT

Options:
    -i --input INPUT            Path to input
    -o --output OUTPUT          Output path for the table of consolidated values
' -> doc

if (! interactive()) {
  library(docopt)
  args <- docopt(doc, version='get duplication rate v 1.0')
  print(args)
} else {
  message('Running in interactive mode. Be sure to specify args manually.')
}

cohort_dir = args[['input']]

ifelse(!dir.exists(args[['output']]),dir.create(args[['output']]),FALSE)

setwd(paste(cohort_dir,"/samples",sep=""))

sample_list <- list.dirs(full.names = FALSE, recursive = FALSE)

numReads = matrix(NA, nrow=length(sample_list), ncol=4)
numReads <- data.frame(numReads)
colnames(numReads) <- c("Sample","Dup_Removed_Reads","All_Reads","Percent_Duplicates")
samples <- c()
numDupReads <- c()
numAllReads <- c()
duplicationRates <- c()

for(sample in sample_list){
  flagsDupRmoved <- read.delim(paste(cohort_dir,"/samples/",sample,"/merged/bwa_mem/deduped_flagstats.txt", sep = ""), header = FALSE, sep = "+")
  flags <- read.delim(paste(cohort_dir,"/samples/",sample,"/merged/bwa_mem/all_flagstats.txt", sep = ""), header = FALSE, sep = "+")
  
  samples <- c(samples, sample)
  numDupReads <- c(numDupReads, as.numeric(flagsDupRmoved$V1[9]))
  numAllReads <- c(numAllReads, as.numeric(flags$V1[9]))
  duplicationRates <- c(duplicationRates, (1-as.numeric(flagsDupRmoved$V1[9])/as.numeric(flags$V1[9]))*100)
}

numReads$Sample <- samples
numReads$Dup_Removed_Reads <- numDupReads
numReads$All_Reads <- numAllReads
numReads$Percent_Duplicates <- duplicationRates

write.table(numReads, file = paste(args[['output']], "flag_stats.tsv", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t", append = F)

rm(list = ls())
gc()