#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

outDir=args[1]
bamQCDir=args[2]
setwd(bamQCDir)

qcList=list.files(pattern = "*_MEDIPS_QCstats.txt$")

qcStats <- data.frame()

for (i in 1:length(qcList)){
  qcFile=qcList[[i]]
  qc <- read.delim(qcFile,sep = '\t')
  qc <- data.frame(qc)
    
  qcStats <- rbind(qcStats,qc)
}
colnames(qcStats) <- c("Sample", "numCpGSites", "numReads_Unique_MEDIPS", "EnrichmentScore_GoGe", "EnrichmentScore_relH", "Percent_CpG_Seq_Coverage_0x", "Percent_CpG_Seq_Coverage_1x", "Percent_CpG_Seq_Coverage_2x", "Percent_CpG_Seq_Coverage_3x", "Percent_CpG_Seq_Coverage_4x", "Percent_CpG_Seq_Coverage_5x", "Percent_CpG_Seq_Coverage_Over5x", "Reads_do_not_cover_CpG","Percent_Reads_do_not_cover_CpG", "Estimated_Saturation_Correlation", "True_Saturation_Correlation" )#relabel the columns by the QC measures

setwd(outDir)

filename="Consolidated_MEDIPS_QC.txt"
write.table(qcStats, file=filename, col.names=T, row.names=F, quote=F, sep="\t", append=F)

rm(list=ls())
gc()