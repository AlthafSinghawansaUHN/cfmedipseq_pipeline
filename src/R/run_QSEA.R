' run_QSEA.R
Run QSEA for counts and beta value estimation and conduct MEDIPS QC.

Usage:
    run_QSEA.R -s SAMPLE -c CHROM -b BAM -o OUTPUT --count Count --beta BETA --qc QCOut [ --group GROUP ]

Options:
    -s --sample SAMPLE          Name of sample
    -c --chrom CHROM            Chromosome
    -b --bam BAM                Path to input BAM file
    -o --output OUTPUT          Output path
    
    --count Count               Output path for count data
    --beta BETA                 Output path for beta methylation estimate
    --qc QCOut                  Output path for qc matrix
    
    --group GROUP               Optional input of whether sample belongs to a group,
                                  such as "treatment" or "control"
' -> doc

if (! interactive()) {
  library(docopt)
  args <- docopt(doc, version='Running QSEA')
  print(args)
} else {
  message('Running in interactive mode. Be sure to specify args manually.')
}

library(GenomicRanges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(IRanges)
library(qsea)
library(tidyverse)
library(BiocParallel)

register(MulticoreParam(workers=4))

BIN_WIDTH = 300
chrom = args[['chrom']]
BSgenome = 'BSgenome.Hsapiens.UCSC.hg38'
mapq = 30 

if (!is.null(args[['group']])) {
  sample_group = args[['group']]
} else {
  sample_group = "unspecified"
}

sample_info <- data.frame(
  sample_name = args[['sample']],
  file_name = args[['bam']],
  group = sample_group
)

qseaset <- createQseaSet(
  sampleTable = sample_info,
  BSgenome = BSgenome,
  window_size = BIN_WIDTH,
  chr.select = chrom,
)

qseaset = addCoverage(qseaset, uniquePos = TRUE, paired = TRUE, parallel = TRUE, minMapQual = mapq)
qseaset = addPatternDensity(qseaset, "CG", name = "CpG")
qseaset = addLibraryFactors(qseaset)
qseaset = addOffset(qseaset, enrichmentPattern = "CpG")

wd = which(getRegions(qseaset)$CpG_density>1 &
           getRegions(qseaset)$CpG_density<15)
signal = (15-getRegions(qseaset)$CpG_density[wd])*.55/15+.25
signal = matrix(signal,nrow=length(signal),ncol=length(getSampleNames(qseaset)))

qseaenrichment <- tryCatch({
  qseaset = addEnrichmentParameters(
    qseaset,
    enrichmentPattern="CpG", 
    windowIdx=wd,
    signal=signal
  ) 
  return(TRUE)
}, error = function(e){
  message('Error: unable to create enrichment paramaters')
  return(FALSE)
}
)

if(qseaenrichment){
  output_beta <- makeTable(
    qseaset,
    norm_methods = c("beta"),
    samples = getSampleNames(qseaset)
  )
  output_counts <- makeTable(
    qseaset,
    norm_methods = c("counts","rpm"),
    samples = getSampleNames(qseaset)
  )
}else{
  output_counts <- makeTable(
    qseaset,
    norm_methods = c("counts","rpm"),
    samples = getSampleNames(qseaset)
  )
  output_beta <- output_counts[,1:4]
  output_beta$beta <- rep(NA, nrow(output_beta))
  colnames(output_beta) <- c("chr","window_start","window_end","CpG_density",paste(args[['sample']],"_beta",sep = ""))
}

qseaset_percentfragsbackground = getOffset(qseaset) * 100

QCstats = data.frame(numReads_Unique_QSEA = qseaset@libraries$file_name[1,"valid_fragments"],
                     QSEA_Percent_Fragments_due_Background = qseaset_percentfragsbackground,
                     QSEA_Enrichment = qseaenrichment)

## write out

setwd(args[['output']])

write_tsv(output_counts, file = args[['count']], col_names = TRUE)

write_tsv(output_beta, file = args[['beta']], col_names = TRUE)



write_tsv(QCstats, file = args[['qc']],  col_names = TRUE) #save QC metrics

if(qseaenrichment){
  png(file = paste("EnrichmentProfile",args[['chrom']],".png", sep = ""), width = 480, height = 480, units = "px")
  plotEPmatrix(qseaset)
  dev.off()
}

rm(list = ls())
gc()