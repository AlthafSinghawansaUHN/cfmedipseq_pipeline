' run_QSEA_MEDIPSQC.R
Run QSEA for counts and beta value estimation and conduct MEDIPS QC.

Usage:
    run_QSEA.R -s SAMPLE -b BAM -o OUTPUT -qc QCOUT [ --group GROUP ]

Options:
    -s --sample SAMPLE          Name of sample
    -b --bam BAM                Path to input BAM file
    -o --output OUTPUT          Output path (RDS file)
    -qc --qcout QCOUT           Path to output QC results of sample
    
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
library(qsea)
library(MEDIPS)
library(tidyverse)

BIN_WIDTH = 300
allmainchrs = paste0('chr', c(1:22))
BSgenome = 'BSgenome.Hsapiens.UCSC.hg38'

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
  chr.select = allmainchrs,
)

qseaset = addCoverage(qseaset, uniquePos = TRUE, paired = TRUE)
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
  output_counts <- makeTable(
    qseaset,
    norm_methods = c("counts","rpm","beta"),
    samples = getSampleNames(qseaset)
  )
}else{
  output_counts <- makeTable(
    qseaset,
    norm_methods = c("counts","rpm"),
    samples = getSampleNames(qseaset)
  )
}

write_tsv(output_counts, file = args[['output']], col_names = TRUE)

qseaset_percentfragsbackground = getOffset(qseaset) * 100

## MEDIPS QC

medipsenrichment <- tryCatch({
  medips_enrichment = MEDIPS.CpGenrich(file = args[['bam']],
                                       BSgenome = BSgenome,
                                       extend = 0,
                                       shift = 0,
                                       uniq = 1,
                                       paired = TRUE,
                                       chr.select = allmainchrs)
  return(TRUE)
}, error = function(e){
  message('Error: unable to create medips enrichment paramaters')
  return(FALSE)
}
)

medips_coverage = MEDIPS.seqCoverage(file = args[['bam']],
                                     pattern = "CG",
                                     BSgenome = BSgenome,
                                     extend = 0,
                                     shift = 0,
                                     uniq = 1,
                                     paired = TRUE,
                                     chr.select = allmainchrs)

medips_saturation = MEDIPS.saturation(file= args[['bam']],
                                      BSgenome = BSgenome,
                                      extend = 0,
                                      shift = 0,
                                      uniq = 1,
                                      window_size = BIN_WIDTH,
                                      nit = 10,
                                      nrit = 1,
                                      empty_bins = TRUE,
                                      rank = FALSE,
                                      chr.select = allmainchrs,
                                      paired = TRUE)

#generating the seqCoverage just on the unique reads
cov.level = c(0, 1, 2, 3, 4, 5)
cov.res = medips_coverage$cov.res
numberReads = medips_coverage$numberReads
numberReadsWO = medips_coverage$numberReadsWO
numberReadsWO_percentage = round((numberReadsWO/numberReads * 100), digits = 2)

results = NULL
for (j in 1:length(cov.level)) {
  if (j == 1) {
    results = c(results, length(cov.res[cov.res <= cov.level[j]])/length(cov.res) * 100)
  }
  else {
    results = c(results, length(cov.res[cov.res > cov.level[j - 1] & cov.res <= cov.level[j]])/length(cov.res) * 100)
  }
}
results = c(results, length(cov.res[cov.res > cov.level[length(cov.level)]])/length(cov.res) * 100)

if(medipsenrichment){
  MEDIPS_EnrichmentScore_GoGe = medips_enrichment$enrichment.score.GoGe
  MEDIPS_EnrichmentScore_relH = medips_enrichment$enrichment.score.relH
}else{
  MEDIPS_EnrichmentScore_GoGe = NA
  MEDIPS_EnrichmentScore_relH = NA
}

QCstats = data.frame(Sample = args[['sample']], 
                     numReads_Unique_QSEA = qseaset@libraries$file_name[1,"valid_fragments"],
                     QSEA_Percent_Fragments_due_Background = qseaset_percentfragsbackground,
                     QSEA_Enrichment = qseaenrichment,
                     numReads_Unique_MEDIPS = medips_coverage$numberReads,
                     MEDIPS_Enrichment = medipsenrichment,
                     EnrichmentScore_GoGe = MEDIPS_EnrichmentScore_GoGe,
                     EnrichmentScore_relH = MEDIPS_EnrichmentScore_relH,
                     Percent_CpG_Seq_Coverage_0x = results[1],
                     Percent_CpG_Seq_Coverage_1x = results[2],
                     Percent_CpG_Seq_Coverage_2x = results[3],
                     Percent_CpG_Seq_Coverage_3x = results[4],
                     Percent_CpG_Seq_Coverage_4x = results[5],
                     Percent_CpG_Seq_Coverage_5x = results[6],
                     Percent_CpG_Seq_Coverage_Over5x = results[7], 
                     Reads_do_not_cover_CpG = medips_coverage$numberReadsWO,
                     Percent_Reads_do_not_cover_CpG = numberReadsWO_percentage,
                     Estimated_Saturation_Correlation = medips_saturation$maxEstCor[2],
                     True_Saturation_Correlation = medips_saturation$maxTruCor[2])

setwd(args[['qcout']])

write_tsv(QCstats, file="QCStats_matrix.tsv",  col_names = TRUE) #save QC metrics

if(qseaenrichment){
  png(file = "EnrichmentProfile.png", width = 480, height = 480, units = "px")
  plotEPmatrix(qseaset)
  dev.off()
}

rm(list = ls())
gc()