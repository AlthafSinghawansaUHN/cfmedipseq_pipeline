' run_MEDIPS.R
Run MEDIPS for counts and conduct MEDIPS QC.

Usage:
    run_MEDIPS.R -b BAM -o OUTPUT -q QCOUT -p PAIRED

Options:
    -b --bam BAM                Path to input BAM file
    -o --output OUTPUT          Output path (RDS file)
    -q --qcout QCOUT            Path to output QC results of sample
    -p --paired PAIRED          Sample is paired end or single end sqeuncing based on cohort
' -> doc

if (! interactive()) {
  library(docopt)
  args <- docopt(doc, version='Run MEDIPS v 1.0')
  print(args)
} else {
  message('Running in interactive mode. Be sure to specify args manually.')
}

library(GenomicRanges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(IRanges)
library(MEDIPS)
library(tidyverse)

BIN_WIDTH = 300
allmainchrs = paste0('chr', c(1:22))
BSgenome = 'BSgenome.Hsapiens.UCSC.hg38'
paired_val = (args[['paired']] == "True")

medips_set = MEDIPS.createSet(file = args[['bam']],
                                     BSgenome = BSgenome,
                                     extend = 0,
                                     shift = 0,
                                     uniq = 1,
                                     window_size = BIN_WIDTH,
                                     paired = paired_val,
                                     chr.select = allmainchrs)

chr.select = medips_set@chr_names
window_size = window_size(medips_set)
chr_lengths = unname( seqlengths(BSgenome.Hsapiens.UCSC.hg38)[ seqnames(BSgenome.Hsapiens.UCSC.hg38@seqinfo)%in%chr.select ] )
no_chr_windows = ceiling(chr_lengths/window_size)
supersize_chr = cumsum(no_chr_windows)
chromosomes = chr.select

all.Granges.genomeVec = MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chromosomes, chr_lengths, window_size)
all.Granges.genomeVec$counts = medips_set@genome_count
all.Granges.genomeVec$cpm = (medips_set@genome_count/medips_set@number_regions)*1000000

count_df <- as.data.frame(all.Granges.genomeVec)
colnames(count_df) <- c("bin_chr","bin_start","bin_end","bin_width","strand","bin_counts","bin_cpm")
count_df = count_df[, c("bin_chr","bin_start","bin_end","bin_counts","bin_cpm")]

write_tsv(count_df, file = args[['output']],  col_names = TRUE) 

medipsenrichment <- tryCatch({
  medips_enrichment = MEDIPS.CpGenrich(file = args[['bam']],
                                       BSgenome = BSgenome,
                                       extend = 0,
                                       shift = 0,
                                       uniq = 1,
                                       paired = paired_val,
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
                                     paired = paired_val,
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
                                      paired = paired_val)

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

QCstats = data.frame(numReads_Unique_MEDIPS = medips_coverage$numberReads,
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

 
write_tsv(QCstats, file=args[['qcout']],  col_names = TRUE) #save QC metrics

rm(list = ls())
gc()