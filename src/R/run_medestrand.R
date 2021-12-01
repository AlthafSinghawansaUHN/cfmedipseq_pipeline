' run_medestrand.R
Run MeDEStrand.

Usage:
    run_medestrand.R -b BAM -o OUTPUT [ -m MEDESTRAND ]

Options:
    -b --bam BAM                Path to input BAM file
    -o --output OUTPUT          Output path (RDS file)
    -m --medestrand MEDESTRAND  Path to MeDEStrand Package
' -> doc

if (! interactive()) {
    library(docopt)
    args <- docopt(doc, version='Running MeDEStrand')
    print(args)
} else {
    message('Running in interactive mode. Be sure to specify args manually.')
}

if (is.null(args[['medestrand']])) {
    library(MeDEStrand)
} else {
    devtools::load_all(args[['medestrand']])
}
library(GenomicRanges)
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)


BIN_WIDTH = 300
allmainchrs = paste0('chr', c(1:22))
BSgenome = 'BSgenome.Hsapiens.UCSC.hg38'

methylset <- MeDEStrand.createSet(
    file = args[['bam']],
    BSgenome = BSgenome,
    uniq = 1,
    extend = 0,
    shift = 0,
    window_size = BIN_WIDTH,
    chr.select = allmainchrs,
    paired = T
)

CS = MeDEStrand.countCG(pattern='CG', refObj=methylset)

absolute_methylation <- MeDEStrand.binMethyl(MSetInput = methylset, CSet = CS, Granges = FALSE)

MSet = methylset[[1]]
chr.select = MSet@chr_names
window_size = window_size(MSet)
chr_lengths = unname( seqlengths(BSgenome.Hsapiens.UCSC.hg38)[ seqnames(BSgenome.Hsapiens.UCSC.hg38@seqinfo)%in%chr.select ] )
no_chr_windows = ceiling(chr_lengths/window_size)
supersize_chr = cumsum(no_chr_windows)
chromosomes=chr.select

all.Granges.genomeVec = MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chromosomes, chr_lengths, window_size)
all.Granges.genomeVec$CF = CS@genome_CF
all.Granges.genomeVec$binMethyl= absolute_methylation

absolute_methylation_df <- as.data.frame(all.Granges.genomeVec)
colnames(absolute_methylation_df) <- c("bin_chr","bin_start","bin_end","bin_width","strand","cpg_count","bin_methyl")
absolute_methylation_df = absolute_methylation_df[, c("bin_chr","bin_start","bin_end","cpg_count","bin_methyl")]

write_tsv(absolute_methylation_df, file = args[['output']], col_names = TRUE)

rm(list = ls())
gc()