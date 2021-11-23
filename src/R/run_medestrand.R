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
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyr)


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

absolute_methylation = MeDEStrand.binMethyl(MSetInput = methylset, CSet = CS, Granges = TRUE)

absolute_methylation_df <- as.data.frame(absolute_methylation)

write_tsv(absolute_methylation_df, file = args[['output']], col_names = TRUE)
