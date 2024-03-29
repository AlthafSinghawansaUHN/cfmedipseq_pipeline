# file: 00-filtered_regions.R
# author: Derek Wong, Ph.D
# date: June 16th, 2021

library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(tidyverse)
library(httr)
library(plyr)

setwd("F:/Bratman Work/Scripts/cfmedipseq_pipeline/assets/filtered_regions")

### Generate gaps file
genome <- "hg38"
mySession <- browserSession()
genome(mySession) <- genome
gaps <- getTable(ucscTableQuery(mySession, track="gap"))
centromeres <- read_tsv("GRCh38_centromeres.bed",
                            col_names=c("chrom", "chromStart", "chromEnd", "band", "type"))
centromeres$type <- c("centromere")
centromeres <- centromeres[, c("chrom", "chromStart", "chromEnd", "type")]
gaps$chromStart <- gaps$chromStart + 1
centromeres$chromStart <- centromeres$chromStart + 1

gaps <- rbind.fill(gaps, centromeres)
gaps.hg38 <- GRanges(gaps$chrom, IRanges(gaps$chromStart,
                                         gaps$chromEnd),
                     type=gaps$type)
gaps.hg38 <- keepSeqlevels(gaps.hg38, paste0("chr", c(1:22, "X", "Y")),
                           pruning.mode="coarse")
hsapiens <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
seqinfo(gaps.hg38) <- seqinfo(hsapiens)[seqlevels(gaps.hg38),]
save(gaps.hg38, file="gaps.hg38.rda")

### Generate blacklist file
#hg38 blacklist file can be found: "https://www.encodeproject.org/files/ENCFF356LFX/"
#hg38 blacklist file can be found: "https://github.com/Boyle-Lab/Blacklist/tree/master/lists"
#hg19 bkacklist file can be found: "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz"))
blacklisted.tib <- read_tsv("GRCh38_unified_blacklist.bed",
                            col_names=c("seqnames", "start", "end"))
blacklisted.tib <- blacklisted.tib %>% mutate(start=start+1)
filters.hg38 <- makeGRangesFromDataFrame(blacklisted.tib,
                                         keep.extra.columns=TRUE)
filters.hg38 <- keepSeqlevels(filters.hg38, paste0("chr", c(1:22, "X", "Y")),
                              pruning.mode="coarse")
seqinfo(filters.hg38) <- seqinfo(Hsapiens)[seqlevels(filters.hg38),]
save(filters.hg38, file="filters.hg38.rda")

### Generate VNTR list file
VNTRs.tib <- read_tsv("VNTRs.hg38.bed",
                      col_names = c("seqnames", "start", "end", "sequence"))
VNTRs.tib <- VNTRs.tib %>% mutate(start=start)
VNTRs.hg38 <- makeGRangesFromDataFrame(VNTRs.tib,
                                       keep.extra.columns=FALSE)
VNTRs.hg38 <- keepSeqlevels(VNTRs.hg38, paste0("chr", c(1:22)),
                            pruning.mode="coarse")
seqinfo(VNTRs.hg38) <- seqinfo(Hsapiens)[seqlevels(VNTRs.hg38),]
save(VNTRs.hg38, file="VNTRs.hg38.rda")
