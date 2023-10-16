library(GenomicAlignments)
library(tidyverse)

## load gaps
load("D:/Bratman Work/Scripts/cfmedipseq_pipeline/assets/filtered_regions/gaps.hg38.rda")

Hsapiens <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
AB <- read.table("D:/Bratman Work/Scripts/cfmedipseq_pipeline/assets/filtered_regions/hg38_tiles.bed", col.names = c("chrom", "chromStart", "chromEnd", "Seqlength"))
AB$chromStart <- AB$chromStart + 1
AB$chromEnd <- AB$chromEnd + 1

## filter out gaps
AB <- makeGRangesFromDataFrame(AB, keep.extra.columns=TRUE)

chromosomes <- GRanges(paste0("chr", 1:22),
                       IRanges(1, seqlengths(Hsapiens)[1:22]))
tcmeres <- gaps.hg38[grepl("centromere|telomere", gaps.hg38$type)]

## remove centromere|telomere
arms <- GenomicRanges::setdiff(chromosomes, tcmeres)

arms <- arms[-c(25,27,29,41,43)]   ## removing specific arms
armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")

arms$arm <- armlevels
AB <- AB[-queryHits(findOverlaps(AB, gaps.hg38))]
AB <- AB[queryHits(findOverlaps(AB, arms))]   
AB$arm <- armlevels[subjectHits(findOverlaps(AB, arms))]

seqinfo(AB) <- seqinfo(Hsapiens)[seqlevels(seqinfo(AB))]
AB <- trim(AB)

saveRDS(AB, file = "D:/Bratman Work/Scripts/cfmedipseq_pipeline/assets/windows_granges/Filtered_100kb_Windows_Grange.rds")

DESeq2_HyperMeth_Features_top300_Annotated_Windows <- readRDS(file = "D:/Bratman Work/Scripts/OCTANE_BreastAnalysis/Results/Signatures/OCTANE_vs_Healthy_Annotated_Windows_DMR_DESeq2_Top300_Hyper_FDR_0.01.rds")

DESeq2_HyperMeth_Features_top300_Annotated_Windows$genomewindow <- rownames(DESeq2_HyperMeth_Features_top300_Annotated_Windows)
DESeq2_HyperMeth_Features_top300_Annotated_Windows_Grange <- separate(DESeq2_HyperMeth_Features_top300_Annotated_Windows,
                                                                      genomewindow,
                                                                      into = c("chrom","chromStart","chromEnd"))
DESeq2_HyperMeth_Features_top300_Annotated_Windows_Grange$chromStart <- as.numeric(DESeq2_HyperMeth_Features_top300_Annotated_Windows_Grange$chromStart)
DESeq2_HyperMeth_Features_top300_Annotated_Windows_Grange$chromEnd <- as.numeric(DESeq2_HyperMeth_Features_top300_Annotated_Windows_Grange$chromEnd)
DESeq2_HyperMeth_Features_top300_Annotated_Windows_Grange <- GRanges(DESeq2_HyperMeth_Features_top300_Annotated_Windows_Grange$chrom,
                                                                     IRanges(DESeq2_HyperMeth_Features_top300_Annotated_Windows_Grange$chromStart,
                                                                             DESeq2_HyperMeth_Features_top300_Annotated_Windows_Grange$chromEnd),
                                                                     type=DESeq2_HyperMeth_Features_top300_Annotated_Windows_Grange$type)
DESeq2_HyperMeth_Features_top300_Annotated_Windows_Grange <- keepSeqlevels(DESeq2_HyperMeth_Features_top300_Annotated_Windows_Grange, paste0("chr", c(1:22)),
                                                                           pruning.mode="coarse")
saveRDS(DESeq2_HyperMeth_Features_top300_Annotated_Windows_Grange,
        file = "D:/Bratman Work/Scripts/cfmedipseq_pipeline/assets/windows_granges/DESeq2_HyperMeth_Features_top300_Annotated_Windows_Grange.rds")


Eric_TCGA_Tumor_Specific <- read.csv("D:/Bratman Work/Scripts/OCTANE_BreastAnalysis/Results/Signatures/Eric_tcga_cancer_signature_details.csv")
Eric_TCGA_BRCA_Tumor_Specific <- Eric_TCGA_Tumor_Specific %>% filter(hyper_cohort == "TCGA-BRCA")
Eric_TCGA_BRCA_Tumor_Specific$GenomeWindow <- paste(Eric_TCGA_BRCA_Tumor_Specific$bin_chr,Eric_TCGA_BRCA_Tumor_Specific$bin_start,Eric_TCGA_BRCA_Tumor_Specific$bin_end, sep = ".")
Eric_TCGA_BRCA_Tumor_Signature <- Eric_TCGA_BRCA_Tumor_Specific[,c("hyper_cohort","immune_group","GenomeWindow")]
Eric_TCGA_BRCA_Tumor_Signature <- unique(Eric_TCGA_BRCA_Tumor_Signature)
Eric_TCGA_BRCA_Tumor_Signature_Immune <- Eric_TCGA_BRCA_Tumor_Signature %>% filter(immune_group == "Immune")
Eric_TCGA_BRCA_Tumor_Signature_Non_Immune <- Eric_TCGA_BRCA_Tumor_Signature %>% filter(immune_group == "Non-immune")
Eric_TCGA_BRCA_Tumor_Signature <- Eric_TCGA_BRCA_Tumor_Signature[!(Eric_TCGA_BRCA_Tumor_Signature$GenomeWindow %in% Eric_TCGA_BRCA_Tumor_Signature_Immune$GenomeWindow), ]
Eric_TCGA_BRCA_Tumor_Signature_Grange <- separate(Eric_TCGA_BRCA_Tumor_Signature,
                                                  GenomeWindow,
                                                  into = c("chrom","chromStart","chromEnd"))
Eric_TCGA_BRCA_Tumor_Signature_Grange$chromStart <- as.numeric(Eric_TCGA_BRCA_Tumor_Signature_Grange$chromStart)
Eric_TCGA_BRCA_Tumor_Signature_Grange$chromEnd <- as.numeric(Eric_TCGA_BRCA_Tumor_Signature_Grange$chromEnd)
Eric_TCGA_BRCA_Tumor_Signature_Grange <- GRanges(Eric_TCGA_BRCA_Tumor_Signature_Grange$chrom,
                                                 IRanges(Eric_TCGA_BRCA_Tumor_Signature_Grange$chromStart,
                                                         Eric_TCGA_BRCA_Tumor_Signature_Grange$chromEnd),
                                                 type=Eric_TCGA_BRCA_Tumor_Signature_Grange$type)
saveRDS(Eric_TCGA_BRCA_Tumor_Signature_Grange,
        file = "D:/Bratman Work/Scripts/cfmedipseq_pipeline/assets/windows_granges/Eric_TCGA_BRCA_Tumor_Signature_Grange.rds")


Nick_signature_HG38 <- readRDS("D:/Bratman Work/Scripts/OCTANE_BreastAnalysis/Results/Signatures/Nick_signature_hg38_300bp_DESeq2.rds")
allmainchrs = paste0('chr', c(1:22))
Nick_signature_HG38 <- Nick_signature_HG38[Nick_signature_HG38$chrom %in% allmainchrs,]
Nick_signature_HG38_Grange <- GRanges(Nick_signature_HG38$chrom,
                                      IRanges(Nick_signature_HG38$start,
                                              Nick_signature_HG38$end),
                                      type=Nick_signature_HG38$type)
saveRDS(Nick_signature_HG38_Grange,
        file = "D:/Bratman Work/Scripts/cfmedipseq_pipeline/assets/windows_granges/Nick_signature_HG38_Grange.rds")

Positive.windows.HER2Pos <- readRDS(file = "D:/Bratman Work/Scripts/OCTANE_BreastAnalysis/Results/SubtypeClassification/Positive.windows.HER2Pos.rds")
Positive.windows.HER2Pos <- data.frame(Windows = Positive.windows.HER2Pos)
Positive.windows.HER2Pos.Grange <- separate(Positive.windows.HER2Pos,
                                            Windows,
                                            into = c("chrom","chromStart","chromEnd"))
Positive.windows.HER2Pos.Grange$chromStart <- as.numeric(Positive.windows.HER2Pos.Grange$chromStart)
Positive.windows.HER2Pos.Grange$chromEnd <- as.numeric(Positive.windows.HER2Pos.Grange$chromEnd)
Positive.windows.HER2Pos.Grange <- GRanges(Positive.windows.HER2Pos.Grange$chrom,
                                                                     IRanges(Positive.windows.HER2Pos.Grange$chromStart,
                                                                             Positive.windows.HER2Pos.Grange$chromEnd),
                                                                     type=Positive.windows.HER2Pos.Grange$type)
Positive.windows.HER2Pos.Grange <- keepSeqlevels(Positive.windows.HER2Pos.Grange, paste0("chr", c(1:22)),
                                                                           pruning.mode="coarse")
saveRDS(Positive.windows.HER2Pos.Grange,
        file = "D:/Bratman Work/Scripts/cfmedipseq_pipeline/assets/windows_granges/Positive_windows_HER2Pos_Grange.rds")

Positive.windows.ERPos <- readRDS(file = "D:/Bratman Work/Scripts/OCTANE_BreastAnalysis/Results/SubtypeClassification/Positive.windows.ERPos.rds")
Positive.windows.ERPos <- data.frame(Windows = Positive.windows.ERPos)
Positive.windows.ERPos.Grange <- separate(Positive.windows.ERPos,
                                            Windows,
                                            into = c("chrom","chromStart","chromEnd"))
Positive.windows.ERPos.Grange$chromStart <- as.numeric(Positive.windows.ERPos.Grange$chromStart)
Positive.windows.ERPos.Grange$chromEnd <- as.numeric(Positive.windows.ERPos.Grange$chromEnd)
Positive.windows.ERPos.Grange <- GRanges(Positive.windows.ERPos.Grange$chrom,
                                           IRanges(Positive.windows.ERPos.Grange$chromStart,
                                                   Positive.windows.ERPos.Grange$chromEnd),
                                           type=Positive.windows.ERPos.Grange$type)
Positive.windows.ERPos.Grange <- keepSeqlevels(Positive.windows.ERPos.Grange, paste0("chr", c(1:22)),
                                                 pruning.mode="coarse")
saveRDS(Positive.windows.ERPos.Grange,
        file = "D:/Bratman Work/Scripts/cfmedipseq_pipeline/assets/windows_granges/Positive_windows_ERPos_Grange.rds")

Positive.windows.TNBC <- readRDS(file = "D:/Bratman Work/Scripts/OCTANE_BreastAnalysis/Results/SubtypeClassification/Positive.windows.TNBC.rds")
Positive.windows.TNBC <- data.frame(Windows = Positive.windows.TNBC)
Positive.windows.TNBC.Grange <- separate(Positive.windows.TNBC,
                                          Windows,
                                          into = c("chrom","chromStart","chromEnd"))
Positive.windows.TNBC.Grange$chromStart <- as.numeric(Positive.windows.TNBC.Grange$chromStart)
Positive.windows.TNBC.Grange$chromEnd <- as.numeric(Positive.windows.TNBC.Grange$chromEnd)
Positive.windows.TNBC.Grange <- GRanges(Positive.windows.TNBC.Grange$chrom,
                                         IRanges(Positive.windows.TNBC.Grange$chromStart,
                                                 Positive.windows.TNBC.Grange$chromEnd),
                                         type=Positive.windows.TNBC.Grange$type)
Positive.windows.TNBC.Grange <- keepSeqlevels(Positive.windows.TNBC.Grange, paste0("chr", c(1:22)),
                                               pruning.mode="coarse")
saveRDS(Positive.windows.TNBC.Grange,
        file = "D:/Bratman Work/Scripts/cfmedipseq_pipeline/assets/windows_granges/Positive_windows_TNBC_Grange.rds")

HG38_FANTOM5EnhancerPromoter_CpG_annotated_Regions_8CpGDense_UniqueWindowID <- readRDS("D:/Bratman Work/Data/Methylation Windows/HG38_FANTOM5EnhancerPromoter_CpG_annotated_Regions_8CpGDense_UniqueWindowID.rds")
HG38_FANTOM5EnhancerPromoter_CpG_annotated_Regions_8CpGDense_UniqueWindowID <- data.frame(Windows = HG38_FANTOM5EnhancerPromoter_CpG_annotated_Regions_8CpGDense_UniqueWindowID)
HG38_FANTOM5EnhancerPromoter_CpG_annotated_Regions_8CpGDense_UniqueWindowID.Grange <- separate(HG38_FANTOM5EnhancerPromoter_CpG_annotated_Regions_8CpGDense_UniqueWindowID,
                                         Windows,
                                         into = c("chrom","chromStart","chromEnd"))
HG38_FANTOM5EnhancerPromoter_CpG_annotated_Regions_8CpGDense_UniqueWindowID.Grange$chromStart <- as.numeric(HG38_FANTOM5EnhancerPromoter_CpG_annotated_Regions_8CpGDense_UniqueWindowID.Grange$chromStart)
HG38_FANTOM5EnhancerPromoter_CpG_annotated_Regions_8CpGDense_UniqueWindowID.Grange$chromEnd <- as.numeric(HG38_FANTOM5EnhancerPromoter_CpG_annotated_Regions_8CpGDense_UniqueWindowID.Grange$chromEnd)
HG38_FANTOM5EnhancerPromoter_CpG_annotated_Regions_8CpGDense_UniqueWindowID.Grange <- GRanges(HG38_FANTOM5EnhancerPromoter_CpG_annotated_Regions_8CpGDense_UniqueWindowID.Grange$chrom,
                                        IRanges(HG38_FANTOM5EnhancerPromoter_CpG_annotated_Regions_8CpGDense_UniqueWindowID.Grange$chromStart,
                                                HG38_FANTOM5EnhancerPromoter_CpG_annotated_Regions_8CpGDense_UniqueWindowID.Grange$chromEnd),
                                        type=HG38_FANTOM5EnhancerPromoter_CpG_annotated_Regions_8CpGDense_UniqueWindowID.Grange$type)
HG38_FANTOM5EnhancerPromoter_CpG_annotated_Regions_8CpGDense_UniqueWindowID.Grange <- keepSeqlevels(HG38_FANTOM5EnhancerPromoter_CpG_annotated_Regions_8CpGDense_UniqueWindowID.Grange, paste0("chr", c(1:22)),
                                              pruning.mode="coarse")
saveRDS(HG38_FANTOM5EnhancerPromoter_CpG_annotated_Regions_8CpGDense_UniqueWindowID.Grange,
        file = "D:/Bratman Work/Scripts/cfmedipseq_pipeline/assets/windows_granges/HG38_FANTOM5EnhancerPromoter_CpG_annotated_Regions_8CpGDense_UniqueWindowID_Grange.rds")
