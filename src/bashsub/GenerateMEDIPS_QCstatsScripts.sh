#!/bin/sh
#$ -cwd

#SBATCH -p all
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 100
#SBATCH -t 01:00:00
#SBATCH -o slurm.%N.%j.out 
#SBATCH -e slurm.%N.%j.err

## Set Directories ##
#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/Justin_HN_Norm_cfDNA_PBL/samples
#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/OCTANE/samples
#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/Justin_HN_Norm_cfDNA_PBL/MEDIPSQC
#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/OCTANE/MEDIPSQC
sampleDir=$1
outDir=$2
[ ! -d $outDir ] && mkdir $outDir

bamQCDir=$outDir/bamMEDIPSQC
[ ! -d $bamQCDir ] && mkdir $bamQCDir
qsubDir=$outDir/qsub
[ ! -d $qsubDir ] && mkdir $qsubDir

cd $sampleDir

for sample in *
do

## Create MEDIPS QC .R scripts 
bamDir=$sampleDir/${sample}/merged/bwa_mem/
echo "## Load Libraries and inputs ##" > $qsubDir/${sample}_MEDIPS_QC.R
echo "library(MEDIPS)" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "library(BSgenome.Hsapiens.UCSC.hg38)\n" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "bamDir=\"$bamDir\"" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "bamQCDir=\"$bamQCDir\"" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "bamFile=\"aligned.sorted.markdup.bam\"" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "setwd(bamDir)\n" >> $qsubDir/${sample}_MEDIPS_QC.R
echo "## Set Global Values ##" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "BSgenome=\"BSgenome.Hsapiens.UCSC.hg38\"" >> $qsubDir/${sample}_MEDIPS_QC.R
echo "extend=0" >> $qsubDir/${sample}_MEDIPS_QC.R
echo "shift=0" >> $qsubDir/${sample}_MEDIPS_QC.R
echo "ws=300" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "AllMainChrs <- c(\"chr1\",\"chr2\",\"chr3\",\"chr4\",\"chr5\",\"chr6\",\"chr7\",\"chr8\",\"chr9\",\"chr10\",\"chr11\",\"chr12\",\"chr13\",\"chr14\",\"chr15\",\"chr16\",\"chr17\",\"chr18\",\"chr19\",\"chr20\",\"chr21\",\"chr22\")\n" >> $qsubDir/${sample}_MEDIPS_QC.R
echo "## Initialize QC Matrix ##" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "QCstats_Mat=matrix(NA, nrow=1, ncol=16)\n" >> $qsubDir/${sample}_MEDIPS_QC.R
echo "## Run MEDIPS QC ##" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "er=MEDIPS.CpGenrich(file = bamFile, BSgenome = BSgenome, extend = extend, shift = shift, uniq = 1, chr.select=AllMainChrs, paired=TRUE)" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "print(\"Calculated CpG Enrichment\")" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "cr=MEDIPS.seqCoverage(file = bamFile, pattern = \"CG\", BSgenome = BSgenome, extend = extend, shift = shift, uniq = 1, chr.select=AllMainChrs, paired=TRUE)" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "print(\"Calculated base coverage values\")" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "sr=MEDIPS.saturation(file= bamFile, BSgenome = BSgenome, extend = extend, shift = shift, uniq = 1, window_size = ws, nit=10, nrit=1, empty_bins=TRUE, rank=FALSE, chr.select=AllMainChrs, paired=TRUE)" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "print(\"Calculated Saturation\")" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "print (\"Start Coverage\")" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "#generating the seqCoverage just on the unique reads" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "cov.level = c(0, 1, 2, 3, 4, 5)" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "cov.res = cr\$cov.res" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "numberReads = cr\$numberReads" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "numberReadsWO = cr\$numberReadsWO" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "numberReadsWO_percentage = round((numberReadsWO/numberReads * 100), digits = 2)" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "results = NULL" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "for (j in 1:length(cov.level)) {" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "    if (j == 1) {" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "        results = c(results, length(cov.res[cov.res <= cov.level[j]])/length(cov.res) * 100)" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "    }" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "    else {" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "        results = c(results, length(cov.res[cov.res > cov.level[j - 1] & cov.res <= cov.level[j]])/length(cov.res) * 100)" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "    }" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "}" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "results = c(results, length(cov.res[cov.res > cov.level[length(cov.level)]])/length(cov.res) * 100)" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "print(\"Got the results\")" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "i=1" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "QCstats_Mat[i,1]=\"${sample}\"" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "QCstats_Mat[i,2]=length(cov.res)" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "QCstats_Mat[i,3]=cr\$numberReads" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "QCstats_Mat[i,4]=er\$enrichment.score.GoGe" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "QCstats_Mat[i,5]=er\$enrichment.score.relH" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "QCstats_Mat[i,6]=results[1]" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "QCstats_Mat[i,7]=results[2]" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "QCstats_Mat[i,8]=results[3]" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "QCstats_Mat[i,9]=results[4]" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "QCstats_Mat[i,10]=results[5]" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "QCstats_Mat[i,11]=results[6]" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "QCstats_Mat[i,12]=results[7]" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "QCstats_Mat[i,13]=cr\$numberReadsWO" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "QCstats_Mat[i,14]=numberReadsWO_percentage" >> $qsubDir/${sample}_MEDIPS_QC.R 
echo -e "QCstats_Mat[i,15]=sr\$maxEstCor[2]" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "QCstats_Mat[i,16]=sr\$maxTruCor[2]" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "colnames(QCstats_Mat)=c(\"Sample\", \"numCpGSites\", \"numReads_Unique_MEDIPS\", \"EnrichmentScore_GoGe\", \"EnrichmentScore_relH\", \"Percent_CpG_Seq_Coverage_0x\", \"Percent_CpG_Seq_Coverage_1x\", \"Percent_CpG_Seq_Coverage_2x\", \"Percent_CpG_Seq_Coverage_3x\", \"Percent_CpG_Seq_Coverage_4x\", \"Percent_CpG_Seq_Coverage_5x\", \"Percent_CpG_Seq_Coverage_Over5x\", \"Reads_do_not_cover_CpG\",\"Percent_Reads_do_not_cover_CpG\", \"Estimated_Saturation_Correlation\", \"True_Saturation_Correlation\" )#relabel the columns by the QC measures\n" >> $qsubDir/${sample}_MEDIPS_QC.R
echo -e "## Output ##" >> $qsubDir/${sample}_MEDIPS_QC.R 
echo -e "setwd(bamQCDir)" >> $qsubDir/${sample}_MEDIPS_QC.R 
echo -e "QC_file_name=(\"${sample}_MEDIPS_QCstats.txt\")" >> $qsubDir/${sample}_MEDIPS_QC.R 
echo -e "write.table(QCstats_Mat, file=QC_file_name, col.names=T, row.names=F, quote=F, sep=\"\\t\", append=F)" >> $qsubDir/${sample}_MEDIPS_QC.R 

## Create Submission .sh scripts

echo -e "#!/bin/sh" > $qsubDir/${sample}_Submit_MEDIPS_QC.sh
echo -e "#$ -cwd\n" >> $qsubDir/${sample}_Submit_MEDIPS_QC.sh
echo -e "#SBATCH -p all" >> $qsubDir/${sample}_Submit_MEDIPS_QC.sh
echo -e "#SBATCH -N 1" >> $qsubDir/${sample}_Submit_MEDIPS_QC.sh
echo -e "#SBATCH -n 1" >> $qsubDir/${sample}_Submit_MEDIPS_QC.sh
echo -e "#SBATCH --mem=30G" >> $qsubDir/${sample}_Submit_MEDIPS_QC.sh
echo -e "#SBATCH -t 12:00:00" >> $qsubDir/${sample}_Submit_MEDIPS_QC.sh
echo -e "#SBATCH -o slurm.%N.%j.out" >> $qsubDir/${sample}_Submit_MEDIPS_QC.sh 
echo -e "#SBATCH -e slurm.%N.%j.err\n" >> $qsubDir/${sample}_Submit_MEDIPS_QC.sh
echo -e "module load R/3.5.0" >> $qsubDir/${sample}_Submit_MEDIPS_QC.sh 
echo -e "R CMD BATCH ${sample}_MEDIPS_QC.R  ${sample}_MEDIPS_QC.Rout" >> $qsubDir/${sample}_Submit_MEDIPS_QC.sh

done


