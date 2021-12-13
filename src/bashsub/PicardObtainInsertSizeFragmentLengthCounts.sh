#!/bin/sh
#$ -cwd

#SBATCH -p all
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=30G
#SBATCH -t 1-0:00:00
#SBATCH -o slurm.%N.%j.out 
#SBATCH -e slurm.%N.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=althaf.singhawansa@uhnresearch.ca

#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/OCTANE/samples
#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/OCTANE/picardIS

sampleDir=$1
outDir=$2
[ ! -d $outDir ] && mkdir $outDir

cd $sampleDir

## load the modules
module load picard samtools R

# this is the correct way to run it
for sample in *
do

## Create MEDIPS QC .R scripts 
file=$sampleDir/${sample}/merged/bwa_mem/aligned.sorted.markdup.bam

## running picard tools
java -jar $picard_dir/picard.jar CollectInsertSizeMetrics I=$file O=$outDir/${sample}_IS_metrics.txt H=$outDir/${sample}_IS_hitstogram.pdf M=0.5

done

Rscript /cluster/home/asinghaw/git/cfmedipseq_pipeline/src/bashsub/PicardInsertSizeFragmentLengthAnalysis.R $outDir