#!/bin/sh
#$ -cwd

#SBATCH -p all
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH -t 1-0:00:00
#SBATCH -o slurm.%N.%j.out 
#SBATCH -e slurm.%N.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=althaf.singhawansa@uhnresearch.ca

#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/OCTANE/results/bam_dedup
#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/OCTANE/results/PICARD_IS
#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_1/results/bam_dedup
#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_2/results/bam_dedup
#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_3/results/bam_dedup
#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_1/results/PICARD_IS
#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_2/results/PICARD_IS
#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_3/results/PICARD_IS

bamDir=$1
outDir=$2
[ ! -d $outDir ] && mkdir $outDir

cd $bamDir

## load the modules
module load picard samtools R

# this is the correct way to run it
for bam in *.bam
do

sample=$(echo $bam | awk '{split($0,arr,".aligned");print arr[1]}')

file=$bamDir/$bam

## running picard tools
java -jar $picard_dir/picard.jar CollectInsertSizeMetrics I=$file O=$outDir/${sample}_IS_metrics.txt H=$outDir/${sample}_IS_hitstogram.pdf M=0.5

done
