#!/bin/sh
#$ -cwd

#SBATCH -p all
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH -t 1-0:00:00
#SBATCH -o slurm.%N.%j.out 
#SBATCH -e slurm.%N.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=althaf.singhawansa@uhnresearch.ca

#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/OCTANE/results/bam_dedup
#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/OCTANE/results/delfi
#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_1/results/bam_dedup
#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_1/results/delfi
#/cluster/projects/scottgroup/people/althaf/assets/filtered_regions

bamDir=$1
outDir=$2
dataDir=$3
[ ! -d $outDir ] && mkdir $outDir

cd $bamDir

# this is the correct way to run it
for bam in *.bam
do

sample=$(echo $bam | awk '{split($0,arr,".aligned");print arr[1]}')

file=$bamDir/$bam

## running medips
sbatch /cluster/home/asinghaw/git/cfmedipseq_pipeline/src/bashsub/Submit_delfi_fragment_profile.sh $file $sample $outDir $dataDir

done

#sbatch Submit_all_delfi_fragment_profile.sh /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_1/results/bam_dedup /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_1/results/delfi /cluster/projects/scottgroup/people/althaf/assets/gaps_filters