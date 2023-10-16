#!/bin/sh
#$ -cwd

#SBATCH -p all
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J Submit_fragment_score_scripts
#SBATCH --mem=4G
#SBATCH -t 1:00:00
#SBATCH -o slurm.%x.%N.%j.out
#SBATCH -e slurm.%x.%N.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=althaf.singhawansa@uhnresearch.ca

#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/OCTANE/results/bam_dedup
#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/OCTANE/results/fragment_score
#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_1/results/bam_dedup
#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_1/results/fragment_score
#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/Justin_HN_Norm_cfDNA_PBL/results/bam_dedup
#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/Justin_HN_Norm_cfDNA_PBL/results/fragment_score
#/cluster/projects/scottgroup/people/althaf/assets/fragment_score_refences/vessies_reference_set.txt
#/cluster/projects/scottgroup/people/althaf/assets/filtered_regions
#/cluster/projects/scottgroup/people/althaf/assets/windows_granges

bamDir=$1
outDir=$2
reference=$3
dataDir=$4
windowDir=$5

[ ! -d $outDir ] && mkdir $outDir

cd $bamDir

# this is the correct way to run it
for bam in *.bam
do

sample=$(echo $bam | awk '{split($0,arr,".aligned");print arr[1]}')

file=$bamDir/$bam

## running medips
sbatch /cluster/home/asinghaw/git/cfmedipseq_pipeline/src/bashsub/Submit_fragment_score.sh $file $sample $outDir $reference $dataDir $windowDir

done

#sbatch Submit_all_fragment_score.sh /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/OCTANE/results/bam_dedup /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/OCTANE/results/fragment_score /cluster/projects/scottgroup/people/althaf/assets/fragment_score_refences/vessies_reference_set.txt /cluster/projects/scottgroup/people/althaf/assets/filtered_regions /cluster/projects/scottgroup/people/althaf/assets/windows_granges
#sbatch Submit_all_fragment_score.sh /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_1/results/bam_dedup /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_1/results/fragment_score /cluster/projects/scottgroup/people/althaf/assets/fragment_score_refences/vessies_reference_set.txt /cluster/projects/scottgroup/people/althaf/assets/filtered_regions /cluster/projects/scottgroup/people/althaf/assets/windows_granges
#sbatch Submit_all_fragment_score.sh /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_2/results/bam_dedup /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_2/results/fragment_score /cluster/projects/scottgroup/people/althaf/assets/fragment_score_refences/vessies_reference_set.txt /cluster/projects/scottgroup/people/althaf/assets/filtered_regions /cluster/projects/scottgroup/people/althaf/assets/windows_granges
#sbatch Submit_all_fragment_score.sh /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_3/results/bam_dedup /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_3/results/fragment_score /cluster/projects/scottgroup/people/althaf/assets/fragment_score_refences/vessies_reference_set.txt /cluster/projects/scottgroup/people/althaf/assets/filtered_regions /cluster/projects/scottgroup/people/althaf/assets/windows_granges
#sbatch Submit_all_fragment_score.sh /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/Justin_HN_Norm_cfDNA_PBL/results/bam_dedup /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/Justin_HN_Norm_cfDNA_PBL/results/fragment_score /cluster/projects/scottgroup/people/althaf/assets/fragment_score_refences/vessies_reference_set.txt /cluster/projects/scottgroup/people/althaf/assets/filtered_regions /cluster/projects/scottgroup/people/althaf/assets/windows_granges