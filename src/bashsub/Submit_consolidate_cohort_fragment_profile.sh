#!/bin/sh
#$ -cwd

#SBATCH -p all
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -J Consolidate_Fragment_Profile
#SBATCH --mem=30G
#SBATCH -t 8:00:00
#SBATCH -o slurm.%x.%N.%j.out
#SBATCH -e slurm.%x.%N.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=althaf.singhawansa@uhnresearch.ca

outDir=$1
dataDir=$2
windowName=$3

[ ! -d $outDir ] && mkdir $outDir

#module load R
source ~/software/conda/bin/activate cfmedip_r

Rscript /cluster/home/asinghaw/git/cfmedipseq_pipeline/src/R/consolidate_cohort_fragment_profile.R  -o $outDir -d $dataDir -n $windowName

#sbatch /cluster/home/asinghaw/git/cfmedipseq_pipeline/src/bashsub/Submit_consolidate_cohort_fragment_profile.sh /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_3/results/fragment_profile/summary2 /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_3/results/fragment_profile Imformative_8CpG_300bp_Windows