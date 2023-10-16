#!/bin/sh
#$ -cwd

#SBATCH -p all
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -J delfi_fragment_profile
#SBATCH --mem=30G
#SBATCH -t 12:00:00
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=althaf.singhawansa@uhnresearch.ca

input=$1
sample_name=$2
outDir=$3
dataDir=$4

#module load R
source ~/software/conda/bin/activate cfmedip_r

Rscript /cluster/home/asinghaw/git/cfmedipseq_pipeline/src/R/delfi_fragment_profile.R -b $input -s $sample_name -o $outDir -d $dataDir -p True 