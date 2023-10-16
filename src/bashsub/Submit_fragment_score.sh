#!/bin/sh
#$ -cwd

#SBATCH -p all
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -J Fragment_score
#SBATCH --mem=30G
#SBATCH -t 4:00:00
#SBATCH -o slurm.%x.%N.%j.out
#SBATCH -e slurm.%x.%N.%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=althaf.singhawansa@uhnresearch.ca

input=$1
sample_name=$2
outDir=$3
reference=$4
dataDir=$5
windowDir=$6

#module load R
source ~/software/conda/bin/activate cfmedip_r

Rscript /cluster/home/asinghaw/git/cfmedipseq_pipeline/src/R/patient_fragment_score.R -b $input -s $sample_name -o $outDir -r $reference -d $dataDir -w $windowDir -p True 