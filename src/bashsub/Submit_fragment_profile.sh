#!/bin/sh
#$ -cwd

#SBATCH -p himem
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -J Fragment_profile
#SBATCH --mem=60G
#SBATCH -t 24:00:00
#SBATCH -o slurm.%x.%N.%j.out
#SBATCH -e slurm.%x.%N.%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=althaf.singhawansa@uhnresearch.ca

input=$1
sample_name=$2
outDir=$3
dataDir=$4
windows_path=$5
windows_name=$6
reference=$7

#module load R
source ~/software/conda/bin/activate cfmedip_r

Rscript /cluster/home/asinghaw/git/cfmedipseq_pipeline/src/R/fragment_profile.R -b $input -s $sample_name -p True -o $outDir -d $dataDir -w $windows_path -n $windows_name -r $reference