#!/bin/sh
#$ -cwd

#SBATCH -p all
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 24:00:00
#SBATCH	-o slurm.%N.%j.out
#SBATCH	-o slurm.%N.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=althaf.singhawansa@uhnresearch.ca

#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/Justin_HN_Norm_cfDNA_PBL/MEDIPSQC
#/cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/OCTANE/MEDIPSQC

module load R
Rscript MEDIPS_QC_2Table.R $1