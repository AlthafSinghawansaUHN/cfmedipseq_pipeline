#! /bin/bash -login
#SBATCH -J snakemake-submission
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=althaf.singhawansa@uhnresearch.ca

. env.sh

echo 'Running on H4H cluster'
snakemake -p --profile slurm
