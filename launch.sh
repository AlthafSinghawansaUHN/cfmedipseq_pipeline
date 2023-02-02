#! /bin/bash -login
#SBATCH -J snakemake-submission
#SBATCH -t 5-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=althaf.singhawansa@uhnresearch.ca

. env.sh

echo 'Running on H4H cluster'
snakemake -np --profile slurm --unlock
snakemake -p --profile slurm