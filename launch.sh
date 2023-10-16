#! /bin/bash -login
#SBATCH -J snakemake-submission
#SBATCH -t 5-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=althaf.singhawansa@uhnresearch.ca

. env.sh

echo 'Running on H4H cluster'
snakemake -np --profile slurm --unlock
snakemake -p --profile slurm
#snakemake --cluster "sbatch -p {resources.partition} -t {resources.time_min} --mem={resources.mem_mb} -c {resources.cpus} -J {rule}_{wildcards}_job -o slurm_log/{rule}_{wildcards}_out -e slurm_log/{rule}_{wildcards}_err --mail-type=FAIL --mail-user= althaf.singhawansa@uhnresearch.ca " --lantency-wait 60 --jobs 100 -- max-jobs-per-second 1 -- restart-times 2 -- rerun-incomplete True -- keep-going True -- use-conda True -- default-resources  [partition='all', cpus=1, mem_mb=2000, time_min=60] -p