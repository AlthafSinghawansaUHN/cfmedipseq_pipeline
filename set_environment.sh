. env.sh

#edit config.yaml before running to deactivate all paths

snakemake --profile slurm -p --use-conda --conda-create-envs-only -j1 conda_env/environment_is_set