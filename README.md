# cfmedipseq_pipeline
Post-processing pipeline for next-generation circulating methylome data generated by cfMeDIP-seq

## Dependencies

The only up-front dependency is Anaconda.

Key Anaconda package dependencies:

- pyyaml
- snakemake (bioconda)
- picard (bioconda)

Example environment is provided as conda_env.yaml

## Snakemake profiles

For a guide on how to create a Snakemake profile for your cluster setup, see https://www.sichong.site/2020/02/25/snakemake-and-slurm-how-to-manage-workflow-with-resource-constraint-on-hpc/

Create a slurm_log directory when running snakemake to place log files based on the out location of the config.yaml, example config.yaml is provided as slurm/config.yaml

place the config.yaml file withing /cluster/home/username/.config/snakemake/slurm

# Running

To run the pipeline you must activate a base conda environment containing the key anaconda dependencies, this is stored under the name "pipeline-cfMeDIP-core", if you want to use the same name and environement it can be done so by installing that environment using

$ cd MEDIPIPE
$ conda activate base
$ mamba env create --file conda_env.yaml

The launch.sh and set environment bash files sources conda to activate this using the env.sh bash file, so if you need to change path dependencies and names edit those files

Prior to running will require setting of environments on SLURM with `bash set_environment.sh`

To run the pipeline on SLURM, submit the launch.sh file with `sbatch launch.sh`.
