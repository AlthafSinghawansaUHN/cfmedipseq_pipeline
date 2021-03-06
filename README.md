# cfmedipseq_pipeline
Post-processing pipeline for next-generation circulating methylome data generated by cfMeDIP-seq

## Dependencies

The only up-front dependency is Anaconda.

Key Anaconda package dependencies:

- pyyaml
- mamba (conda-forge)
- snakemake (bioconda)
- picard (bioconda)
- samtools (bioconda)
- bwa (bioconda)
- R packages: dplyr, data.table

CountsReg R package not handled by Conda and will have to be hacked into the conda environment through R install.packages() 

## Snakemake profiles

For a guide on how to create a Snakemake profile for your cluster setup, see https://www.sichong.site/2020/02/25/snakemake-and-slurm-how-to-manage-workflow-with-resource-constraint-on-hpc/

Create a slurm_log directory when running snakemake to place log files based on the out location of the config.yaml, example config.yaml is provided as snakemake_slurm_config.yaml

place the config.yaml file withing /cluster/home/username/.config/snakemake/slurm

# Running

Prior to running will require setting of environments on SLURM with `bash set_environment.sh`

To run the pipeline on SLURM, submit the launch.sh file with `sbatch launch.sh`.
