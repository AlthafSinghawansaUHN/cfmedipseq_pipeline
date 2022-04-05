#! /bin/bash -login
#SBATCH -J snakemake-submission
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=althaf.singhawansa@uhnresearch.ca

. env.sh

echo 'Running on H4H cluster'
snakemake --use-conda --profile slurm -p /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/OCTANE/samples/OCT_010565/merged/bin_stats/bin_stats_fit_nbglm.tsv 
#snakemake --use-conda --profile slurm -np /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/test/samples/4518_Norm_testsample/merged/QSEA/qsea_output.tsv
#snakemake --use-conda --profile slurm -np /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/OCTANE/samples/OCT_010565/merged/bin_stats/bin_stats_fit_nbglm.tsv
#snakemake --use-conda --profile slurm -np /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/Justin_HN_Norm_cfDNA_PBL/samples/JustinHN_2496_HN_PBL/merged/bwa_mem/all_flagstats.txt
#snakemake --use-conda --profile slurm -np /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/Justin_HN_Norm_cfDNA_PBL/samples/JustinHN_2496_HN_PBL/merged/libs_cleaned.txt
#snakemake --use-conda --profile slurm -np /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/Justin_HN_Norm_cfDNA_PBL/samples/JustinHN_2496_HN_PBL/merged/bin_stats/chrom_bin_stats_libs_cleaned.txt