#!/bin/sh
#$ -cwd

#SBATCH -p all
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=30G
#SBATCH -t 10:00:00
#SBATCH	-o slurm.%N.%j.out
#SBATCH	-o slurm.%N.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=althaf.singhawansa@uhnresearch.ca

module load bwa/0.7.15

cp /cluster/tools/data/genomes/Arabidopsis/TAIR10/iGenomes/Sequence/Chromosomes/1.fa /cluster/projects/scottgroup/people/althaf/assets/reference/genomes/arabidopsis/Arabidopsis1.fa
cp /cluster/tools/data/genomes/Arabidopsis/TAIR10/iGenomes/Sequence/Chromosomes/3.fa /cluster/projects/scottgroup/people/althaf/assets/reference/genomes/arabidopsis/Arabidopsis3.fa


cd /cluster/projects/scottgroup/people/althaf/assets/reference/genomes/hg38_arabidopsis_chr_1_3/
cat /cluster/tools/data/genomes/human/hg38/iGenomes/Sequence/WholeGenomeFasta/genome.fa  /cluster/projects/scottgroup/people/althaf/assets/reference/genomes/arabidopsis/Arabidopsis1.fa /cluster/projects/scottgroup/people/althaf/assets/reference/genomes/arabidopsis/Arabidopsis3.fa > hg38_arabidopsis_chr_1_3.fa

bwa index hg38_arabidopsis_chr_1_3.fa hg38_arabidopsis_chr_1_3
mkdir /cluster/projects/scottgroup/people/althaf/assets/reference/genomes/hg38_arabidopsis_chr_1_3/BWA_index
mv hg38_arabidopsis_chr_1_3.fa.* /cluster/projects/scottgroup/people/althaf/assets/reference/genomes/hg38_arabidopsis_chr_1_3/BWA_index