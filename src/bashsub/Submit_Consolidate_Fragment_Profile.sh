#!/bin/sh
#$ -cwd

#SBATCH -p all
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -J Submit_Consolidate_Fragment_Profile
#SBATCH --mem=2G
#SBATCH -t 10:00
#SBATCH -o slurm.%x.%N.%j.out
#SBATCH -e slurm.%x.%N.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=althaf.singhawansa@uhnresearch.ca

sbatch /cluster/home/asinghaw/git/cfmedipseq_pipeline/src/bashsub/Submit_consolidate_cohort_fragment_profile.sh /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/OCTANE/results/fragment_profile/summary /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/OCTANE/results/fragment_profile Filtered_100kb_Windows
sbatch /cluster/home/asinghaw/git/cfmedipseq_pipeline/src/bashsub/Submit_consolidate_cohort_fragment_profile.sh /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/OCTANE/results/fragment_profile/summary /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/OCTANE/results/fragment_profile Imformative_8CpG_300bp_Windows

sbatch /cluster/home/asinghaw/git/cfmedipseq_pipeline/src/bashsub/Submit_consolidate_cohort_fragment_profile.sh /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_1/results/fragment_profile/summary /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_1/results/fragment_profile Filtered_100kb_Windows
sbatch /cluster/home/asinghaw/git/cfmedipseq_pipeline/src/bashsub/Submit_consolidate_cohort_fragment_profile.sh /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_1/results/fragment_profile/summary /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_1/results/fragment_profile Imformative_8CpG_300bp_Windows

sbatch /cluster/home/asinghaw/git/cfmedipseq_pipeline/src/bashsub/Submit_consolidate_cohort_fragment_profile.sh /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_2/results/fragment_profile/summary /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_2/results/fragment_profile Filtered_100kb_Windows
sbatch /cluster/home/asinghaw/git/cfmedipseq_pipeline/src/bashsub/Submit_consolidate_cohort_fragment_profile.sh /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_2/results/fragment_profile/summary /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_2/results/fragment_profile Imformative_8CpG_300bp_Windows

sbatch /cluster/home/asinghaw/git/cfmedipseq_pipeline/src/bashsub/Submit_consolidate_cohort_fragment_profile.sh /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_3/results/fragment_profile/summary /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_3/results/fragment_profile Filtered_100kb_Windows
sbatch /cluster/home/asinghaw/git/cfmedipseq_pipeline/src/bashsub/Submit_consolidate_cohort_fragment_profile.sh /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_3/results/fragment_profile/summary /cluster/projects/scottgroup/people/althaf/pipelines/cfmedipseq_pipeline/NickCheng_PreDiagnosis_BreastCancer_CancerFree_3/results/fragment_profile Imformative_8CpG_300bp_Windows