import yaml
import pandas as pd
import gzip
import statistics

# ------------ #
#  Base Setup  #
# ------------ #

configfile: "config.yml"

path_to_data = config['data']['base_path']

excluded_cases = config['data']['excluded_cases']
if excluded_cases:
    print('Excluding cases:')
    for case in excluded_cases:
        print(case)
else:
    excluded_cases = []

# ---------------------------- #
#  Important Shared Functions  #
# ---------------------------- #

def get_active_cohorts():
    """Returns a list of active cohorts.
    Active cohorts are defined in config.yml data > cohorts > [cohortname] > active: True
    """
    return(c for c in config['data']['cohorts'] if config['data']['cohorts'][c]['active'])

def get_cohort_data(cohort):
    """Parses the samplesheet for a specific cohort.
    """
    
    samplesheet = pd.read_csv(config['data']['cohorts'][cohort]['samplesheet'], comment='#').drop_duplicates()
    samplesheet = samplesheet[~samplesheet.sample_name.isin(excluded_cases)]
    return samplesheet

def get_cohort_config(cohort):
    """Returns the configuration details for specific cohort.
    Importantly, retrieving a cohort using this function (rather than directly
    accessing it via the config dictionary) will also inject the default settings
    (config > data > defaults) wherever a specific setting is not specified.
    """
    
    config_data = dict(config['data']['defaults'])
    if 'settings' in config['data']['cohorts'][cohort]:
        for key in config['data']['cohorts'][cohort]['settings']:
            if key in config_data:
                config_data[key] = config['data']['cohorts'][cohort]['settings'][key]
            else:
                raise(Exception('Setting {} does not exist for cohort {}. Available settings: {}'.format(
                    key, cohort, ', '.join(config_data.keys())
                )))
    return(config_data)	

def get_fastq_path(cohort, sample, library, read_in_pair=1):
    """Retrieves the path to the fastq file.
    Keyword arguments:
        cohort -- name of the cohort whose samplesheet should be accessed.
        sample -- identifier of the sample as specified in the samplesheet.
        library -- integer representing the library index as specified in the samplesheet.
        read_in_pair -- 1 or 2 - representing read 1 or read 2 in paired end data.
    """
    
    library = int(library)
    cohort_data = get_cohort_data(cohort)
    sample_line = cohort_data[
        (cohort_data.sample_name == sample) &
        (cohort_data.library_index == library) &
        (cohort_data.read_in_pair == read_in_pair)
    ]
    return sample_line.path.to_list()[0]

def get_all_samples(cohort=None):
    """Retrieves all samples to be processed.
    Does so by calling get_cohort_data, and therefore filters out excluded_cases.
    Keyword arguments:
        cohort -- Name of a cohort, OPTIONAL. If not specified, returns all samples
                  across all cohorts.
    """
    
    all_samples = pd.concat([
        get_cohort_data(cohort_name).assign(cohort_name = cohort_name)
        for cohort_name
        in config['data']['cohorts']
        if config['data']['cohorts'][cohort_name]['active']
        ])

    if cohort is None:
        return(all_samples)
    else:
        return(all_samples[all_samples.cohort_name == cohort])

def get_all_samples_list(cohort=None):
    """Returns all samples in list format.
    By calling get_all_samples()
    """
    
    return get_all_samples(cohort).sample_name.unique().tolist()

def get_all_samples_with_cohorts():
    """Returns a list of tuples with cohort name and sample name."""
    samples = get_all_samples()[['cohort_name', 'sample_name']]
    return(zip(
        samples.drop_duplicates().cohort_name.tolist(),
        samples.drop_duplicates().sample_name.tolist()
    ))

def clean(command):
    command = command.replace('\n', ' ').replace('\t', ' ')
    while '  ' in command:
        command = command.replace('  ', ' ')
    return command.strip()

# ------------------------------ #
#  Beginning of Snakemake Rules  #
# ------------------------------ #

rule all:
    input:
        # Create aligned bam files 
        expand(
            [path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam'.format(
                    cohort=v[0],
                    sample=v[1]
                ) for v in get_all_samples_with_cohorts()
            ]
        ),
        # Run Consolidated QC
        expand(
            path_to_data + '/{cohort}/results/qc/flagstats/consolidated/flag_stats.tsv',
            cohort = get_active_cohorts()
        ),
        # Run Consolidated MEDIPS
        expand(
            path_to_data + '/{cohort}/results/MEDIPS/consolidated/MEDIPS_Counts.tsv',
            cohort = get_active_cohorts()
        ),
        # Run Consolidated MeDEStrand
        expand(
            path_to_data + '/{cohort}/results/MeDEStrand/consolidated/MeDEStrand_AbsMethyl.tsv',
            cohort = get_active_cohorts()
        ),
        #run MeD-ReMix
        #expand(
        #    path_to_data + '/samples/{sample}/merged/bin_stats/bin_stats_fit_nbglm.tsv',
        #    sample = get_all_samples_list()
        #),


# ------------------------------------------ #
#  Install conda dependencies automatically  #
# ------------------------------------------ #
# This section is only necessary to support air-gapped environments,
# where the cluster nodes have no direct access to internet.
# In this setting, you must first load all necessary packages from
# a build node with internet access. This can be done simply
# by running 
# bash set_environment.sh
# from a node with internet access, which loads install_dependencies
# rule below to install all of the necessary conda dependencies
# in advance. Any time the git repo is pulled or if you modify
# any of the conda environments in conda_env, this needs to be
# re-run.

rule set_environment:
    input:
        'conda_env/biopython_installed',
        'conda_env/cfmedip_r_installed',
        'conda_env/fastqc_installed',
        'conda_env/samtools_installed',
        'conda_env/trimgalore_installed'
    output:
        'conda_env/environment_is_set'
    shell:
        'touch {output}'

rule biopython:
    output:
        'conda_env/biopython_installed'
    conda: 'conda_env/biopython.yml'
    shell:
        'touch {output}'

rule cfmedip_r:
    output:
        'conda_env/cfmedip_r_installed'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'touch {output}'

rule fastqc:
    output:
        'conda_env/fastqc_installed'
    conda: 'conda_env/fastqc.yml'
    shell:
        'touch {output}'

rule samtools:
    output:
        'conda_env/samtools_installed'
    conda: 'conda_env/samtools.yml'
    shell:
        'touch {output}'

rule trimgalore:
    output:
        'conda_env/trimgalore_installed'
    conda: 'conda_env/trimgalore.yml'
    shell:
        'touch {output}'

# -------------------------- #
#  Pre-process input FASTQs  #
# -------------------------- #

rule gunzip_fastq:
    input:
        lambda wildcards: get_fastq_path(wildcards.cohort, wildcards.sample, int(wildcards.lib), int(wildcards.read))
    output:
        temp(path_to_data + '/{cohort}/tmp/gunzip_fastq/{sample}_lib{lib}_R{read}.fastq')
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    shell:
        'gunzip -c {input} > {output}'
		
# QC of input FASTQs using FASTQC
rule fastqc_fastq:
    input:
        path_to_data + '/{cohort}/tmp/gunzip_fastq/{sample}_lib{lib}_R{read}.fastq'
    output:
        html=path_to_data + '/{cohort}/results/qc_input/{sample}_lib{lib}_R{read}_fastqc.html',
        zipfile=path_to_data + '/{cohort}/results/qc_input/{sample}_lib{lib}_R{read}_fastqc.zip'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    params:
        outdir = lambda wildcards, output: '/'.join(output.html.split('/')[0:-1])
    conda: 'conda_env/fastqc.yml'
    shell:
        'fastqc --outdir {params.outdir} {input}'

# Extract Barcodes using ConsensusCruncher
# Pulls the path to extract_barcodes.py from config > paths > dependencies > extract_barcodes_path
# Runs extract_barcodes.py based on cohort settings config > data > defaults or
# config > [cohort] > settings
rule extract_barcodes:
    input:
        R1 =    path_to_data + '/{cohort}/tmp/gunzip_fastq/{sample}_lib{lib}_R1.fastq',
        R2 =    path_to_data + '/{cohort}/tmp/gunzip_fastq/{sample}_lib{lib}_R2.fastq',
        R1_qc =    path_to_data + '/{cohort}/results/qc_input/{sample}_lib{lib}_R1_fastqc.html',
        R2_qc =    path_to_data + '/{cohort}/results/qc_input/{sample}_lib{lib}_R2_fastqc.html'
    output:
        R1 =    temp(path_to_data + '/{cohort}/tmp/extract_barcodes/{sample}_lib{lib}/{sample}_lib{lib}_extract_barcode_R1.fastq'),
        R2 =    temp(path_to_data + '/{cohort}/tmp/extract_barcodes/{sample}_lib{lib}/{sample}_lib{lib}_extract_barcode_R2.fastq'),
        stats =   path_to_data + '/{cohort}/results/qc_input_extract_barcodes/{sample}_lib{lib}_barcode_stats.txt'
    params:
        outprefix = lambda wildcards, output: output.R1.split('_barcode_')[0],
        bpattern = lambda wildcards: get_cohort_config(wildcards.cohort)['bpattern'],
        blist = lambda wildcards: get_cohort_config(wildcards.cohort)['barcodes'],
        extract_barcodes = config['paths']['dependencies']['extract_barcodes_path'],
        barcode_stats_tmp_path = path_to_data + '/{cohort}/tmp/extract_barcodes/{sample}_lib{lib}_barcode_stats.txt',
        barcode_stats_out_path = path_to_data + '/{cohort}/results/qc_input_extract_barcodes/'
    resources: cpus=1, mem_mb=16000, time_min='1-00:00:00'
    run:
        if params.bpattern is None:
            if params.blist is None:
                shell("cp {input.R1} > {output.R1}")
                shell("cp {input.R2} > {output.R2}")
                shell("touch {output.stats}")
            else:
                shell("python {params.extract_barcodes} --read1 {input.R1} --read2 {input.R2} --outfile {params.outprefix} --blist {params.blist}")
                shell("cp {params.barcode_stats_tmp_path} {params.barcode_stats_out_path}")
        else:
            shell("python {params.extract_barcodes} --read1 {input.R1} --read2 {input.R2} --outfile {params.outprefix} --bpattern {params.bpattern} --blist {params.blist}")
            shell("cp {params.barcode_stats_tmp_path} {params.barcode_stats_out_path}")
			
# Trims FASTQ using trimgalore to remove barcode sequences
# By default, trims 10 base pairs from the 5' end, which seems to be correct for OICR cfMeDIP-seq output.
# This can be configured in the config.yml under data > cohorts > settings > trimgalore.
rule trim_fastq_paired:
    input:
        R1 = path_to_data + '/{cohort}/tmp/extract_barcodes/{sample}_lib{lib}/{sample}_lib{lib}_extract_barcode_R1.fastq',
        R2 = path_to_data + '/{cohort}/tmp/extract_barcodes/{sample}_lib{lib}/{sample}_lib{lib}_extract_barcode_R2.fastq'
    output:
        trimmed_1 = temp(path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_extract_barcode_R1_val_1.fq'),
        trimmed_2 = temp(path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_extract_barcode_R2_val_2.fq'),
        report_1 = path_to_data + '/{cohort}/results/qc/trim_report/{sample}_lib{lib}_extract_barcode_R1.fastq_trimming_report.txt',
        report_2 = path_to_data + '/{cohort}/results/qc/trim_report/{sample}_lib{lib}_extract_barcode_R2.fastq_trimming_report.txt'
    params:
        outdir = lambda wildcards, output: '/'.join(output.trimmed_1.split('/')[0:-1]),
        trimgalore_settings = lambda wildcards: get_cohort_config(wildcards.cohort)['trimgalore']
    resources: cpus=4, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/trimgalore.yml'
    shell:
        'trim_galore --cores 4 --dont_gzip --paired {params.trimgalore_settings} --output_dir {params.outdir} {input.R1} {input.R2} && cp ' + path_to_data + '/{wildcards.cohort}/tmp/trim_fastq/{wildcards.sample}_lib{wildcards.lib}_extract_barcode_R*.fastq_trimming_report.txt ' + path_to_data + '/{wildcards.cohort}/results/qc/trim_report/'

rule trim_fastq_single:
    input:
        R1 = path_to_data + '/{cohort}/tmp/gunzip_fastq/{sample}_lib{lib}_R1.fastq',
        R1_qc =    path_to_data + '/{cohort}/results/qc_input/{sample}_lib{lib}_R1_fastqc.html',
    output:
        trimmed_1 = temp(path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_R1_trimmed.fq'),
        report_1 = path_to_data + '/{cohort}/results/qc/trim_report/{sample}_lib{lib}_R1.fastq_trimming_report.txt',
    params:
        outdir = lambda wildcards, output: '/'.join(output.trimmed_1.split('/')[0:-1]),
        trimgalore_settings = lambda wildcards: get_cohort_config(wildcards.cohort)['trimgalore']
    resources: cpus=4, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/trimgalore.yml'
    shell:
        'trim_galore --cores 4 --dont_gzip {params.trimgalore_settings} --output_dir {params.outdir} {input.R1}  && cp ' + path_to_data + '/{wildcards.cohort}/tmp/trim_fastq/{wildcards.sample}_lib{wildcards.lib}_R1.fastq_trimming_report.txt ' + path_to_data + '/{wildcards.cohort}/results/qc/trim_report/'

# --------------------------- #
#  Align FASTQs to reference  #
# --------------------------- #

def get_read_group_from_fastq(fastq_file, sample_name):
    """Extracts the read group data from FASTQs automatically.
    This information is used in the shell script of rule bwa_mem.
    """
    
    with gzip.open(fastq_file, 'rt') as fastq:
        header = next(fastq)
        (instrument, run_number, flowcell, lane, tile, xpos, ypos) = header.split(' ')[0].split(':')
        lib_value = header.strip().split(' ')[1].split(':')[3]
        rg_line = r"@RG\tID:{flowcell}_{lane}\tSM:{sample}\tPL:Illumina\tPU:.\tLB:{lib_value}".format(
            flowcell = flowcell,
            lane = lane,
            sample = sample_name,
            lib_value = lib_value
        )
        return rg_line

# Run BWA mem on FASTQs after extracting barcodes and trimming for paired end data
rule bwa_mem_paired:
    input:
        R1 = path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_extract_barcode_R1_val_1.fq',
        R2 = path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_extract_barcode_R2_val_2.fq',
    output:
        temp(path_to_data + '/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.paired.bam')
    resources: cpus=4, mem_mb=16000, time_min='72:00:00'
    params:
        read_group = lambda wildcards, input: get_read_group_from_fastq(
            fastq_file = get_fastq_path(wildcards.cohort, wildcards.sample, wildcards.lib),
            sample_name = wildcards.sample
        ),
        bwa_index = lambda wildcards: get_cohort_config(wildcards.cohort)['bwa_index']
    conda: 'conda_env/samtools.yml'
    shell:
        "bwa mem -M -t4 -R'{params.read_group}' {params.bwa_index} {input.R1} {input.R2} | samtools view -bSo {output}"

#Run BWA mem on FASTQs of single end data (no umi extraction or trimming)
rule bwa_mem_single:
    input:
        R1 =    path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_R1_trimmed.fq'
    output:
        temp(path_to_data + '/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.single.bam')
    resources: cpus=4, mem_mb=16000, time_min='72:00:00'
    params:
       read_group = lambda wildcards, input: get_read_group_from_fastq(
          fastq_file = get_fastq_path(wildcards.cohort, wildcards.sample, wildcards.lib),
            sample_name = wildcards.sample
        ),
        bwa_index = lambda wildcards: get_cohort_config(wildcards.cohort)['bwa_index']
    conda: 'conda_env/samtools.yml'
    shell:
        "bwa mem -M -t4 -R'{params.read_group}' {params.bwa_index} {input.R1} | samtools view -bSo {output}"			

# Sort Bam file and remove unmapped reads
rule bam_to_sorted_bam_paired:
    input:
        path_to_data + '/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.paired.bam'
    output:
        bam = temp(path_to_data + '/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.paired.sorted.bam'),
        index = temp(path_to_data + '/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.paired.sorted.bam.bai'),
    resources: cpus=32, mem_mb=30000, time_min='72:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        clean(r'''
        samtools view -buS -f 2 -F 4 -@32 {input} |
        samtools fixmate -m - - |
        samtools sort -@32 -o {output.bam} && samtools index {output.bam} 
        ''')

rule bam_to_sorted_bam_single:
    input:
        path_to_data + '/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.single.bam'
    output:
        bam = temp(path_to_data + '/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.single.sorted.bam'),
        index = temp(path_to_data + '/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.single.sorted.bam.bai'),
    resources: cpus=32, mem_mb=30000, time_min='72:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        clean(r'''
        samtools view -buS -F 4 -@32 {input} |
        samtools fixmate -m - - |
        samtools sort -@32 -o {output.bam} && samtools index {output.bam} 
        ''')

def get_libraries_of_sample(sample):
    """Returns all library indices of a sample based on samplesheet.
    """
    
    filtered_table = get_all_samples()[get_all_samples().sample_name == sample]
    return list(set(filtered_table.library_index.to_list()))

# If there are multiple libraries for a given sample, as specified in samplesheet,
# these libraries are automatically merged at this step into a single unified BAM.
rule merge_bam_paired:
    input:
        lambda wildcards: expand(
                path_to_data + '/' + wildcards.cohort + '/tmp/bwa_mem/' + wildcards.sample + '_lib{lib}.paired.sorted.bam',
                lib=get_libraries_of_sample(wildcards.sample)
        )
    output:
        bam = temp(path_to_data + '/{cohort}/tmp/merge_bam/{sample}.paired.aligned.sorted.bam'),
        index = temp(path_to_data + '/{cohort}/tmp/merge_bam/{sample}.paired.aligned.sorted.bam.bai')
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        'samtools merge {output.bam} {input} && samtools index {output}'
		
rule merge_bam_single:
    input:
        lambda wildcards: expand(
                path_to_data + '/' + wildcards.cohort + '/tmp/bwa_mem/' + wildcards.sample + '_lib{lib}.single.sorted.bam',
                lib=get_libraries_of_sample(wildcards.sample)
        )
    output:
        bam = temp(path_to_data + '/{cohort}/tmp/merge_bam/{sample}.single.aligned.sorted.bam'),
        index = temp(path_to_data + '/{cohort}/tmp/merge_bam/{sample}.single.aligned.sorted.bam.bai')
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        'samtools merge {output.bam} {input} && samtools index {output}'

def get_paired_or_single_end_sorted_merged_bam(wildcards):
    paired_data = get_cohort_config(wildcards.cohort)['paired']
    if paired_data:
        return path_to_data + '/' + wildcards.cohort + '/tmp/merge_bam/' + wildcards.sample + '.paired.aligned.sorted.bam'
    else:
        return path_to_data + '/' + wildcards.cohort + '/tmp/merge_bam/' + wildcards.sample + '.single.aligned.sorted.bam'
	
# Bam markdup and create index. 
# This step finalizes the definitive BAM file.
rule bam_markdup:
    input:
        lambda wildcards: get_paired_or_single_end_sorted_merged_bam(wildcards)
    output:
        bam = path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam',
        index = path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam.bai'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        "samtools markdup -r {input} {output.bam} && samtools index {output.bam}"

# ----------------- #
#  QC of BAM files  #
# ----------------- #

# Run FASTQC on final BAM files.
rule fastqc_bam:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam',
    output:
        html = path_to_data + '/{cohort}/results/qc/bam_fastqc/{sample}.aligned.sorted.markdup_fastqc.html',
        zipfile = path_to_data + '/{cohort}/results/qc/bam_fastqc/{sample}.aligned.sorted.markdup_fastqc.zip'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    params:
        outdir = lambda wildcards, output: '/'.join(output.html.split('/')[0:-1])
    conda: 'conda_env/fastqc.yml'
    shell:
        'fastqc --outdir {params.outdir} {input}'
		
# Run Flagstat to get basic stats on final BAM files.
rule bam_flagstat:
    input:
        aligned = lambda wildcards: get_paired_or_single_end_sorted_merged_bam(wildcards),
        deduped = path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam'
    output:
        aligned = path_to_data + '/{cohort}/results/qc/flagstats/{sample}.aligned.sorted.bam.flagstat.tsv',
        deduped = path_to_data + '/{cohort}/results/qc/flagstats/{sample}.aligned.sorted.markdup.bam.flagstat.tsv'
    params:
        qc_path = path_to_data + '/{cohort}/results/qc/flagstats/'
    resources: cpus=1, mem_mb=8000, time_min='1:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        "samtools flagstat {input.aligned} -O tsv > {output.aligned} && samtools flagstat {input.deduped} -O tsv > {output.deduped}"

def get_paired_or_single_end_unsorted_bam(wildcards):
    paired_data = get_cohort_config(wildcards.cohort)['paired']
    if paired_data:
        return path_to_data + '/' + wildcards.cohort + '/tmp/bwa_mem/' + wildcards.sample + '_lib' + wildcards.lib + '.paired.bam'
    else:
        return path_to_data + '/' + wildcards.cohort + '/tmp/bwa_mem/' + wildcards.sample + '_lib' + wildcards.lib + '.single.bam'
		
rule unfiltered_bam_flagstat:
    input:
        lambda wildcards : get_paired_or_single_end_unsorted_bam(wildcards)
    output:
        path_to_data + '/{cohort}/results/qc/flagstats/{sample}_lib{lib}.bam.flagstat.tsv'
    resources: cpus=1, mem_mb=8000, time_min='1:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        "samtools flagstat {input} -O tsv > {output}"

def get_all_samples_with_libraries(cohort=None):
    """Returns a list of tuples with cohort name and sample name."""
    samples = get_all_samples(cohort)[['sample_name', 'library_index']]
    return(zip(
        samples.drop_duplicates().sample_name.tolist(),
        samples.drop_duplicates().library_index.tolist()
    ))
	
# Unified QC rule which runs all of the above QCs
rule consolidate_cohort_bam_qc_info:
    input:
        fastqc = lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/qc/bam_fastqc/{sample}.aligned.sorted.markdup_fastqc.html',
            sample = get_all_samples_list(wildcards.cohort)
        ),
        aligned = lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/qc/flagstats/{sample}.aligned.sorted.bam.flagstat.tsv',
            sample = get_all_samples_list(wildcards.cohort)
        ),
        deduped = lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/qc/flagstats/{sample}.aligned.sorted.markdup.bam.flagstat.tsv',
            sample = get_all_samples_list(wildcards.cohort)
        ),
        lib_stat_files = lambda wildcards: expand(
            [path_to_data + '/' + wildcards.cohort + '/results/qc/flagstats/{sample}_lib{lib}.bam.flagstat.tsv'.format(
                    sample = v[0],
                    lib = v[1]
                ) for v in get_all_samples_with_libraries(wildcards.cohort)
            ]
        )
    output:
        path_to_data + '/{cohort}/results/qc/flagstats/consolidated/flag_stats.tsv'
    params:
        in_path = path_to_data + '/{cohort}/results/qc/flagstats/',
        out_path = path_to_data + '/{cohort}/results/qc/flagstats/consolidated/'
    resources: cpus=1, mem_mb=30000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/get_read_count_info.R -i {params.in_path} -o {params.out_path}'
	
# ------------------------------------------ #
#  Compute BAM bin stats including coverage  #
# ------------------------------------------ #

def get_bsgenome_chrom(species, chrom):
    chrom_map = {'1': 'Chr1', '3': 'Chr3'}
    if species == 'human':
        return chrom
    elif species == 'arabidopsis':
        return chrom_map[chrom]

rule bam_bin_stats:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam',
    output:
        binstat = temp(path_to_data + '/{cohort}/tmp/bam_bin_stats/bin_stats_{sample}_{species}_{chrom}.tsv'),
        filtered = path_to_data + '/{cohort}/results/bam_bin_stats/removed_bins_{sample}_{species}_{chrom}.tsv',
    params:
        bsgenome = lambda wildcards: get_cohort_config(wildcards.cohort)['bsgenome'][wildcards.species],
        bsgenome_chr = lambda wildcards: get_bsgenome_chrom(wildcards.species, wildcards.chrom)
    resources: cpus=1, mem_mb=30000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        clean('''
        Rscript src/R/bin_stats.R
            -b {input}
            -g {params.bsgenome}
            -c {wildcards.chrom}
            -o {output.binstat}
            --filtered {output.filtered}
            --bsgchr {params.bsgenome_chr}
        ''')

# Preload chromosome values for below rule
chromosomes = {}
for c in config['data']['cohorts']:
    if config['data']['cohorts'][c]['active']:
        chr_data = get_cohort_config(c)['chromosomes']
        chromosome_tuples = [(species, chrom) for species in chr_data for chrom in chr_data[species].split(',')]
        chromosomes[c] = chromosome_tuples
		
rule merge_bin_stats:
    input:
        lambda wildcards: [path_to_data + '/' + wildcards.cohort + '/tmp/bam_bin_stats/bin_stats_{{sample}}_{species}_{chrom}.tsv'.format(species=a[0], chrom=a[1]) for a in chromosomes[wildcards.cohort]],
    output:
        path_to_data + '/{cohort}/results/merge_bin_stats/bin_stats_{sample}.feather'
    params:
        paths = lambda wildcards, input: ','.join(input)
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/row_bind_tables.R -p "{params.paths}" -o {output} --in-tsv --out-feather --omit-paths'
		
# ----------------------------------------------------------------------------- #
#  Fit cfMeDIP-seq coverage stats to infer absolute methylation using MedReMix  #
# ----------------------------------------------------------------------------- #

# This is the currently used negative binomial GLM approach to fitting
rule cfmedip_nbglm:
    input:
        path_to_data + '/{cohort}/results/merge_bin_stats/bin_stats_{sample}.feather'
    output:
        fit=path_to_data + '/{cohort}/results/cfmedip_nbglm/{sample}_fit_nbglm.feather',
        model=path_to_data + '/{cohort}/results/cfmedip_nbglm/{sample}_fit_nbglm_model.Rds',
    resources: cpus=1, time_min='5-00:00:00', mem_mb=lambda wildcards, attempt: 16000 if attempt == 1 else 30000
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/cfmedip_nbglm.R -i {input} -o {output.fit} --modelout {output.model}'

# Below is the full bayesian approach for fitting, which remains under development
rule fit_bin_stats:
    input:
        path_to_data + '/{cohort}/results/merge_bin_stats/bin_stats_{sample}.feather'
    output:
        path_to_data + '/{cohort}/results/fit_bin_stats/{sample}_fit_{method}.tsv'
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/fit_cpg_bins.R -i {input} -o {output} --method {wildcards.method}'

# ------------------------------------------------ #
#  MEDIPS Method of methylation count calculation  #
# ------------------------------------------------ #

rule run_MEDIPS:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam'
    output:
        medips_count = path_to_data + '/{cohort}/results/MEDIPS/{sample}.medips_output.tsv',
        medips_qc = path_to_data + '/{cohort}/results/MEDIPS/{sample}.QCStats_matrix.tsv'
    params:
	    paired_val = lambda wildcards: get_cohort_config(wildcards.cohort)['paired']
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/run_MEDIPS.R -b {input} -o {output.medips_count} -q {output.medips_qc} -p {params.paired_val}'

rule consolidate_cohort_MEDIPS:
    input:
        lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/MEDIPS/{sample}.medips_output.tsv',
            sample = get_all_samples_list(wildcards.cohort)
        )
    output:
        path_to_data + '/{cohort}/results/MEDIPS/consolidated/MEDIPS_Counts.tsv',
        path_to_data + '/{cohort}/results/MEDIPS/consolidated/MEDIPS_CPM.tsv'
    params:
        data = 'MEDIPS',
        in_path = path_to_data + '/{cohort}/results/MEDIPS/',
        out_path = path_to_data + '/{cohort}/results/MEDIPS/consolidated/'
    resources: cpus=1, mem_mb=30000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/consolidate_cohort_samples.R -i {params.in_path} -o {params.out_path} -d {params.data}'

# --------------------------------------------- #
#  MeDEStrand Method of methylation correction  #
# --------------------------------------------- #

rule run_MeDEStrand:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam'
    output:
        path_to_data + '/{cohort}/results/MeDEStrand/{sample}.medestrand_output.tsv'
    params:
        medestrand = config['paths']['dependencies']['medestrand_path'],
	    paired_val = lambda wildcards: get_cohort_config(wildcards.cohort)['paired']
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/run_medestrand.R -b {input} -o {output} -p {params.paired_val} -m {params.medestrand}'

rule consolidate_cohort_MeDEStrand:
    input:
        lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/MeDEStrand/{sample}.medestrand_output.tsv',
            sample = get_all_samples_list(wildcards.cohort)
        )
    output:
        path_to_data + '/{cohort}/results/MeDEStrand/consolidated/MeDEStrand_AbsMethyl.tsv'
    params:
        data = 'MeDEStrand',
        in_path = path_to_data + '/{cohort}/results/MeDEStrand/',
        out_path = path_to_data + '/{cohort}/results/MeDEStrand/consolidated/'
    resources: cpus=1, mem_mb=30000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/consolidate_cohort_samples.R -i {params.in_path} -o {params.out_path} -d {params.data}'


rule run_QSEA:
    input:
        path_to_data + '/samples/{sample}/merged/bwa_mem/aligned.sorted.markdup.bam'
    output:
        qsea_count = path_to_data + '/samples/{sample}/merged/QSEA/by_chromosome/qsea_count_{chrom}_output.tsv',
        qsea_beta = path_to_data + '/samples/{sample}/merged/QSEA/by_chromosome/qsea_beta_{chrom}_output.tsv',
        qsea_qc = path_to_data + '/samples/{sample}/merged/QSEA/by_chromosome/QCStats_{chrom}_matrix.tsv'
    params:
        out = path_to_data + '/samples/{sample}/merged/QSEA/by_chromosome/'
    resources: cpus=4, mem_mb=30000, time_min='3-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/run_QSEA.R -s {wildcards.sample} -c {wildcards.chrom} -b {input} -o {params.out} --count {output.qsea_count} --beta {output.qsea_beta} --qc {output.qsea_qc}'

rule merge_QSEA_count:
    input:
        lambda wildcards: [path_to_data + '/samples/{{sample}}/merged/QSEA/by_chromosome/qsea_count_{chrom}_output.tsv'.format(chrom=a) for a in chromosomes[wildcards.cohort]]
    output:
        path_to_data + '/samples/{sample}/merged/QSEA/qsea_count_output.tsv'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    run:
        for i, input_file in enumerate(input):
            input_data = pd.read_csv(input_file, delimiter='\t', comment='#')
            if i == 0:
                input_data.to_csv(output[0], header=True, sep='\t', index=False)
            else:
                input_data.to_csv(output[0], header=False, sep='\t', index=False, mode='a')

rule merge_QSEA_beta:
    input:
        lambda wildcards: [path_to_data + '/samples/{{sample}}/merged/QSEA/by_chromosome/qsea_beta_{chrom}_output.tsv'.format(chrom=a) for a in chromosomes[wildcards.cohort]]
    output:
        path_to_data + '/samples/{sample}/merged/QSEA/qsea_beta_output.tsv'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    run:
        for i, input_file in enumerate(input):
            input_data = pd.read_csv(input_file, delimiter='\t', comment='#')
            if i == 0:
                input_data.to_csv(output[0], header=True, sep='\t', index=False)
            else:
                input_data.to_csv(output[0], header=False, sep='\t', index=False, mode='a')

rule merge_QSEA_QCStats:
    input:
        lambda wildcards: [path_to_data + '/samples/{{sample}}/merged/QSEA/by_chromosome/QCStats_{chrom}_matrix.tsv'.format(chrom=a) for a in chromosomes[wildcards.cohort]]
    output:
        path_to_data + '/samples/{sample}/merged/QSEA/QCStats_matrix.tsv'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    run:
        for i, input_file in enumerate(input):
            input_data = pd.read_csv(input_file, delimiter='\t', comment='#')
            if i == 0:
                input_data.to_csv(output[0], header=True, sep='\t', index=False)
            else:
                input_data.to_csv(output[0], header=False, sep='\t', index=False, mode='a')

