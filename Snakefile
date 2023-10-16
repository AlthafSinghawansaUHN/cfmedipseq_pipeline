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
           [path_to_data + '/{cohort}/results/bam_dedup/{sample}.aligned.sorted.dedup.bam.bai'.format(
                    cohort=v[0],
                    sample=v[1]
                ) for v in get_all_samples_with_cohorts()
            ]
        ),
        # Run Consolidated QC
        expand(
            path_to_data + '/{cohort}/results/qc/flagstats_bam/consolidated/flag_stats.tsv',
            cohort = get_active_cohorts()
        ),
        # Run Consolidated MEDIPS
        expand(
            path_to_data + '/{cohort}/results/MEDIPS/consolidated/MEDIPS_Counts.tsv',
            cohort = get_active_cohorts()
        ),
        expand(
            path_to_data + '/{cohort}/results/MEDIPS/consolidated/MEDIPS_QCStats.tsv',
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
        ## NICK pipeline
        # Create aligned bam files 
        #expand(
        #    [path_to_data + '/{cohort}/results/bowtie2_bam_dedup/{sample}.aligned.sorted.markdup.bam.bai'.format(
        #            cohort=v[0],
        #            sample=v[1]
        #        ) for v in get_all_samples_with_cohorts()
        #    ]
        #),
        # Run Consolidated QC
        #expand(
        #    path_to_data + '/{cohort}/results/qc/bowtie2_flagstats/consolidated/flag_stats.tsv',
        #    cohort = get_active_cohorts()
        #),
        # Run Consolidated MEDIPS
        #expand(
        #    path_to_data + '/{cohort}/results/MEDIPS_Nick/consolidated/MEDIPS_Counts.tsv',
        #    cohort = get_active_cohorts()
        #),
        #expand(
        #   path_to_data + '/{cohort}/results/MEDIPS_Nick/consolidated/MEDIPS_QCStats.tsv',
        #    cohort = get_active_cohorts()
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
        'conda_env/trimgalore_installed',
        'conda_env/umi_tools_installed',
        'conda_env/trimmomatic_installed',
        'conda_env/bowtie2_installed'
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

rule umi_tools:
    output:
        'conda_env/umi_tools_installed'
    conda: 'conda_env/umi_tools.yml'
    shell:
        'touch {output}'

rule trimmomatic:
    output:
        'conda_env/trimmomatic_installed'
    conda: 'conda_env/trimmomatic.yml'
    shell:
        'touch {output}'

rule bowtie2:
    output:
        'conda_env/bowtie2_installed'
    conda: 'conda_env/bowtie2.yml'
    shell:
        'touch {output}'

# -------------------------- #
#  Pre-process input FASTQs  #
# -------------------------- #

rule move_fastq:
    input:
        lambda wildcards: get_fastq_path(wildcards.cohort, wildcards.sample, int(wildcards.lib), int(wildcards.read))
    output:
        temp(path_to_data + '/{cohort}/tmp/fastq/{sample}_lib{lib}_R{read}.fastq.gz')
    resources: cpus=1, mem_mb=16000, time_min='24:00:00'
    shell:
        'cp {input} {output}'

# QC of input FASTQs using FASTQC
rule fastqc_fastq:
    input:
        path_to_data + '/{cohort}/tmp/fastq/{sample}_lib{lib}_R{read}.fastq.gz'
    output:
        html=path_to_data + '/{cohort}/results/qc_input/{sample}_lib{lib}_R{read}_fastqc.html',
        zipfile=path_to_data + '/{cohort}/results/qc_input/{sample}_lib{lib}_R{read}_fastqc.zip'
    resources: cpus=1, mem_mb=16000, time_min='24:00:00'
    params:
        outdir = lambda wildcards, output: '/'.join(output.html.split('/')[0:-1])
    conda: 'conda_env/fastqc.yml'
    shell:
        'fastqc --outdir {params.outdir} {input}'

# Extract Barcodes using UMI_tools
# Using umi_tools extract function to extract barcodes, requires a barcode pattern
# will set up a barcode pattern (separate setting from the bpattern supplied to extract_barcodes)
# using cohort settings config > data > defaults > umi_tools or config > [cohort] > settings > umi_tools
    ###umi_tools extract --extract-method=regex \
    ###                  --bc-pattern="(?P<umi_1>^[ACGT]{3}[ACG])(?P<discard_1>T)|(?P<umi_1>^[ACGT]{3})(?P<discard_1>T)" \
    ###                  --bc-pattern2="(?P<umi_2>^[ACGT]{3}[ACG])(?P<discard_2>T)|(?P<umi_2>^[ACGT]{3})(?P<discard_2>T)" \
    ###                  --stdin=cfMe_PDAC_PCSI_1010_Ct_T_1000_R1.fastq.gz  \
    ###                  --read2-in=cfMe_PDAC_PCSI_1010_Ct_T_1000_R2.fastq.gz  \
    ###                  --stdout=cfMe_PDAC_PCSI_1010_Ct_T_1000_R1_barcoded.fastq.gz \
    ###                  --read2-out=cfMe_PDAC_PCSI_1010_Ct_T_1000_R2_barcoded.fastq.gz \
    ###                  --log=log.txt
rule umi_tools_extract_pe:
    input:
        R1 =    path_to_data + '/{cohort}/tmp/fastq/{sample}_lib{lib}_R1.fastq.gz',
        R2 =    path_to_data + '/{cohort}/tmp/fastq/{sample}_lib{lib}_R2.fastq.gz',
        R1_qc =    path_to_data + '/{cohort}/results/qc_input/{sample}_lib{lib}_R1_fastqc.html',
        R2_qc =    path_to_data + '/{cohort}/results/qc_input/{sample}_lib{lib}_R2_fastqc.html'
    output:
        R1 =    temp(path_to_data + '/{cohort}/tmp/umi_tools_extract_pe/{sample}_lib{lib}/{sample}_lib{lib}_umi_extract_R1.fastq.gz'),
        R2 =    temp(path_to_data + '/{cohort}/tmp/umi_tools_extract_pe/{sample}_lib{lib}/{sample}_lib{lib}_umi_extract_R2.fastq.gz'),
        logs =   path_to_data + '/{cohort}/results/qc_input_umi_tools_extract_pe/{sample}_lib{lib}_barcode_logs.txt'
    params:
        umi_tools_settings = lambda wildcards: get_cohort_config(wildcards.cohort)['umi_tools']
    resources: cpus=1, mem_mb=30000, time_min='1-00:00:00'
    conda: 'conda_env/umi_tools.yml'
    shell:
        'umi_tools extract {params.umi_tools_settings} -I {input.R1} --read2-in={input.R2} --log={output.logs} --stdout={output.R1} --read2-out={output.R2}'

rule umi_tools_extract_se:
    input:
        R1 =	path_to_data + '/{cohort}/tmp/fastq/{sample}_lib{lib}_R1.fastq.gz',
        R1_qc =    path_to_data + '/{cohort}/results/qc_input/{sample}_lib{lib}_R1_fastqc.html'
    output:
        R1 =    temp(path_to_data + '/{cohort}/tmp/umi_tools_extract_se/{sample}_lib{lib}/{sample}_lib{lib}_umi_extract_R1.fastq.gz'),
        logs =   path_to_data + '/{cohort}/results/qc_input_umi_tools_extract_se/{sample}_lib{lib}_barcode_logs.txt'
    params:
        umi_tools_settings = lambda wildcards: get_cohort_config(wildcards.cohort)['umi_tools']
    resources: cpus=1, mem_mb=30000, time_min='1-00:00:00'
    conda: 'conda_env/umi_tools.yml'
    shell:
        'umi_tools extract {params.umi_tools_settings} -I {input.R1} --log={output.logs} --stdout={output.R1}'

# Trims FASTQ using trimgalore to remove barcode sequences
# By default, trims 10 base pairs from the 5' end, which seems to be correct for OICR cfMeDIP-seq output.
# This can be configured in the config.yml under data > cohorts > settings > trimgalore.
rule trim_fastq_umi_tools_extract_pe:
    input:
        R1 = path_to_data + '/{cohort}/tmp/umi_tools_extract_pe/{sample}_lib{lib}/{sample}_lib{lib}_umi_extract_R1.fastq.gz',
        R2 = path_to_data + '/{cohort}/tmp/umi_tools_extract_pe/{sample}_lib{lib}/{sample}_lib{lib}_umi_extract_R2.fastq.gz'
    output:
        trimmed_1 = temp(path_to_data + '/{cohort}/tmp/trim_fastq_umi_tools_extract_pe/{sample}_lib{lib}_umi_extract_R1_val_1.fq.gz'),
        trimmed_2 = temp(path_to_data + '/{cohort}/tmp/trim_fastq_umi_tools_extract_pe/{sample}_lib{lib}_umi_extract_R2_val_2.fq.gz'),
        report_1 = path_to_data + '/{cohort}/results/qc/trim_report/{sample}_lib{lib}_umi_extract_R1.fastq.gz_trimming_report.txt',
        report_2 = path_to_data + '/{cohort}/results/qc/trim_report/{sample}_lib{lib}_umi_extract_R2.fastq.gz_trimming_report.txt'
    params:
        outdir = lambda wildcards, output: '/'.join(output.trimmed_1.split('/')[0:-1]),
        trimgalore_settings = lambda wildcards: get_cohort_config(wildcards.cohort)['trimgalore']
    resources: cpus=12, mem_mb=30000, time_min='24:00:00'
    conda: 'conda_env/trimgalore.yml'
    shell:
        'trim_galore --cores 12 --paired {params.trimgalore_settings} --output_dir {params.outdir} {input.R1} {input.R2} && mv ' + path_to_data + '/{wildcards.cohort}/tmp/trim_fastq_umi_tools_extract_pe/{wildcards.sample}_lib{wildcards.lib}_umi_extract_R*.fastq.gz_trimming_report.txt ' + path_to_data + '/{wildcards.cohort}/results/qc/trim_report/'

rule trim_fastq_umi_tools_extract_se:
    input:
        R1 = path_to_data + '/{cohort}/tmp/umi_tools_extract_se/{sample}_lib{lib}/{sample}_lib{lib}_umi_extract_R1.fastq.gz'
    output:
        trimmed_1 = temp(path_to_data + '/{cohort}/tmp/trim_fastq_umi_tools_extract_se/{sample}_lib{lib}_umi_extract_R1_trimmed.fq.gz'),
        report_1 = path_to_data + '/{cohort}/results/qc/trim_report/{sample}_lib{lib}_umi_extract_R1.fastq.gz_trimming_report.txt',
    params:
        outdir = lambda wildcards, output: '/'.join(output.trimmed_1.split('/')[0:-1]),
        trimgalore_settings = lambda wildcards: get_cohort_config(wildcards.cohort)['trimgalore']
    resources: cpus=12, mem_mb=30000, time_min='24:00:00'
    conda: 'conda_env/trimgalore.yml'
    shell:
        'trim_galore --cores 12 {params.trimgalore_settings} --output_dir {params.outdir} {input.R1} && mv ' + path_to_data + '/{wildcards.cohort}/tmp/trim_fastq_umi_tools_extract_se/{wildcards.sample}_lib{wildcards.lib}_umi_extract_R*.fastq.gz_trimming_report.txt ' + path_to_data + '/{wildcards.cohort}/results/qc/trim_report/'

rule trim_fastq_pe:
    input:
        R1 = lambda wildcards: get_fastq_path(wildcards.cohort, wildcards.sample, int(wildcards.lib), int(1)),
        R2 = lambda wildcards: get_fastq_path(wildcards.cohort, wildcards.sample, int(wildcards.lib), int(2))
    output:
        trimmed_1 = temp(path_to_data + '/{cohort}/tmp/trim_fastq_pe/{sample}_lib{lib}_R1_val_1.fq.gz'),
        report_1 = path_to_data + '/{cohort}/results/qc/trim_report/{sample}_lib{lib}_R1.fastq.gz_trimming_report.txt',
        trimmed_2 = temp(path_to_data + '/{cohort}/tmp/trim_fastq_pe/{sample}_lib{lib}_R2_val_2.fq.gz'),
        report_2 = path_to_data + '/{cohort}/results/qc/trim_report/{sample}_lib{lib}_R2.fastq.gz_trimming_report.txt',
    params:
        outdir = lambda wildcards, output: '/'.join(output.trimmed_1.split('/')[0:-1]),
        trimgalore_settings = lambda wildcards: get_cohort_config(wildcards.cohort)['trimgalore']
    resources: cpus=12, mem_mb=30000, time_min='24:00:00'
    conda: 'conda_env/trimgalore.yml'
    shell:
        'trim_galore --cores 12 --paired {params.trimgalore_settings} --output_dir {params.outdir} {input.R1} {input.R2} && mv ' + path_to_data + '/{wildcards.cohort}/tmp/trim_fastq_pe/{wildcards.sample}_lib{wildcards.lib}_umi_extract_R*.fastq.gz_trimming_report.txt ' + path_to_data + '/{wildcards.cohort}/results/qc/trim_report/'

rule trim_fastq_se:
    input:
        R1 = lambda wildcards: get_fastq_path(wildcards.cohort, wildcards.sample, int(wildcards.lib), int(1)),
        R1_qc =    path_to_data + '/{cohort}/results/qc_input/{sample}_lib{lib}_R1_fastqc.html',
    output:
        trimmed_1 = temp(path_to_data + '/{cohort}/tmp/trim_fastq_se/{sample}_lib{lib}_R1_trimmed.fq.gz'),
        report_1 = path_to_data + '/{cohort}/results/qc/trim_report/{sample}_lib{lib}_R1.fastq.gz_trimming_report.txt',
    params:
        outdir = lambda wildcards, output: '/'.join(output.trimmed_1.split('/')[0:-1]),
        trimgalore_settings = lambda wildcards: get_cohort_config(wildcards.cohort)['trimgalore']
    resources: cpus=12, mem_mb=30000, time_min='24:00:00'
    conda: 'conda_env/trimgalore.yml'
    shell:
        'trim_galore --cores 12 {params.trimgalore_settings} --output_dir {params.outdir} {input.R1}  && mv ' + path_to_data + '/{wildcards.cohort}/tmp/trim_fastq_se/{wildcards.sample}_lib{wildcards.lib}_R1.fastq.gz_trimming_report.txt ' + path_to_data + '/{wildcards.cohort}/results/qc/trim_report/'

# Extract Barcodes using ConsensusCruncher
# Pulls the path to extract_barcodes.py from config > paths > dependencies > extract_barcodes_path
# Runs extract_barcodes.py based on cohort settings config > data > defaults or
# config > [cohort] > settings
rule extract_barcodes:
    input:
        R1 =    lambda wildcards: get_fastq_path(wildcards.cohort, wildcards.sample, int(wildcards.lib), int(1)),
        R2 =    lambda wildcards: get_fastq_path(wildcards.cohort, wildcards.sample, int(wildcards.lib), int(2)),
        R1_qc =    path_to_data + '/{cohort}/results/qc_input/{sample}_lib{lib}_R1_fastqc.html',
        R2_qc =    path_to_data + '/{cohort}/results/qc_input/{sample}_lib{lib}_R2_fastqc.html'
    output:
        R1 =    temp(path_to_data + '/{cohort}/tmp/extract_barcodes/{sample}_lib{lib}/{sample}_lib{lib}_extract_barcode_R1.fastq.gz'),
        R2 =    temp(path_to_data + '/{cohort}/tmp/extract_barcodes/{sample}_lib{lib}/{sample}_lib{lib}_extract_barcode_R2.fastq.gz'),
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

rule trim_fastq_extract_barcodes_pe:
    input:
        R1 = path_to_data + '/{cohort}/tmp/extract_barcodes/{sample}_lib{lib}/{sample}_lib{lib}_extract_barcode_R1.fastq.gz',
        R2 = path_to_data + '/{cohort}/tmp/extract_barcodes/{sample}_lib{lib}/{sample}_lib{lib}_extract_barcode_R2.fastq.gz'
    output:
        trimmed_1 = temp(path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_extract_barcode_R1_val_1.fastq.gz'),
        trimmed_2 = temp(path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_extract_barcode_R2_val_2.fastq.gz'),
        report_1 = path_to_data + '/{cohort}/results/qc/trim_report/{sample}_lib{lib}_extract_barcode_R1.fastq_trimming_report.txt',
        report_2 = path_to_data + '/{cohort}/results/qc/trim_report/{sample}_lib{lib}_extract_barcode_R2.fastq_trimming_report.txt'
    params:
        outdir = lambda wildcards, output: '/'.join(output.trimmed_1.split('/')[0:-1]),
        trimgalore_settings = lambda wildcards: get_cohort_config(wildcards.cohort)['trimgalore']
    resources: cpus=12, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/trimgalore.yml'
    shell:
        'trim_galore --cores 12 --dont_gzip --paired {params.trimgalore_settings} --output_dir {params.outdir} {input.R1} {input.R2} && cp ' + path_to_data + '/{wildcards.cohort}/tmp/trim_fastq/{wildcards.sample}_lib{wildcards.lib}_extract_barcode_R*.fastq.gz_trimming_report.txt ' + path_to_data + '/{wildcards.cohort}/results/qc/trim_report/'


# Trims FASTQ using trimmomatics to remove barcode sequences
# By default, trims TRAILING:3 MINLEN:30 SLIDINGWINDOW:4:7
# This can be configured in the config.yml under data > cohorts > settings > trimmomatic.
rule trimmomatic_fastq_paired:
    input:
        R1 = path_to_data + '/{cohort}/tmp/umi_tools_extract/{sample}_lib{lib}/{sample}_lib{lib}_umi_extract_R1.fastq.gz',
        R2 = path_to_data + '/{cohort}/tmp/umi_tools_extract/{sample}_lib{lib}/{sample}_lib{lib}_umi_extract_R2.fastq.gz'
    output:
        trimmed_1 = temp(path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_umi_extract_trimmed_R1.fastq.gz'),
        trimmed_2 = temp(path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_umi_extract_trimmed_R2.fastq.gz'),
        trimmed_unpaired_1 = temp(path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_umi_extract_unpaired_trimmed_R1.fastq.gz'),
        trimmed_unpaired_2 = temp(path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_umi_extract_unpaired_trimmed_R2.fastq.gz'),
        report = path_to_data + '/{cohort}/results/qc/trim_report/{sample}_lib{lib}_umi_extract.fastq_trimmomatic_trimming_report.txt'
    params:
        outdir = lambda wildcards, output: '/'.join(output.trimmed_1.split('/')[0:-1]),
        trimgalore_settings = lambda wildcards: get_cohort_config(wildcards.cohort)['trimmomatic']
    resources: cpus=4, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/trimmomatic.yml'
    shell:
        'trimmomatic PE -threads 4 -trimlog {output.report} {input.R1} {input.R2} {output.trimmed_1} {output.trimmed_unpaired_1} {output.trimmed_2} {output.trimmed_unpaired_2} TRAILING:3 MINLEN:30 SLIDINGWINDOW:4:7'


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

#Run BWA mem on FASTQs after extracting barcodes and trimming for paired end data
rule bwa_mem_umi_tools_extract_pe:
    input:
        R1 = path_to_data + '/{cohort}/tmp/trim_fastq_umi_tools_extract_pe/{sample}_lib{lib}_umi_extract_R1_val_1.fq.gz',
        R2 = path_to_data + '/{cohort}/tmp/trim_fastq_umi_tools_extract_pe/{sample}_lib{lib}_umi_extract_R2_val_2.fq.gz',
    output:
        temp(path_to_data + '/{cohort}/tmp/bwa_mem_umi_tools_extract_pe/{sample}_lib{lib}.bam')
    resources: partition='himem', cpus=12, mem_mb=60000, time_min='72:00:00'
    params:
        read_group = lambda wildcards, input: get_read_group_from_fastq(
            fastq_file = get_fastq_path(wildcards.cohort, wildcards.sample, wildcards.lib),
            sample_name = wildcards.sample
        ),
        bwa_index = lambda wildcards: get_cohort_config(wildcards.cohort)['bwa_index']
    conda: 'conda_env/samtools.yml'
    shell:
        "bwa mem -M -t12 -R'{params.read_group}' {params.bwa_index} {input.R1} {input.R2} | samtools view -bSo {output}"

rule bwa_mem_umi_tools_extract_se:
    input:
        R1 = path_to_data + '/{cohort}/tmp/trim_fastq_umi_tools_extract_se/{sample}_lib{lib}_umi_extract_R1_trimmed.fq.gz'
    output:
        temp(path_to_data + '/{cohort}/tmp/bwa_mem_umi_tools_extract_se/{sample}_lib{lib}.bam')
    resources: partition='himem', cpus=12, mem_mb=60000, time_min='72:00:00'
    params:
        read_group = lambda wildcards, input: get_read_group_from_fastq(
            fastq_file = get_fastq_path(wildcards.cohort, wildcards.sample, wildcards.lib),
            sample_name = wildcards.sample
        ),
        bwa_index = lambda wildcards: get_cohort_config(wildcards.cohort)['bwa_index']
    conda: 'conda_env/samtools.yml'
    shell:
        "bwa mem -M -t12 -R'{params.read_group}' {params.bwa_index} {input.R1} | samtools view -bSo {output}"

#Run BWA mem on FASTQs of no umi extraction data
rule bwa_mem_pe:
    input:
        R1 = path_to_data + '/{cohort}/tmp/trim_fastq_pe/{sample}_lib{lib}_R1_val_1.fq.gz',
        R2 = path_to_data + '/{cohort}/tmp/trim_fastq_pe/{sample}_lib{lib}_R2_val_2.fq.gz'
    output:
        temp(path_to_data + '/{cohort}/tmp/bwa_mem_pe/{sample}_lib{lib}.bam')
    resources: partition='himem', cpus=12, mem_mb=60000, time_min='72:00:00'
    params:
        read_group = lambda wildcards, input: get_read_group_from_fastq(
            fastq_file = get_fastq_path(wildcards.cohort, wildcards.sample, wildcards.lib),
            sample_name = wildcards.sample
        ),
        bwa_index = lambda wildcards: get_cohort_config(wildcards.cohort)['bwa_index']
    conda: 'conda_env/samtools.yml'
    shell:
        "bwa mem -M -t12 -R'{params.read_group}' {params.bwa_index} {input.R1} {input.R2} | samtools view -bSo {output}"

rule bwa_mem_se:
    input:
        R1 = path_to_data + '/{cohort}/tmp/trim_fastq_se/{sample}_lib{lib}_R1_trimmed.fq.gz'
    output:
        temp(path_to_data + '/{cohort}/tmp/bwa_mem_se/{sample}_lib{lib}.bam')
    resources: partition='himem', cpus=12, mem_mb=60000, time_min='72:00:00'
    params:
       read_group = lambda wildcards, input: get_read_group_from_fastq(
          fastq_file = get_fastq_path(wildcards.cohort, wildcards.sample, wildcards.lib),
            sample_name = wildcards.sample
        ),
        bwa_index = lambda wildcards: get_cohort_config(wildcards.cohort)['bwa_index']
    conda: 'conda_env/samtools.yml'
    shell:
        "bwa mem -M -t12 -R'{params.read_group}' {params.bwa_index} {input.R1} | samtools view -bSo {output}"			

# Sort Bam file and remove unmapped reads
rule bam_to_sorted_bam_umi_tools_extract_pe:
    input:
        path_to_data + '/{cohort}/tmp/bwa_mem_umi_tools_extract_pe/{sample}_lib{lib}.bam'
    output:
        bam = temp(path_to_data + '/{cohort}/tmp/bwa_mem_umi_tools_extract_pe/{sample}_lib{lib}.sorted.bam'),
        index = temp(path_to_data + '/{cohort}/tmp/bwa_mem_umi_tools_extract_pe/{sample}_lib{lib}.sorted.bam.bai'),
    resources: cpus=32, mem_mb=30000, time_min='72:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        clean(r'''
        samtools view -buS -f 2 -F 2828 -@32 {input} |
        samtools fixmate -m - - |
        samtools sort -@32 -o {output.bam} && samtools index {output.bam} 
        ''')

rule bam_to_sorted_bam_umi_tools_extract_se:
    input:
        path_to_data + '/{cohort}/tmp/bwa_mem_umi_tools_extract_se/{sample}_lib{lib}.bam'
    output:
        bam = temp(path_to_data + '/{cohort}/tmp/bwa_mem_umi_tools_extract_se/{sample}_lib{lib}.sorted.bam'),
        index = temp(path_to_data + '/{cohort}/tmp/bwa_mem_umi_tools_extract_se/{sample}_lib{lib}.sorted.bam.bai'),
    resources: cpus=32, mem_mb=30000, time_min='72:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        clean(r'''
        samtools view -buS -F 2820 -@32 {input} |
        samtools fixmate -m - - |
        samtools sort -@32 -o {output.bam} && samtools index {output.bam} 
        ''')

rule bam_to_sorted_bam_pe:
    input:
        path_to_data + '/{cohort}/tmp/bwa_mem_pe/{sample}_lib{lib}.bam'
    output:
        bam = temp(path_to_data + '/{cohort}/tmp/bwa_mem_pe/{sample}_lib{lib}.sorted.bam'),
        index = temp(path_to_data + '/{cohort}/tmp/bwa_mem_pe/{sample}_lib{lib}.sorted.bam.bai'),
    resources: cpus=32, mem_mb=30000, time_min='72:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        clean(r'''
        samtools view -buS -f 2 -F 2828 -@32 {input} |
        samtools fixmate -m - - |
        samtools sort -@32 -o {output.bam} && samtools index {output.bam} 
        ''')

rule bam_to_sorted_bam_se:
    input:
        path_to_data + '/{cohort}/tmp/bwa_mem_se/{sample}_lib{lib}.bam'
    output:
        bam = temp(path_to_data + '/{cohort}/tmp/bwa_mem_se/{sample}_lib{lib}.sorted.bam'),
        index = temp(path_to_data + '/{cohort}/tmp/bwa_mem_se/{sample}_lib{lib}.sorted.bam.bai'),
    resources: cpus=32, mem_mb=30000, time_min='72:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        clean(r'''
        samtools view -buS -F 2820 -@32 {input} |
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
rule merge_bam_umi_tools_extract_pe:
    input:
        lambda wildcards: expand(
                path_to_data + '/' + wildcards.cohort + '/tmp/bwa_mem_umi_tools_extract_pe/' + wildcards.sample + '_lib{lib}.sorted.bam',
                lib=get_libraries_of_sample(wildcards.sample)
        )
    output:
        bam = temp(path_to_data + '/{cohort}/tmp/merge_bam_umi_tools_extract_pe/{sample}.aligned.sorted.bam'),
        index = temp(path_to_data + '/{cohort}/tmp/merge_bam_umi_tools_extract_pe/{sample}.aligned.sorted.bam.bai')
    resources: cpus=1, mem_mb=30000, time_min='24:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        'samtools merge {output.bam} {input} && samtools index {output.bam}'

rule merge_bam_umi_tools_extract_se:
    input:
        lambda wildcards: expand(
                path_to_data + '/' + wildcards.cohort + '/tmp/bwa_mem_umi_tools_extract_se/' + wildcards.sample + '_lib{lib}.sorted.bam',
                lib=get_libraries_of_sample(wildcards.sample)
        )
    output:
        bam = temp(path_to_data + '/{cohort}/tmp/merge_bam_umi_tools_extract_se/{sample}.aligned.sorted.bam'),
        index = temp(path_to_data + '/{cohort}/tmp/merge_bam_umi_tools_extract_se/{sample}.aligned.sorted.bam.bai')
    resources: cpus=1, mem_mb=30000, time_min='24:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        'samtools merge {output.bam} {input} && samtools index {output.bam}'

rule merge_bam_pe:
    input:
        lambda wildcards: expand(
                path_to_data + '/' + wildcards.cohort + '/tmp/bwa_mem_pe/' + wildcards.sample + '_lib{lib}.sorted.bam',
                lib=get_libraries_of_sample(wildcards.sample)
        )
    output:
        bam = temp(path_to_data + '/{cohort}/tmp/merge_bam_pe/{sample}.aligned.sorted.bam'),
        index = temp(path_to_data + '/{cohort}/tmp/merge_bam_pe/{sample}.aligned.sorted.bam.bai')
    resources: cpus=1, mem_mb=30000, time_min='24:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        'samtools merge {output.bam} {input} && samtools index {output.bam}'

rule merge_bam_se:
    input:
        lambda wildcards: expand(
                path_to_data + '/' + wildcards.cohort + '/tmp/bwa_mem_se/' + wildcards.sample + '_lib{lib}.sorted.bam',
                lib=get_libraries_of_sample(wildcards.sample)
        )
    output:
        bam = temp(path_to_data + '/{cohort}/tmp/merge_bam_se/{sample}.aligned.sorted.bam'),
        index = temp(path_to_data + '/{cohort}/tmp/merge_bam_se/{sample}.aligned.sorted.bam.bai')
    resources: cpus=1, mem_mb=30000, time_min='24:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        'samtools merge {output.bam} {input} && samtools index {output.bam}'

# Deduplication
def get_merged_bam(wildcards):
    paired_data = get_cohort_config(wildcards.cohort)['paired']
    umi_data = get_cohort_config(wildcards.cohort)['UMI']
    if paired_data and umi_data:
        return path_to_data + '/' + wildcards.cohort + '/tmp/merge_bam_umi_tools_extract_pe/' + wildcards.sample + '.aligned.sorted.bam'
    elif paired_data == False and umi_data:
        return path_to_data + '/' + wildcards.cohort + '/tmp/merge_bam_umi_tools_extract_se/' + wildcards.sample + '.aligned.sorted.bam'
    elif paired_data and umi_data == False:
        return path_to_data + '/' + wildcards.cohort + '/tmp/merge_bam_pe/' + wildcards.sample + '.aligned.sorted.bam'
    else:
        return path_to_data + '/' + wildcards.cohort + '/tmp/merge_bam_se/' + wildcards.sample + '.aligned.sorted.bam'

def get_merged_index(wildcards):
    paired_data = get_cohort_config(wildcards.cohort)['paired']
    umi_data = get_cohort_config(wildcards.cohort)['UMI']
    if paired_data and umi_data:
        return path_to_data + '/' + wildcards.cohort + '/tmp/merge_bam_umi_tools_extract_pe/' + wildcards.sample + '.aligned.sorted.bam.bai'
    elif paired_data == False and umi_data:
        return path_to_data + '/' + wildcards.cohort + '/tmp/merge_bam_umi_tools_extract_se/' + wildcards.sample + '.aligned.sorted.bam.bai'
    elif paired_data and umi_data == False:
        return path_to_data + '/' + wildcards.cohort + '/tmp/merge_bam_pe/' + wildcards.sample + '.aligned.sorted.bam.bai'
    else:
        return path_to_data + '/' + wildcards.cohort + '/tmp/merge_bam_se/' + wildcards.sample + '.aligned.sorted.bam.bai'

# Bam markdup and create index. 
rule bam_markdup_pe:
    input:
        bam = lambda wildcards: get_merged_bam(wildcards),
        index = lambda wildcards: get_merged_index(wildcards)
    output:
        bam = temp(path_to_data + '/{cohort}/tmp/bam_markdup_pe/{sample}.aligned.sorted.dedup.bam')
    resources: partition='himem', cpus=1, mem_mb=61440, time_min='24:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        "samtools markdup -r {input.bam} {output.bam}"

rule bam_markdup_se:
    input:
        bam = lambda wildcards: get_merged_bam(wildcards),
        index = lambda wildcards: get_merged_index(wildcards)
    output:
        bam = temp(path_to_data + '/{cohort}/tmp/bam_markdup_se/{sample}.aligned.sorted.dedup.bam')
    resources: partition='himem', cpus=1, mem_mb=61440, time_min='24:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        "samtools markdup -r {input.bam} {output.bam}"

# Umi_tools dedup. 
rule umi_tools_dedup_pe:
    input:
        bam = lambda wildcards: get_merged_bam(wildcards),
        index = lambda wildcards: get_merged_index(wildcards)
    output:
        bam = temp(path_to_data + '/{cohort}/tmp/umi_tools_dedup_pe/{sample}.aligned.sorted.dedup.bam'),
        stats_per_umi = path_to_data + '/{cohort}/results/umi_tools_dedup_pe/stats/{sample}_per_umi.tsv',
        stats_per_umi_per_position = path_to_data + '/{cohort}/results/umi_tools_dedup_pe/stats/{sample}_per_umi_per_position.tsv',
        stats_edit_distance = path_to_data + '/{cohort}/results/umi_tools_dedup_pe/stats/{sample}_edit_distance.tsv'
    params:
        stats = path_to_data + '/{cohort}/results/umi_tools_dedup_pe/stats/{sample}'
    resources: partition='himem', cpus=1, mem_mb=61440, time_min='24:00:00'
    conda: 'conda_env/umi_tools.yml'
    shell:
        "umi_tools dedup --paired --umi-separator='_' -I {input.bam} --output-stats={params.stats} -S {output.bam}"

rule umi_tools_dedup_se:
    input:
        bam = lambda wildcards: get_merged_bam(wildcards),
        index = lambda wildcards: get_merged_index(wildcards)
    output:
        bam = temp(path_to_data + '/{cohort}/tmp/umi_tools_dedup_se/{sample}.aligned.sorted.dedup.bam'),
        stats_per_umi = path_to_data + '/{cohort}/results/umi_tools_dedup_se/stats/{sample}_per_umi.tsv',
        stats_per_umi_per_position = path_to_data + '/{cohort}/results/umi_tools_dedup_se/stats/{sample}_per_umi_per_position.tsv',
        stats_edit_distance = path_to_data + '/{cohort}/results/umi_tools_dedup_se/stats/{sample}_edit_distance.tsv'
    params:
        stats = path_to_data + '/{cohort}/results/umi_tools_dedup_se/stats/{sample}'
    resources: partition='himem', cpus=1, mem_mb=61440, time_min='24:00:00'
    conda: 'conda_env/umi_tools.yml'
    shell:
        "umi_tools dedup --umi-separator='_' -I {input.bam} --output-stats={params.stats} -S {output.bam}"

def get_final_bam(wildcards):
    paired_data = get_cohort_config(wildcards.cohort)['paired']
    umi_data = get_cohort_config(wildcards.cohort)['UMI']
    if paired_data and umi_data:
        return path_to_data + '/' + wildcards.cohort + '/tmp/umi_tools_dedup_pe/' + wildcards.sample + '.aligned.sorted.dedup.bam'
    elif paired_data == False and umi_data:
        return path_to_data + '/' + wildcards.cohort + '/{cohort}/tmp/umi_tools_dedup_se/' + wildcards.sample + '.aligned.sorted.dedup.bam'
    elif paired_data and umi_data == False:
        return path_to_data + '/' + wildcards.cohort + '/tmp/bam_markdup_pe/' + wildcards.sample + '.aligned.sorted.dedup.bam'
    else:
        return path_to_data + '/' + wildcards.cohort + '/tmp/bam_markdup_se/' + wildcards.sample + '.aligned.sorted.dedup.bam'

# Create index. 
# This step finalizes the definitive BAM file.
rule samtools_index:
    input:
        bam = lambda wildcards: get_final_bam(wildcards)
    output:
        bam = path_to_data + '/{cohort}/results/bam_dedup/{sample}.aligned.sorted.dedup.bam',
        index = path_to_data + '/{cohort}/results/bam_dedup/{sample}.aligned.sorted.dedup.bam.bai'
    resources: cpus=1, mem_mb=30000, time_min='24:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        "mv {input.bam} {output.bam} && samtools index {output.bam}"

# Extract Barcodes method
def get_paired_or_single_end_sorted_merged_bam(wildcards):
    paired_data = get_cohort_config(wildcards.cohort)['paired']
    if paired_data:
        return path_to_data + '/' + wildcards.cohort + '/tmp/merge_bam_extract_barcodes/' + wildcards.sample + '.paired.aligned.sorted.bam'
    else:
        return path_to_data + '/' + wildcards.cohort + '/tmp/merge_bam/' + wildcards.sample + '.single.aligned.sorted.bam'

rule bwa_mem_extract_barcodes_pe:
    input:
        R1 = path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_extract_barcode_R1_val_1.fastq.gz',
        R2 = path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_extract_barcode_R2_val_2.fastq.gz',
    output:
        temp(path_to_data + '/{cohort}/tmp/bwa_mem_extract_barcodes/{sample}_lib{lib}.paired.bam')
    resources: cpus=12, mem_mb=30000, time_min='72:00:00'
    params:
        read_group = lambda wildcards, input: get_read_group_from_fastq(
            fastq_file = get_fastq_path(wildcards.cohort, wildcards.sample, wildcards.lib),
            sample_name = wildcards.sample
        ),
        bwa_index = lambda wildcards: get_cohort_config(wildcards.cohort)['bwa_index']
    conda: 'conda_env/samtools.yml'
    shell:
        "bwa mem -M -t12 -R'{params.read_group}' {params.bwa_index} {input.R1} {input.R2} | samtools view -bSo {output}"

rule bam_to_sorted_bam_paired_extract_barcodes:
    input:
        path_to_data + '/{cohort}/tmp/bwa_mem_extract_barcodes/{sample}_lib{lib}.paired.bam'
    output:
        bam = temp(path_to_data + '/{cohort}/tmp/bwa_mem_extract_barcodes/{sample}_lib{lib}.paired.sorted.bam'),
        index = temp(path_to_data + '/{cohort}/tmp/bwa_mem_extract_barcodes/{sample}_lib{lib}.paired.sorted.bam.bai'),
    resources: cpus=32, mem_mb=30000, time_min='72:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        clean(r'''
        samtools view -buS -f 2 -F 2828 -@32 {input} |
        samtools fixmate -m - - |
        samtools sort -@32 -o {output.bam} && samtools index {output.bam} 
        ''')

rule merge_bam_extract_barcodes_pe:
    input:
        lambda wildcards: expand(
                path_to_data + '/' + wildcards.cohort + '/tmp/bwa_mem_extract_barcodes/' + wildcards.sample + '_lib{lib}.paired.sorted.bam',
                lib=get_libraries_of_sample(wildcards.sample)
        )
    output:
        bam = temp(path_to_data + '/{cohort}/tmp/merge_bam_extract_barcodes/{sample}.paired.aligned.sorted.bam'),
        index = temp(path_to_data + '/{cohort}/tmp/merge_bam_extract_barcodes/{sample}.paired.aligned.sorted.bam.bai')
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        'samtools merge {output.bam} {input} && samtools index {output.bam}'

rule bam_markdup_extract_barcodes:
    input:
        path_to_data + '/{cohort}/tmp/merge_bam_extract_barcodes/{sample}.paired.aligned.sorted.bam'
    output:
        bam = path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.dedup.bam',
        index = path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.dedup.bam.bai'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        "samtools markdup -r {input} {output.bam} && samtools index {output.bam}"

#bowtie2 -p $nt --rg-id $ID --rg $SM --rg $PU --rg $PL --rg $LB -x $bt2Index -1 ${sample}_R1.fastq.gz -2 ${sample}_R2.fastq.gz | samtools view -bSo $bamOutDir/${sample}.bam
# Run bowtie2 on FASTQs after extracting barcodes and trimming for paired end data
rule bowtie2_paired:
    input:
        R1 = path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_umi_extract_trimmed_R1.fastq.gz',
        R2 = path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_umi_extract_trimmed_R2.fastq.gz'
    output:
        temp(path_to_data + '/{cohort}/tmp/bowtie2/{sample}_lib{lib}.paired.bam')
    resources: cpus=4, mem_mb=30000, time_min='72:00:00'
    params:
        read_group = lambda wildcards, input: get_read_group_from_fastq(
            fastq_file = get_fastq_path(wildcards.cohort, wildcards.sample, wildcards.lib),
            sample_name = wildcards.sample
        ),
        bt2_index = lambda wildcards: get_cohort_config(wildcards.cohort)['bt2_index']
    conda: 'conda_env/bowtie2.yml'
    shell:
        "bowtie2 -p 4 -R'{params.read_group}' -x {params.bt2_index} -1 {input.R1} -2 {input.R2} | samtools view -bSo {output}"

rule bam_to_sorted_bam_bowtie2_paired:
    input:
        path_to_data + '/{cohort}/tmp/bowtie2/{sample}_lib{lib}.paired.bam'
    output:
        bam = temp(path_to_data + '/{cohort}/tmp/bowtie2/{sample}_lib{lib}.paired.sorted.bam'),
        index = temp(path_to_data + '/{cohort}/tmp/bowtie2/{sample}_lib{lib}.paired.sorted.bam.bai'),
    resources: cpus=32, mem_mb=30000, time_min='72:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        clean(r'''
        samtools view -buS -F 4 -@32 {input} |
        samtools fixmate -m - - |
        samtools sort -@32 -o {output.bam} && samtools index {output.bam} 
        ''')

rule merge_bam_bowtie2_paired:
    input:
        lambda wildcards: expand(
                path_to_data + '/' + wildcards.cohort + '/tmp/bowtie2/' + wildcards.sample + '_lib{lib}.paired.sorted.bam',
                lib=get_libraries_of_sample(wildcards.sample)
        )
    output:
        bam = temp(path_to_data + '/{cohort}/tmp/bowtie2_merge_bam/{sample}.paired.aligned.sorted.bam'),
        index = temp(path_to_data + '/{cohort}/tmp/bowtie2_merge_bam/{sample}.paired.aligned.sorted.bam.bai')
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        'samtools merge {output.bam} {input} && samtools index {output.bam}'

rule umi_tools_dedup_bowtie2:
    input:
        bam = path_to_data + '/{cohort}/tmp/bowtie2_merge_bam/{sample}.paired.aligned.sorted.bam',
        index = path_to_data + '/{cohort}/tmp/bowtie2_merge_bam/{sample}.paired.aligned.sorted.bam.bai'
    output:
        bam = path_to_data + '/{cohort}/results/bowtie2_bam_dedup/{sample}.aligned.sorted.dedup.bam',
        stats_per_umi = path_to_data + '/{cohort}/results/bowtie2_bam_dedup/stats/{sample}_per_umi.tsv',
        stats_per_umi_per_position = path_to_data + '/{cohort}/results/bowtie2_bam_dedup/stats/{sample}_per_umi_per_position.tsv',
        stats_edit_distance = path_to_data + '/{cohort}/results/bowtie2_bam_dedup/stats/{sample}_edit_distance.tsv'
    params:
        stats = path_to_data + '/{cohort}/results/bowtie2_bam_dedup/stats/{sample}'
    resources: cpus=1, mem_mb=30000, time_min='24:00:00'
    conda: 'conda_env/umi_tools.yml'
    shell:
        "umi_tools dedup -I {input.bam} --output-stats={params.stats} -S {output.bam}"

rule samtools_index_bowtie2:
    input:
        bam = path_to_data + '/{cohort}/results/bowtie2_bam_dedup/{sample}.aligned.sorted.dedup.bam'
    output:
        index = path_to_data + '/{cohort}/results/bowtie2_bam_dedup/{sample}.aligned.sorted.dedup.bam.bai'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        "samtools index {input.bam}"

# ----------------- #
#  QC of BAM files  #
# ----------------- #
def get_paired_or_single_end_unsorted_bam(wildcards):
    paired_data = get_cohort_config(wildcards.cohort)['paired']
    if paired_data:
        return path_to_data + '/' + wildcards.cohort + '/tmp/bwa_mem_extract_barcodes/' + wildcards.sample + '_lib' + wildcards.lib + '.paired.bam'
    else:
        return path_to_data + '/' + wildcards.cohort + '/tmp/bwa_mem/' + wildcards.sample + '_lib' + wildcards.lib + '.single.bam'

def get_all_samples_with_libraries(cohort=None):
    """Returns a list of tuples with cohort name and sample name."""
    samples = get_all_samples(cohort)[['sample_name', 'library_index']]
    return(zip(
        samples.drop_duplicates().sample_name.tolist(),
        samples.drop_duplicates().library_index.tolist()
    ))

def get_aligned_bam(wildcards):
    paired_data = get_cohort_config(wildcards.cohort)['paired']
    umi_data = get_cohort_config(wildcards.cohort)['UMI']
    if paired_data and umi_data:
        return path_to_data + '/' + wildcards.cohort + '/tmp/bwa_mem_umi_tools_extract_pe/' + wildcards.sample + '_lib' + wildcards.lib + '.bam'
    elif paired_data == False and umi_data:
        return path_to_data + '/' + wildcards.cohort + '/tmp/bwa_mem_umi_tools_extract_se/' + wildcards.sample + '_lib' + wildcards.lib + '.bam'
    elif paired_data and umi_data == False:
        return path_to_data + '/' + wildcards.cohort + '/tmp/bwa_mem_pe/' + wildcards.sample + '_lib' + wildcards.lib + '.bam'
    else:
        return path_to_data + '/' + wildcards.cohort + '/tmp/bwa_mem_se/' + wildcards.sample + '_lib' + wildcards.lib + '.bam'
	
# Run FASTQC on final BAM files.
rule fastqc_bam:
    input:
        path_to_data + '/{cohort}/results/bam_dedup/{sample}.aligned.sorted.dedup.bam'
    output:
        html = path_to_data + '/{cohort}/results/qc/bam_dedup_fastqc/{sample}.aligned.sorted.dedup_fastqc.html',
        zipfile = path_to_data + '/{cohort}/results/qc/bam_dedup_fastqc/{sample}.aligned.sorted.dedup_fastqc.zip'
    resources: cpus=1, mem_mb=16000, time_min='24:00:00'
    params:
        outdir = lambda wildcards, output: '/'.join(output.html.split('/')[0:-1])
    conda: 'conda_env/fastqc.yml'
    shell:
        'fastqc --outdir {params.outdir} {input}'

# Run Flagstat to get basic stats on final BAM files.
rule bam_flagstat:
    input:
        aligned = lambda wildcards : get_merged_bam(wildcards),
        deduped = path_to_data + '/{cohort}/results/bam_dedup/{sample}.aligned.sorted.dedup.bam'
    output:
        aligned = path_to_data + '/{cohort}/results/qc/flagstats_bam/{sample}.aligned.sorted.bam.flagstat.tsv',
        deduped = path_to_data + '/{cohort}/results/qc/flagstats_bam/{sample}.aligned.sorted.dedup.bam.flagstat.tsv'
    params:
        qc_path = path_to_data + '/{cohort}/results/qc/flagstats_bam/'
    resources: cpus=1, mem_mb=16000, time_min='1:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        "samtools flagstat {input.aligned} -O tsv > {output.aligned} && samtools flagstat {input.deduped} -O tsv > {output.deduped}"

# Run Flagstat on the unfiltered BAM files
rule unfiltered_bam_flagstat:
    input:
        lambda wildcards : get_aligned_bam(wildcards)
    output:
        path_to_data + '/{cohort}/results/qc/flagstats_bam/{sample}_lib{lib}.bam.flagstat.tsv'
    resources: cpus=1, mem_mb=16000, time_min='1:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        "samtools flagstat {input} -O tsv > {output}"

# Unified QC rule which runs all of the above QCs
rule consolidate_cohort_bam_qc_info:
    input:
        fastqc = lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/qc/bam_dedup_fastqc/{sample}.aligned.sorted.dedup_fastqc.html',
            sample = get_all_samples_list(wildcards.cohort)
        ),
        aligned = lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/qc/flagstats_bam/{sample}.aligned.sorted.bam.flagstat.tsv',
            sample = get_all_samples_list(wildcards.cohort)
        ),
        deduped = lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/qc/flagstats_bam/{sample}.aligned.sorted.dedup.bam.flagstat.tsv',
            sample = get_all_samples_list(wildcards.cohort)
        ),
        lib_stat_files = lambda wildcards: expand(
            [path_to_data + '/' + wildcards.cohort + '/results/qc/flagstats_bam/{sample}_lib{lib}.bam.flagstat.tsv'.format(
                    sample = v[0],
                    lib = v[1]
                ) for v in get_all_samples_with_libraries(wildcards.cohort)
            ]
        )
    output:
        path_to_data + '/{cohort}/results/qc/flagstats_bam/consolidated/flag_stats.tsv'
    params:
        in_path = path_to_data + '/{cohort}/results/qc/flagstats_bam/',
        out_path = path_to_data + '/{cohort}/results/qc/flagstats_bam/consolidated/'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/get_read_count_info.R -i {params.in_path} -o {params.out_path}'

#### Independent for testing
# Run FASTQC on final BAM files.
rule fastqc_bam_markdup:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.dedup.bam',
    output:
        html = path_to_data + '/{cohort}/results/qc/bam_markdup_fastqc/{sample}.aligned.sorted.markdup_fastqc.html',
        zipfile = path_to_data + '/{cohort}/results/qc/bam_markdup_fastqc/{sample}.aligned.sorted.markdup_fastqc.zip'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    params:
        outdir = lambda wildcards, output: '/'.join(output.html.split('/')[0:-1])
    conda: 'conda_env/fastqc.yml'
    shell:
        'fastqc --outdir {params.outdir} {input}'

rule fastqc_bowtie2_bam:
    input:
        path_to_data + '/{cohort}/results/bowtie2_bam_dedup/{sample}.aligned.sorted.dedup.bam'
    output:
        html = path_to_data + '/{cohort}/results/qc/bowtie2_bam_fastqc/{sample}.aligned.sorted.markdup_fastqc.html',
        zipfile = path_to_data + '/{cohort}/results/qc/bowtie2_bam_fastqc/{sample}.aligned.sorted.markdup_fastqc.zip'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    params:
        outdir = lambda wildcards, output: '/'.join(output.html.split('/')[0:-1])
    conda: 'conda_env/fastqc.yml'
    shell:
        'fastqc --outdir {params.outdir} {input}'

# Run Flagstat to get basic stats on final BAM files.
rule bam_flagstat_markdup:
    input:
        aligned = lambda wildcards: get_paired_or_single_end_sorted_merged_bam(wildcards),
        deduped = path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.dedup.bam'
    output:
        aligned = path_to_data + '/{cohort}/results/qc/flagstats_markdup/{sample}.aligned.sorted.bam.flagstat.tsv',
        deduped = path_to_data + '/{cohort}/results/qc/flagstats_markdup/{sample}.aligned.sorted.dedup.bam.flagstat.tsv'
    params:
        qc_path = path_to_data + '/{cohort}/results/qc/flagstats_markdup/'
    resources: cpus=1, mem_mb=8000, time_min='1:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        "samtools flagstat {input.aligned} -O tsv > {output.aligned} && samtools flagstat {input.deduped} -O tsv > {output.deduped}"


rule bowtie2_bam_flagstat:
    input:
        aligned = path_to_data + '/{cohort}/tmp/bowtie2_merge_bam/{sample}.paired.aligned.sorted.bam',
        deduped = path_to_data + '/{cohort}/results/bowtie2_bam_dedup/{sample}.aligned.sorted.dedup.bam'
    output:
        aligned = path_to_data + '/{cohort}/results/qc/bowtie2_flagstats/{sample}.aligned.sorted.bam.flagstat.tsv',
        deduped = path_to_data + '/{cohort}/results/qc/bowtie2_flagstats/{sample}.aligned.sorted.dedup.bam.flagstat.tsv'
    params:
        qc_path = path_to_data + '/{cohort}/results/qc/bowtie2_flagstats/'
    resources: cpus=1, mem_mb=8000, time_min='1:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        "samtools flagstat {input.aligned} -O tsv > {output.aligned} && samtools flagstat {input.deduped} -O tsv > {output.deduped}"

# Run Flagstat on the unfiltered BAM files
rule unfiltered_bam_flagstat_extract_barcodes:
    input:
        lambda wildcards : get_paired_or_single_end_unsorted_bam(wildcards)
    output:
        path_to_data + '/{cohort}/results/qc/flagstats_markdup/{sample}_lib{lib}.bam.flagstat.tsv'
    resources: cpus=1, mem_mb=8000, time_min='1:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        "samtools flagstat {input} -O tsv > {output}"

rule unfiltered_bowtie2_bam_flagstat:
    input:
        path_to_data + '/{cohort}/tmp/bowtie2/{sample}_lib{lib}.paired.bam'
    output:
        path_to_data + '/{cohort}/results/qc/bowtie2_flagstats/{sample}_lib{lib}.bam.flagstat.tsv'
    resources: cpus=1, mem_mb=8000, time_min='1:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        "samtools flagstat {input} -O tsv > {output}"

# Unified QC rule which runs all of the above QCs
rule consolidate_cohort_bam_qc_info_extract_barcodes:
    input:
        fastqc = lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/qc/bam_markdup_fastqc/{sample}.aligned.sorted.markdup_fastqc.html',
            sample = get_all_samples_list(wildcards.cohort)
        ),
        aligned = lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/qc/flagstats_markdup/{sample}.aligned.sorted.bam.flagstat.tsv',
            sample = get_all_samples_list(wildcards.cohort)
        ),
        deduped = lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/qc/flagstats_markdup/{sample}.aligned.sorted.dedup.bam.flagstat.tsv',
            sample = get_all_samples_list(wildcards.cohort)
        ),
        lib_stat_files = lambda wildcards: expand(
            [path_to_data + '/' + wildcards.cohort + '/results/qc/flagstats_markdup/{sample}_lib{lib}.bam.flagstat.tsv'.format(
                    sample = v[0],
                    lib = v[1]
                ) for v in get_all_samples_with_libraries(wildcards.cohort)
            ]
        )
    output:
        path_to_data + '/{cohort}/results/qc/flagstats_markdup/consolidated/flag_stats.tsv'
    params:
        in_path = path_to_data + '/{cohort}/results/qc/flagstats_markdup/',
        out_path = path_to_data + '/{cohort}/results/qc/flagstats_markdup/consolidated/'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/get_read_count_info.R -i {params.in_path} -o {params.out_path}'

rule consolidate_cohort_bowtie2_bam_qc_info:
    input:
        fastqc = lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/qc/bowtie2_bam_fastqc/{sample}.aligned.sorted.markdup_fastqc.html',
            sample = get_all_samples_list(wildcards.cohort)
        ),
        aligned = lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/qc/bowtie2_flagstats/{sample}.aligned.sorted.bam.flagstat.tsv',
            sample = get_all_samples_list(wildcards.cohort)
        ),
        deduped = lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/qc/bowtie2_flagstats/{sample}.aligned.sorted.dedup.bam.flagstat.tsv',
            sample = get_all_samples_list(wildcards.cohort)
        ),
        lib_stat_files = lambda wildcards: expand(
            [path_to_data + '/' + wildcards.cohort + '/results/qc/bowtie2_flagstats/{sample}_lib{lib}.bam.flagstat.tsv'.format(
                    sample = v[0],
                    lib = v[1]
                ) for v in get_all_samples_with_libraries(wildcards.cohort)
            ]
        )
    output:
        path_to_data + '/{cohort}/results/qc/bowtie2_flagstats/consolidated/flag_stats.tsv'
    params:
        in_path = path_to_data + '/{cohort}/results/qc/bowtie2_flagstats/',
        out_path = path_to_data + '/{cohort}/results/qc/bowtie2_flagstats/consolidated/'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/get_read_count_info.R -i {params.in_path} -o {params.out_path}'

# ------------------------------------------------ #
#  MEDIPS Method of methylation count calculation  #
# ------------------------------------------------ #

rule run_MEDIPS:
    input:
        path_to_data + '/{cohort}/results/bam_dedup/{sample}.aligned.sorted.dedup.bam'
    output:
        medips_count = path_to_data + '/{cohort}/results/MEDIPS/{sample}.medips_output.tsv',
        medips_qc = path_to_data + '/{cohort}/results/MEDIPS/{sample}.QCStats_matrix.tsv'
    params:
        paired_val = lambda wildcards: get_cohort_config(wildcards.cohort)['paired'],
        extnd = lambda wildcards: get_cohort_config(wildcards.cohort)['extend_dist'],
        shift = 0,
        uniq = 0,
        binwidth = lambda wildcards: get_cohort_config(wildcards.cohort)['bin_width'],
        sample_name = '{sample}'
    resources: partition='himem', cpus=1, mem_mb=61440, time_min='5-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/run_MEDIPS.R -b {input} -n {params.sample_name} -o {output.medips_count} -q {output.medips_qc} -p {params.paired_val} -e {params.extnd} -s {params.shift} -u {params.uniq} -w {params.binwidth}'

rule consolidate_cohort_MEDIPS:
    input:
        lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/MEDIPS/{sample}.medips_output.tsv',
            sample = get_all_samples_list(wildcards.cohort)
        )
    output:
        path_to_data + '/{cohort}/results/MEDIPS/consolidated/MEDIPS_Counts.tsv',
        path_to_data + '/{cohort}/results/MEDIPS/consolidated/MEDIPS_CPM.tsv',
        path_to_data + '/{cohort}/results/MEDIPS/consolidated/MEDIPS_Counts.rds',
        path_to_data + '/{cohort}/results/MEDIPS/consolidated/MEDIPS_CPM.rds'
    params:
        data = 'MEDIPS',
        in_path = path_to_data + '/{cohort}/results/MEDIPS/',
        out_path = path_to_data + '/{cohort}/results/MEDIPS/consolidated/'
    resources: cpus=1, mem_mb=30000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/consolidate_cohort_samples.R -i {params.in_path} -o {params.out_path} -d {params.data}'

rule consolidate_cohort_MEDIPS_QC:
    input:
        lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/MEDIPS/{sample}.QCStats_matrix.tsv',
            sample = get_all_samples_list(wildcards.cohort)
        )
    output:
        path_to_data + '/{cohort}/results/MEDIPS/consolidated/MEDIPS_QCStats.tsv'
    params:
        in_path = path_to_data + '/{cohort}/results/MEDIPS/',
        out_path = path_to_data + '/{cohort}/results/MEDIPS/consolidated/'
    resources: cpus=1, mem_mb=16000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/consolidate_cohort_MEDIPS_QC.R -i {params.in_path} -o {params.out_path}'


### Independent for testing
# Extract Barcode markdup
rule run_MEDIPS_markdup:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.dedup.bam'
    output:
        medips_count = path_to_data + '/{cohort}/results/MEDIPS_markdup/{sample}.medips_output.tsv',
        medips_qc = path_to_data + '/{cohort}/results/MEDIPS_markdup/{sample}.QCStats_matrix.tsv'
    params:
        paired_val = lambda wildcards: get_cohort_config(wildcards.cohort)['paired'],
        extnd = lambda wildcards: get_cohort_config(wildcards.cohort)['extend_dist'],
        shift = 0,
        uniq = 0,
        binwidth = lambda wildcards: get_cohort_config(wildcards.cohort)['bin_width'],
        sample_name = '{sample}'
    resources: partition='himem', cpus=1, mem_mb=61440, time_min='5-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/run_MEDIPS.R -b {input} -n {params.sample_name} -o {output.medips_count} -q {output.medips_qc} -p {params.paired_val} -e {params.extnd} -s {params.shift} -u {params.uniq} -w {params.binwidth}'

rule consolidate_cohort_MEDIPS_markdup:
    input:
        lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/MEDIPS_markdup/{sample}.medips_output.tsv',
            sample = get_all_samples_list(wildcards.cohort)
        )
    output:
        path_to_data + '/{cohort}/results/MEDIPS_markdup/consolidated/MEDIPS_Counts.tsv',
        path_to_data + '/{cohort}/results/MEDIPS_markdup/consolidated/MEDIPS_CPM.tsv',
        path_to_data + '/{cohort}/results/MEDIPS_markdup/consolidated/MEDIPS_Counts.rds',
        path_to_data + '/{cohort}/results/MEDIPS_markdup/consolidated/MEDIPS_CPM.rds'
    params:
        data = 'MEDIPS',
        in_path = path_to_data + '/{cohort}/results/MEDIPS_markdup/',
        out_path = path_to_data + '/{cohort}/results/MEDIPS_markdup/consolidated/'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/consolidate_cohort_samples.R -i {params.in_path} -o {params.out_path} -d {params.data}'

rule consolidate_cohort_MEDIPS_QC_markdup:
    input:
        lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/MEDIPS_markdup/{sample}.QCStats_matrix.tsv',
            sample = get_all_samples_list(wildcards.cohort)
        )
    output:
        path_to_data + '/{cohort}/results/MEDIPS_markdup/consolidated/MEDIPS_QCStats.tsv'
    params:
        in_path = path_to_data + '/{cohort}/results/MEDIPS_markdup/',
        out_path = path_to_data + '/{cohort}/results/MEDIPS_markdup/consolidated/'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/consolidate_cohort_MEDIPS_QC.R -i {params.in_path} -o {params.out_path}'

# UMI_tools dedup
rule run_MEDIPS_umi_tools_dedup:
    input:
        path_to_data + '/{cohort}/results/bam_umi_tools_dedup/{sample}.aligned.sorted.dedup.bam'
    output:
        medips_count = path_to_data + '/{cohort}/results/MEDIPS_umi_tools_dedup/{sample}.medips_output.tsv',
        medips_qc = path_to_data + '/{cohort}/results/MEDIPS_umi_tools_dedup/{sample}.QCStats_matrix.tsv'
    params:
        paired_val = lambda wildcards: get_cohort_config(wildcards.cohort)['paired'],
        extnd = lambda wildcards: get_cohort_config(wildcards.cohort)['extend_dist'],
        shift = 0,
        uniq = 0,
        binwidth = lambda wildcards: get_cohort_config(wildcards.cohort)['bin_width'],
        sample_name = '{sample}'
    resources: partition='himem', cpus=1, mem_mb=61440, time_min='5-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/run_MEDIPS.R -b {input} -n {params.sample_name} -o {output.medips_count} -q {output.medips_qc} -p {params.paired_val} -e {params.extnd} -s {params.shift} -u {params.uniq} -w {params.binwidth}'

rule consolidate_cohort_MEDIPS_umi_tools_dedup:
    input:
        lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/MEDIPS_umi_tools_dedup/{sample}.medips_output.tsv',
            sample = get_all_samples_list(wildcards.cohort)
        )
    output:
        path_to_data + '/{cohort}/results/MEDIPS_umi_tools_dedup/consolidated/MEDIPS_Counts.tsv',
        path_to_data + '/{cohort}/results/MEDIPS_umi_tools_dedup/consolidated/MEDIPS_CPM.tsv',
        path_to_data + '/{cohort}/results/MEDIPS_umi_tools_dedup/consolidated/MEDIPS_Counts.rds',
        path_to_data + '/{cohort}/results/MEDIPS_umi_tools_dedup/consolidated/MEDIPS_CPM.rds'
    params:
        data = 'MEDIPS',
        in_path = path_to_data + '/{cohort}/results/MEDIPS_umi_tools_dedup/',
        out_path = path_to_data + '/{cohort}/results/MEDIPS_umi_tools_dedup/consolidated/'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/consolidate_cohort_samples.R -i {params.in_path} -o {params.out_path} -d {params.data}'

rule consolidate_cohort_MEDIPS_QC_umi_tools_dedup:
    input:
        lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/MEDIPS_umi_tools_dedup/{sample}.QCStats_matrix.tsv',
            sample = get_all_samples_list(wildcards.cohort)
        )
    output:
        path_to_data + '/{cohort}/results/MEDIPS_umi_tools_dedup/consolidated/MEDIPS_QCStats.tsv'
    params:
        in_path = path_to_data + '/{cohort}/results/MEDIPS_umi_tools_dedup/',
        out_path = path_to_data + '/{cohort}/results/MEDIPS_umi_tools_dedup/consolidated/'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/consolidate_cohort_MEDIPS_QC.R -i {params.in_path} -o {params.out_path}'

rule run_MEDIPS_Nick:
    input:
        path_to_data + '/{cohort}/results/bowtie2_bam_dedup/{sample}.aligned.sorted.markdup.bam'
    output:
        medips_count = path_to_data + '/{cohort}/results/MEDIPS_Nick/{sample}.medips_output.tsv',
        medips_qc = path_to_data + '/{cohort}/results/MEDIPS_Nick/{sample}.QCStats_matrix.tsv'
    params:
        paired_val = lambda wildcards: get_cohort_config(wildcards.cohort)['paired'],
        extnd = 300,
        shift = 0,
        uniq = 0,
        binwidth = 300,
        sample_name = '{sample}'
    resources: partition='himem', cpus=1, mem_mb=61440, time_min='5-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/run_MEDIPS.R -b {input} -n {params.sample_name} -o {output.medips_count} -q {output.medips_qc} -p {params.paired_val} -e {params.extnd} -s {params.shift} -u {params.uniq} -w {params.binwidth}'

rule consolidate_cohort_MEDIPS_Nick:
    input:
        lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/MEDIPS_Nick/{sample}.medips_output.tsv',
            sample = get_all_samples_list(wildcards.cohort)
        )
    output:
        path_to_data + '/{cohort}/results/MEDIPS_Nick/consolidated/MEDIPS_Counts.tsv',
        path_to_data + '/{cohort}/results/MEDIPS_Nick/consolidated/MEDIPS_CPM.tsv',
        path_to_data + '/{cohort}/results/MEDIPS_Nick/consolidated/MEDIPS_Counts.rds',
        path_to_data + '/{cohort}/results/MEDIPS_Nick/consolidated/MEDIPS_CPM.rds'
    params:
        data = 'MEDIPS',
        in_path = path_to_data + '/{cohort}/results/MEDIPS_Nick/',
        out_path = path_to_data + '/{cohort}/results/MEDIPS_Nick/consolidated/'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/consolidate_cohort_samples.R -i {params.in_path} -o {params.out_path} -d {params.data}'

rule consolidate_cohort_MEDIPS_QC_Nick:
    input:
        lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/MEDIPS_Nick/{sample}.QCStats_matrix.tsv',
            sample = get_all_samples_list(wildcards.cohort)
        )
    output:
        path_to_data + '/{cohort}/results/MEDIPS_Nick/consolidated/MEDIPS_QCStats.tsv'
    params:
        in_path = path_to_data + '/{cohort}/results/MEDIPS_Nick/',
        out_path = path_to_data + '/{cohort}/results/MEDIPS_Nick/consolidated/'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/consolidate_cohort_MEDIPS_QC.R -i {params.in_path} -o {params.out_path}'

# --------------------------------------------- #
#  MeDEStrand Method of methylation correction  #
# --------------------------------------------- #

rule run_MeDEStrand:
    input:
        path_to_data + '/{cohort}/results/bam_dedup/{sample}.aligned.sorted.dedup.bam'
    output:
        path_to_data + '/{cohort}/results/MeDEStrand/{sample}.medestrand_output.tsv'
    params:
        medestrand = config['paths']['dependencies']['medestrand_path'],
        paired_val = lambda wildcards: get_cohort_config(wildcards.cohort)['paired'],
        extnd = lambda wildcards: get_cohort_config(wildcards.cohort)['extend_dist'],
        shift = 0,
        uniq = 0,
        binwidth = lambda wildcards: get_cohort_config(wildcards.cohort)['bin_width'],
        sample_name = '{sample}'
    resources: partition='himem', cpus=1, mem_mb=61440, time_min='5-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/run_medestrand.R -b {input} -n {params.sample_name} -o {output} -p {params.paired_val} -m {params.medestrand} -e {params.extnd} -s {params.shift} -u {params.uniq} -w {params.binwidth}'
    
rule consolidate_cohort_MeDEStrand:
    input:
        lambda wildcards: expand(
            path_to_data + '/' + wildcards.cohort + '/results/MeDEStrand/{sample}.medestrand_output.tsv',
            sample = get_all_samples_list(wildcards.cohort)
        )
    output:
        path_to_data + '/{cohort}/results/MeDEStrand/consolidated/MeDEStrand_AbsMethyl.tsv',
        path_to_data + '/{cohort}/results/MeDEStrand/consolidated/MeDEStrand_AbsMethyl.rds'
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
        path_to_data + '/samples/{sample}/merged/bwa_mem/aligned.sorted.dedup.bam'
    output:
        qsea_count = path_to_data + '/samples/{sample}/merged/QSEA/by_chromosome/qsea_count_{chrom}_output.tsv',
        qsea_beta = path_to_data + '/samples/{sample}/merged/QSEA/by_chromosome/qsea_beta_{chrom}_output.tsv',
        qsea_qc = path_to_data + '/samples/{sample}/merged/QSEA/by_chromosome/QCStats_{chrom}_matrix.tsv'
    params:
        out = path_to_data + '/samples/{sample}/merged/QSEA/by_chromosome/'
    resources: partition='himem', cpus=4, mem_mb=61440, time_min='3-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/run_QSEA.R -s {wildcards.sample} -c {wildcards.chrom} -b {input} -o {params.out} --count {output.qsea_count} --beta {output.qsea_beta} --qc {output.qsea_qc}'

rule merge_QSEA_count:
    input:
        lambda wildcards: [path_to_data + '/samples/{{sample}}/merged/QSEA/by_chromosome/qsea_count_{chrom}_output.tsv'.format(chrom=a) for a in chromosomes[wildcards.cohort]]
    output:
        path_to_data + '/samples/{sample}/merged/QSEA/qsea_count_output.tsv'
    resources: partition='himem', cpus=1, mem_mb=8000, time_min='24:00:00'
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
    resources: partition='himem', cpus=1, mem_mb=8000, time_min='24:00:00'
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
    resources: partition='himem', cpus=1, mem_mb=8000, time_min='24:00:00'
    run:
        for i, input_file in enumerate(input):
            input_data = pd.read_csv(input_file, delimiter='\t', comment='#')
            if i == 0:
                input_data.to_csv(output[0], header=True, sep='\t', index=False)
            else:
                input_data.to_csv(output[0], header=False, sep='\t', index=False, mode='a')

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
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.dedup.bam',
    output:
        binstat = temp(path_to_data + '/{cohort}/tmp/bam_bin_stats/bin_stats_{sample}_{species}_{chrom}.tsv'),
        filtered = path_to_data + '/{cohort}/results/bam_bin_stats/removed_bins_{sample}_{species}_{chrom}.tsv',
    params:
        bsgenome = lambda wildcards: get_cohort_config(wildcards.cohort)['bsgenome'][wildcards.species],
        bsgenome_chr = lambda wildcards: get_bsgenome_chrom(wildcards.species, wildcards.chrom)
    resources: partition='himem', cpus=1, mem_mb=61440, time_min='24:00:00'
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
    resources: partition='himem', cpus=1, mem_mb=8000, time_min='24:00:00'
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
    resources: partition='himem', cpus=1, time_min='5-00:00:00', mem_mb=lambda wildcards, attempt: 61440 if attempt == 1 else 32000
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/cfmedip_nbglm.R -i {input} -o {output.fit} --modelout {output.model}'

# Below is the full bayesian approach for fitting, which remains under development
rule fit_bin_stats:
    input:
        path_to_data + '/{cohort}/results/merge_bin_stats/bin_stats_{sample}.feather'
    output:
        path_to_data + '/{cohort}/results/fit_bin_stats/{sample}_fit_{method}.tsv'
    resources: partition='himem', cpus=1, mem_mb=61440, time_min='5-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/fit_cpg_bins.R -i {input} -o {output} --method {wildcards.method}'

# ------------------------------------------------ #
#        Delfi Fragment Profile Calculation  	   #
# ------------------------------------------------ #

rule run_delfi_fragment_profile:
    input:
        path_to_data + '/{cohort}/results/bam_dedup/{sample}.aligned.sorted.dedup.bam'
    output:
        path_to_data + '/{cohort}/results/delfi/{sample}_100kb_fragment_profile.rds'
    params:
        paired_val = lambda wildcards: get_cohort_config(wildcards.cohort)['paired'],
        out_path = path_to_data + '/{cohort}/results/delfi',
        data_path = lambda wildcards: get_cohort_config(wildcards.cohort)['window_assets_dir'],
        sample_name = '{sample}'
    resources: partition='all', cpus=1, mem_mb=30000, time_min='12:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/delfi_fragment_profile.R -b {input} -s {params.sample_name} -o {params.out_path} -d $dataDir -p {params.paired_val}'