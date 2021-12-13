import yaml
import pandas as pd
import gzip

configfile: "config.yml"

for cohort in config['data']['cohorts']:
    if config['data']['cohorts'][cohort]['active']:
        break

path_to_data = config['data']['base_path'] + '/' + config['data']['cohorts'][cohort]['cohort_name']
b_pattern = config['data']['cohorts'][cohort]['bpattern']
b_list = config['data']['cohorts'][cohort]['barcodes']

def get_cohort_data(cohort_name):
    samplesheet = pd.read_csv(config['data']['cohorts'][cohort_name]['samplesheet'], comment='#').drop_duplicates()
    samplesheet = samplesheet.sort_values(by=['sample_name', 'read_in_pair'])
    return samplesheet

def get_fastq_path(sample, library, read_in_pair=1):
    library = int(library)
    all_samples = get_all_samples()
    sample_line = all_samples[
        (all_samples.sample_name == sample) &
        (all_samples.library_index == library) &
        (all_samples.read_in_pair == read_in_pair)
    ]
    return sample_line.path.to_list()[0]

def get_all_samples():
    all_samples = pd.concat([
        get_cohort_data(cohort_name)
        for cohort_name
        in config['data']['cohorts']
        if config['data']['cohorts'][cohort_name]['active']
        ])
    return all_samples

def get_all_samples_list():
    return get_all_samples().sample_name.unique().tolist()

def clean(command):
    command = command.replace('\n', ' ').replace('\t', ' ')
    while '  ' in command:
        command = command.replace('  ', ' ')
    return command.strip()

rule all:
    input:
        #flag stats of the bam files
        expand(
            path_to_data + '/samples/{sample}/merged/bwa_mem/all_flagstats.txt',
            sample = get_all_samples_list()
        ),
        path_to_data + '/flagstats/flag_stats.tsv',
        #run MeD-ReMix
        expand(
            path_to_data + '/samples/{sample}/merged/bin_stats/bin_stats_fit_nbglm.tsv',
            sample = get_all_samples_list()
        ),
        #run QSEA
        #expand(
        #    path_to_data + '/samples/{sample}/merged/QSEA/qsea_count_output.tsv',
        #    sample = get_all_samples_list()
        #),
        #expand(
        #    path_to_data + '/samples/{sample}/merged/QSEA/qsea_beta_output.tsv',
        #    sample = get_all_samples_list()
        #),
        #expand(
        #    path_to_data + '/samples/{sample}/merged/QSEA/QCStats_matrix.tsv',
        #    sample = get_all_samples_list()
        #),
        #run MeDEStrand
        expand(
            path_to_data + '/samples/{sample}/merged/MeDEStrand/medestrand_output.tsv',
            sample = get_all_samples_list()
        ),
        #consolidate
        path_to_data + '/MeDEStrand/MeDEStrand_AbsMethyl.tsv',
        #run MEDPIS
        expand(
            path_to_data + '/samples/{sample}/merged/MEDIPS/medips_output.tsv',
            sample = get_all_samples_list()
        ),
        #consolidate
        path_to_data + '/MEDIPS/MEDIPS_Counts.tsv',
        ##clean intermediary libraries if they are taking up too much space
        #expand(
        #    path_to_data + '/samples/{sample}/libraries/fastq_libs_cleaned.txt',
        #    sample = get_all_samples_list()
        #),
        #expand(
        #    path_to_data + '/samples/{sample}/libraries/bam_libs_cleaned.txt',
        #    sample = get_all_samples_list()
        #),
        #expand(
        #    path_to_data + '/samples/{sample}/merged/bin_stats/chrom_bin_stats_libs_cleaned.txt',
        #    sample = get_all_samples_list()
        #),
        #expand(
        #    path_to_data + '/samples/{sample}/merged/QSEA/chrom_qsea_libs_cleaned.txt',
        #    sample = get_all_samples_list()
        #)


rule set_environment:
    output:
        'environment_is_set'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'touch {output}'

rule gunzip_fastq:
    input:
        lambda wildcards: get_fastq_path(wildcards.sample, int(wildcards.lib), int(wildcards.read))
    output:
        path_to_data + '/samples/{sample}/libraries/lib_{lib}/fastq/R{read}.fastq'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    shell:
        'gunzip -c {input} > {output}'

rule extract_barcodes:
    input:
        R1 =    path_to_data + '/samples/{sample}/libraries/lib_{lib}/fastq/R1.fastq',
        R2 =    path_to_data + '/samples/{sample}/libraries/lib_{lib}/fastq/R2.fastq'
    output:
        R1 =    path_to_data + '/samples/{sample}/libraries/lib_{lib}/fastq/extract_barcode_R1.fastq.gz',
        R2 =    path_to_data + '/samples/{sample}/libraries/lib_{lib}/fastq/extract_barcode_R2.fastq.gz'
    params:
        barcode_out_R1 = path_to_data + '/samples/{sample}/libraries/lib_{lib}/fastq/extract_barcode_R1.fastq',
        barcode_out_R2 = path_to_data + '/samples/{sample}/libraries/lib_{lib}/fastq/extract_barcode_R2.fastq',
        outprefix = lambda wildcards, output: output.R1.split('_barcode_')[0],
        bpattern = b_pattern,
        blist = b_list
    resources: cpus=1, mem_mb=16000, time_min='1-00:00:00'
    run:
        if params.bpattern is None:
            if params.blist is None:
                shell("gzip -c {input.R1} > {output.R1}")
                shell("gzip -c {input.R2} > {output.R2}")
                shell("rm {input.R1} {input.R2}")
            else:
                shell("python src/ConsensusCruncher/ConsensusCruncher/extract_barcodes.py --read1 {input.R1} --read2 {input.R2} --outfile {params.outprefix} --blist {params.blist}")
                shell("gzip {params.barcode_out_R1} {params.barcode_out_R2}")
                shell("rm {input.R1} {input.R2}")
        else:
            shell("python src/ConsensusCruncher/ConsensusCruncher/extract_barcodes.py --read1 {input.R1} --read2 {input.R2} --outfile {params.outprefix} --bpattern {params.bpattern} --blist {params.blist}")
            shell("gzip {params.barcode_out_R1} {params.barcode_out_R2}")
            shell("rm {input.R1} {input.R2}")

def get_read_group_from_fastq(fastq_file, sample_name):
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

rule bwa_mem:
    input:
        path_to_data + '/samples/{sample}/libraries/lib_{lib}/fastq/extract_barcode_R1.fastq.gz',
        path_to_data + '/samples/{sample}/libraries/lib_{lib}/fastq/extract_barcode_R2.fastq.gz'
    output:
        path_to_data + '/samples/{sample}/libraries/lib_{lib}/bwa_mem/aligned.bam'
    resources: cpus=4, mem_mb=16000, time_min='72:00:00'
    params:
        read_group = lambda wildcards, input: get_read_group_from_fastq(
            fastq_file = get_fastq_path(wildcards.sample, wildcards.lib),
            sample_name = wildcards.sample
        ),
        bwa_index = config['paths']['bwa_index']
    shell:
        "bwa mem -M -t4 -R'{params.read_group}' {params.bwa_index} {input} | samtools view -bSo {output}"

rule clean_fastq_libs:
    output:
        path_to_data + '/samples/{sample}/libraries/fastq_libs_cleaned.txt'
    params:
        fastq_paths = lambda wildcards: expand(
                path_to_data + '/samples/' + wildcards.sample + '/libraries/lib_{lib}/fastq/',
                lib=get_libraries_of_sample(wildcards.sample)
        )
    resources: cpus=1, mem_mb=8000, time_min='10:00:00'
    shell:
        "rm -r {params.fastq_paths} && touch {output}"

rule bam_to_sorted_bam:
    input:
        path_to_data + '/samples/{sample}/libraries/lib_{lib}/bwa_mem/aligned.bam'
    output:
        bam=path_to_data + '/samples/{sample}/libraries/lib_{lib}/bwa_mem/aligned.sorted.bam',
        index=path_to_data + '/samples/{sample}/libraries/lib_{lib}/bwa_mem/aligned.sorted.bam.bai'
    resources: cpus=32, mem_mb=30000, time_min='72:00:00'
    shell:
        # Try it without fixmate for replication purposes
        #"samtools view -buS -f 2 -F 4 -@4 {input} | samtools fixmate -m - - | samtools sort -@4 -o {output.bam} && samtools index {output.bam}"
        clean(r'''
        samtools view -buS -f 2 -F 4 -@32 {input} |
        samtools fixmate -m - - |
        samtools sort -@32 -o {output.bam} && samtools index {output.bam}
        ''')

def get_libraries_of_sample(sample):
    filtered_table = get_all_samples()[get_all_samples().sample_name == sample]
    return(list(set(filtered_table.library_index.to_list())))

rule merge_bam:
    input:
        lambda wildcards: expand(
                path_to_data + '/samples/' + wildcards.sample + '/libraries/lib_{lib}/bwa_mem/aligned.sorted.bam',
                lib=get_libraries_of_sample(wildcards.sample)
        )
    output:
        path_to_data + '/samples/{sample}/merged/bwa_mem/aligned.sorted.bam'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    shell:
        'samtools merge {output} {input} && samtools index {output}'

rule clean_bam_libs:		
    output:
        path_to_data + '/samples/{sample}/libraries/bam_libs_cleaned.txt'
    params:
        bam_paths = lambda wildcards: expand(
                path_to_data + '/samples/' + wildcards.sample + '/libraries/lib_{lib}/bwa_mem/',
                lib=get_libraries_of_sample(wildcards.sample)
        )
    resources: cpus=1, mem_mb=8000, time_min='10:00:00'
    shell:
        "rm -r {params.bam_paths} && touch {output}"

rule bam_markdup:
    input:
        path_to_data + '/samples/{sample}/merged/bwa_mem/aligned.sorted.bam'
    output:
        bam=path_to_data + '/samples/{sample}/merged/bwa_mem/aligned.sorted.markdup.bam',
        index=path_to_data + '/samples/{sample}/merged/bwa_mem/aligned.sorted.markdup.bam.bai'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    resources: mem_mb=8000, time_min='72:00:00'
    shell:
        "samtools markdup -r {input} {output.bam} && samtools index {output.bam}"

rule bam_flagstats:
    input:
        all = path_to_data + '/samples/{sample}/merged/bwa_mem/aligned.sorted.bam',
        deduped = path_to_data + '/samples/{sample}/merged/bwa_mem/aligned.sorted.markdup.bam'
    output:
        all=path_to_data + '/samples/{sample}/merged/bwa_mem/all_flagstats.txt',
        deduped=path_to_data + '/samples/{sample}/merged/bwa_mem/deduped_flagstats.txt'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    shell:
        "samtools flagstat {input.all} > {output.all} && samtools flagstat {input.deduped} > {output.deduped}"

def get_bsgenome_chrom(species, chrom):
    chrom_map = {'1': 'Chr1', '3': 'Chr3'}
    if species == 'human':
        return chrom
    elif species == 'arabidopsis':
        return chrom_map[chrom]

rule bam_bin_stats:
    input:
        path_to_data + '/samples/{sample}/merged/bwa_mem/aligned.sorted.markdup.bam'
    output:
        binstat=path_to_data + '/samples/{sample}/merged/bin_stats/by_chromosome/bin_stats_{species}_{chrom}.tsv',
        filtered=path_to_data + '/samples/{sample}/merged/bin_stats/filtered_out/bin_stats_{species}_{chrom}.tsv'
    params:
        bsgenome = lambda wildcards:config['paths']['bsgenome'][wildcards.species],
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

MAJOR_HUMAN_CHROMOSOMES = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
all_chromosome_tuples = [('human', c) for c in MAJOR_HUMAN_CHROMOSOMES] + [('arabidopsis', '1'), ('arabidopsis', '3')]

rule merge_bin_stats:
    input:
        [path_to_data + '/samples/{{sample}}/merged/bin_stats/by_chromosome/bin_stats_{species}_{chrom}.tsv'.format(species=a[0], chrom=a[1]) for a in all_chromosome_tuples]
    output:
        path_to_data + '/samples/{sample}/merged/bin_stats/bin_stats.tsv'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    run:
        for i, input_file in enumerate(input):
            input_data = pd.read_csv(input_file, delimiter='\t', comment='#')
            if i == 0:
                input_data.to_csv(output[0], header=True, sep='\t', index=False)
            else:
                input_data.to_csv(output[0], header=False, sep='\t', index=False, mode='a')

rule clean_chrom_bin_stats_libs:
    output:
        path_to_data + '/samples/{sample}/merged/bin_stats/chrom_bin_stats_libs_cleaned.txt'
    params:
        binstat = path_to_data + '/samples/{sample}/merged/bin_stats/by_chromosome',
        filtered = path_to_data + '/samples/{sample}/merged/bin_stats/filtered_out'
    resources: cpus=1, mem_mb=8000, time_min='10:00:00'
    shell:
        "rm -r {params.binstat} {params.filtered} && touch {output}"

# This is the currently used negative binomial GLM approach to fitting
rule cfmedip_nbglm:
    input:
        path_to_data + '/samples/{sample}/merged/bin_stats/bin_stats.tsv'
    output:
        fit=path_to_data + '/samples/{sample}/merged/bin_stats/bin_stats_fit_nbglm.tsv',
        model=path_to_data + '/samples/{sample}/merged/bin_stats/bin_stats_model_nbglm.Rds'
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/cfmedip_nbglm.R -i {input} -o {output.fit} --modelout {output.model}'

rule run_QSEA:
    input:
        path_to_data + '/samples/{sample}/merged/bwa_mem/aligned.sorted.markdup.bam'
    output:
        qsea_count=path_to_data + '/samples/{sample}/merged/QSEA/by_chromosome/qsea_count_{chrom}_output.tsv',
        qsea_beta=path_to_data + '/samples/{sample}/merged/QSEA/by_chromosome/qsea_beta_{chrom}_output.tsv',
        qsea_qc=path_to_data + '/samples/{sample}/merged/QSEA/by_chromosome/QCStats_{chrom}_matrix.tsv'
    params:
        out = path_to_data + '/samples/{sample}/merged/QSEA/by_chromosome/'
    resources: cpus=4, mem_mb=30000, time_min='3-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/run_QSEA.R -s {wildcards.sample} -c {wildcards.chrom} -b {input} -o {params.out} --count {output.qsea_count} --beta {output.qsea_beta} --qc {output.qsea_qc}'

rule merge_QSEA_count:
    input:
        [path_to_data + '/samples/{{sample}}/merged/QSEA/by_chromosome/qsea_count_{chrom}_output.tsv'.format(chrom=a) for a in MAJOR_HUMAN_CHROMOSOMES]
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
        [path_to_data + '/samples/{{sample}}/merged/QSEA/by_chromosome/qsea_beta_{chrom}_output.tsv'.format(chrom=a) for a in MAJOR_HUMAN_CHROMOSOMES]
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
        [path_to_data + '/samples/{{sample}}/merged/QSEA/by_chromosome/QCStats_{chrom}_matrix.tsv'.format(chrom=a) for a in MAJOR_HUMAN_CHROMOSOMES]
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

rule clean_chrom_qsea_libs:
    output:
        path_to_data + '/samples/{sample}/merged/QSEA/chrom_qsea_libs_cleaned.txt'
    params:
        chrom_paths = path_to_data + '/samples/{sample}/merged/QSEA/by_chromosome/'
    resources: cpus=1, mem_mb=8000, time_min='10:00:00'
    shell:
        "rm -r {params.chrom_paths} && touch {output}"

rule run_MEDIPS:
    input:
        path_to_data + '/samples/{sample}/merged/bwa_mem/aligned.sorted.markdup.bam'
    output:
        medips_count=path_to_data + '/samples/{sample}/merged/MEDIPS/medips_output.tsv',
        medips_qc=path_to_data + '/samples/{sample}/merged/MEDIPS/QCStats_matrix.tsv'
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/run_MEDIPS.R -b {input} -o {output.medips_count} -q {output.medips_qc}'

rule run_medestrand:
    input:
        path_to_data + '/samples/{sample}/merged/bwa_mem/aligned.sorted.markdup.bam'
    output:
        path_to_data + '/samples/{sample}/merged/MeDEStrand/medestrand_output.tsv'
    params:
        medestrand = config['paths']['dependencies']['medestrand_path']
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/run_medestrand.R -b {input} -o {output} -m {params.medestrand}'

rule consolidate_cohort_MEDIPS:
    input:
        expand(
            path_to_data + '/samples/{sample}/merged/MEDIPS/medips_output.tsv',
            sample = get_all_samples_list()
        )
    output:
        path_to_data + '/MEDIPS/MEDIPS_Counts.tsv',
        path_to_data + '/MEDIPS/MEDIPS_CPM.tsv'
    params:
        data = 'MEDIPS',
        out_path = path_to_data + '/MEDIPS/',
        cohort_path = path_to_data
    resources: cpus=1, mem_mb=30000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/consolidate_cohort_samples.R -i {params.cohort_path} -o {params.out_path} -d {params.data}'

rule consolidate_cohort_MeDEStrand:
    input:
        expand(
            path_to_data + '/samples/{sample}/merged/MeDEStrand/medestrand_output.tsv',
            sample = get_all_samples_list()
        )
    output:
        path_to_data + '/MeDEStrand/MeDEStrand_AbsMethyl.tsv'
    params:
        data = 'MeDEStrand',
        out_path = path_to_data + '/MeDEStrand/',
        cohort_path = path_to_data
    resources: cpus=1, mem_mb=30000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/consolidate_cohort_samples.R -i {params.cohort_path} -o {params.out_path} -d {params.data}'

rule consolidate_cohort_duplication_rate:
    input:
        all = expand(
            path_to_data + '/samples/{sample}/merged/bwa_mem/all_flagstats.txt',
            sample = get_all_samples_list()
        ),
        deduped = expand(
            path_to_data + '/samples/{sample}/merged/bwa_mem/deduped_flagstats.txt',
            sample = get_all_samples_list()
        )
    output:
        path_to_data + '/flagstats/flag_stats.tsv'
    params:
        out_path = path_to_data + '/flagstats/',
        cohort_path = path_to_data
    resources: cpus=1, mem_mb=30000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/get_duplication_rate.R -i {params.cohort_path} -o {params.out_path}'

# Below is the full bayesian approach for fitting, which remains under development
rule fit_bin_stats:
    input:
        path_to_data + '/samples/{sample}/merged/bin_stats/bin_stats.tsv'
    output:
        path_to_data + '/samples/{sample}/merged/bin_stats/bin_stats_fit_{method}.Rds'
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/fit_cpg_bins.R -i {input} -o {output} --method {wildcards.method}'
