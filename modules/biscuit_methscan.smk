# rule all:
#     input:
#         #join(config['outdir'], f'biscuit/p1_A1/p1_A1_HCG.bed') # test pipeline for one sample
#         #[join(config['outidr'], f'biscuit/{sample}/{sample}_HCG.bed') for sn in metadata.index]
#         join(config['outdir'], 'methscan/VMR_matrices'),
#         multiext(join(config['outdir'], f'qc_plots/{config['dataset']}'), '_CpG_count_vs_meth_frac.png')


# htslib creates a FASTA index beside its input; stage it inside the results tree.
rule prepare_biscuit_reference:
    input:
        ancient(config['reference']['genome'])
    output:
        genome = join(config['outdir'], 'reference/genome.fa'),
        fai = join(config['outdir'], 'reference/genome.fa.fai')
    conda: '../envs/biscuit.yaml'
    shell:
        '''
        ln -sf {input:q} {output.genome:q}
        samtools faidx {output.genome:q}
        '''


if not config['reference'].get('biscuit_index'):
    rule prepare_biscuit_index:
        input:
            ancient(rules.prepare_biscuit_reference.output.genome)
        output:
            directory(join(config['outdir'], 'reference/biscuit_index'))
        params:
            prefix = config['biscuit_index_prefix'],
            alg = config['biscuit_index_alg']
        log: join(config['outdir'], 'reference/biscuit_index.log')
        conda: '../envs/biscuit.yaml'
        shell:
            '''
            mkdir -p {output:q}
            biscuit index {input:q} -p {params.prefix:q} {params.alg} > {log:q} 2>&1
            '''


rule trim_adaptors:
    input:
        read1 = ancient(lambda wildcards: fqid_to_reads[wildcards.fqid]['read1']),
        read2 = ancient(lambda wildcards: fqid_to_reads[wildcards.fqid]['read2'])
    output:
        read1 = temp(join(config['outdir'], 'biscuit/{sample}/{fqid}_R1_val_1.fq.gz')),
        read2 = temp(join(config['outdir'], 'biscuit/{sample}/{fqid}_R2_val_2.fq.gz')),
        report1 = join(config['outdir'], 'biscuit/{sample}/{fqid}_R1.fastq.gz_trimming_report.txt'),
        report2 = join(config['outdir'], 'biscuit/{sample}/{fqid}_R2.fastq.gz_trimming_report.txt')
    params:
        outdir = join(config['outdir'], 'biscuit/{sample}')
    threads: 2
    resources:
        mem_mb=4000,
        walltime=120
    log: join(config['outdir'], 'biscuit/{sample}/{fqid}.trim_galore.log')
    conda: '../envs/biscuit.yaml'
    shell:
        'trim_galore --paired {input} --output_dir {params.outdir} 2> {log}'


# BISCUIT accepts one FASTQ pair; keep resequenced runs in identical mate order.
rule merge_trimmed_reads:
    input:
        read1 = lambda wildcards: expand(rules.trim_adaptors.output.read1, sample=wildcards.sample, fqid=sample_to_fqid[wildcards.sample]),
        read2 = lambda wildcards: expand(rules.trim_adaptors.output.read2, sample=wildcards.sample, fqid=sample_to_fqid[wildcards.sample])
    output:
        read1 = temp(join(config['outdir'], 'biscuit/{sample}/merged_R1.fq.gz')),
        read2 = temp(join(config['outdir'], 'biscuit/{sample}/merged_R2.fq.gz'))
    shell:
        '''
        cat {input.read1:q} > {output.read1:q}
        cat {input.read2:q} > {output.read2:q}
        '''


rule align_to_ref:
    input:
        read1 = rules.merge_trimmed_reads.output.read1,
        read2 = rules.merge_trimmed_reads.output.read2,
        index = ancient(config['biscuit_index_inputs'])
    output:
        temp(join(config['outdir'], 'biscuit/{sample}/{sample}.biscuit_aligned.bam'))
    params:
        base = config['biscuit_index_prefix'],
        biscuit_args = config['biscuit_align_args']
    log: join(config['outdir'], 'biscuit/{sample}/{sample}.biscuit_align.log')
    threads: min(4, workflow.cores) 
    conda: '../envs/biscuit.yaml'
    shell:
        'biscuit align -@ {threads} {params.base:q} {params.biscuit_args} {input.read1:q} {input.read2:q} 2> {log:q} | samtools view -b -o {output:q}'


rule deduplicate_and_sort_bam:
    input:
        bam = rules.align_to_ref.output,
        ref_genome = ancient(rules.prepare_biscuit_reference.output.genome),
        fai = ancient(rules.prepare_biscuit_reference.output.fai)
    output:
        bam = join(config['outdir'], 'biscuit/{sample}/{sample}.dedup_sorted.bam'),
        stats = join(config['outdir'], 'biscuit/{sample}/{sample}.dupsifter.stat')
    log: join(config['outdir'], 'biscuit/{sample}/{sample}.dupsifter.log')
    conda: '../envs/biscuit.yaml'
    shell:
        'dupsifter {input.ref_genome} {input.bam} -O {output.stats} 2> {log} | samtools sort -o {output.bam}'


# merging is needed when one sample has multiple FASTQ IDs (e.g. when resequenced)
# rule merge_sort_aligned_reads:
#     input:
#         lambda wildcards: [join(config['outdir'], f'biscuit/{{sample}}/{fqid}_R1R2.dedup.bam') for fqid in metadata['DNA FASTQ ID'][wildcards.sample].split(' ')]
#     output:
#         join(config['outdir'], 'biscuit/{sample}/{sample}_DNA.dedup.sorted.bam')
#     conda: '../envs/biscuit.yaml'
#     shell:
#         'samtools cat {input} | samtools sort -o {output}'


rule index_bam:
    input:
        join(config['outdir'], 'biscuit/{sample}/{sample}.dedup_sorted.bam')
    output:
        join(config['outdir'], 'biscuit/{sample}/{sample}.dedup_sorted.bam.bai')
    conda: '../envs/biscuit.yaml'
    shell:
        'samtools index {input}'


rule extract_variants:
    input:
        bam = rules.deduplicate_and_sort_bam.output.bam,
        bai = rules.deduplicate_and_sort_bam.output.bam + '.bai',
        ref_genome = ancient(rules.prepare_biscuit_reference.output.genome),
        fai = ancient(rules.prepare_biscuit_reference.output.fai)
    output:
        vcf = temp(join(config['outdir'], 'biscuit/{sample}/{sample}_variants.vcf.bgz')),
        meth_average = join(config['outdir'], 'biscuit/{sample}/{sample}_meth_average.tsv')
    params:
        biscuit_args = config['biscuit_pileup_args'],
        stats_prefix = join(config['outdir'], 'biscuit/{sample}/{sample}')
    log: join(config['outdir'], 'biscuit/{sample}/{sample}.biscuit_pileup.log')
    threads: min(4, workflow.cores)
    conda: '../envs/biscuit.yaml'
    shell:
        'biscuit pileup -@ {threads} -N -w {params.stats_prefix:q} {params.biscuit_args} {input.ref_genome} {input.bam} 2> {log} | bgzip -@ {threads} -o {output.vcf}'


rule extract_methylation:
    input:
        rules.extract_variants.output.vcf
    output:
        hcg = join(config['outdir'], 'biscuit/{sample}/{sample}_HCG.bed'),
        gch = join(config['outdir'], 'biscuit/{sample}/{sample}_GCH.bed')
    log: join(config['outdir'], 'biscuit/{sample}/{sample}.biscuit_vcf2bed.log')
    conda: '../envs/biscuit.yaml'
    shell:
        '''
        biscuit vcf2bed -k 1 -t hcg {input} > {output.hcg} 2>> {log}
        biscuit vcf2bed -k 1 -t gch {input} > {output.gch} 2>> {log}
        '''


rule prepare_methscan_data_and_rename_columns:
    input:
        expand(join(config['outdir'], 'biscuit/{sample}/{sample}_HCG.bed'), sample=sample_to_fqid.keys())
    output:
        directory(join(config['outdir'], 'methscan/compact_data'))
    params:
        methscan_args = config['methscan_prepare_args']
    log: join(config['outdir'], 'methscan/methscan_prepare.log')
    conda: '../envs/methscan.yaml'
    run:
        shell('methscan prepare {params.methscan_args} {input} {output} 2> {log}')
        # rename the columns of methscan data to the original sample names
        cell_stats = pd.read_csv(join(str(output), 'cell_stats.csv'))
        cell_stats['cell_name'] = cell_stats['cell_name'].str.removesuffix('_HCG')
        cell_stats['cell_name'].to_csv(join(str(output), 'column_header.txt'), header=False, index=False)
        cell_stats.to_csv(join(str(output), 'cell_stats.csv'), index=False)


rule filter_methscan_data:
    input:
        data_dir = rules.prepare_methscan_data_and_rename_columns.output,
        cell_names = config.get('methscan_filter_cell_names', [])
    output:
        directory(join(config['outdir'], 'methscan/filtered_data'))
    params:
        methscan_args = config['methscan_filter_args']
    log: join(config['outdir'], 'methscan/methscan_filter.log')
    conda: '../envs/methscan.yaml'
    shell:
        'methscan filter {params.methscan_args} {input.data_dir} {output} 2> {log}'


rule find_methscan_VMRs:
    input:
        rules.filter_methscan_data.output
    output:
        join(config['outdir'], f'methscan/{config["dataset"]}_VMRs.bed')
    params:
        smooth_args = config['methscan_smooth_args'],
        scan_args = config['methscan_scan_args']
    log: join(config['outdir'], 'methscan/methscan_scan.log')
    conda: '../envs/methscan.yaml'
    threads: workflow.cores
    shell:
        '''
        methscan smooth {params.smooth_args} {input} 2>> {log}
        methscan scan {params.scan_args} --threads {threads} {input} {output} 2>> {log}
        '''


rule construct_methscan_matrix:
    input:
        vmrs = rules.find_methscan_VMRs.output,
        data_dir = rules.filter_methscan_data.output
    output:
        directory(join(config['outdir'], 'methscan/VMR_matrices'))
    params:
        methscan_args = config['methscan_matrix_args']
    log: join(config['outdir'], 'methscan/methscan_matrix.log')
    conda: '../envs/methscan.yaml'
    threads: workflow.cores
    shell:
        'methscan matrix {params.methscan_args} --threads {threads} {input.vmrs} {input.data_dir} {output} 2> {log}'


rule build_meth_anndata:
    input:
        rules.construct_methscan_matrix.output
    output:
        join(config['outdir'], f'{config["dataset"]}.biscuit_methscan.h5ad')
    conda: '../envs/anndata.yaml'
    shell:
        'python3 scripts/summarize_methscan_matrix.py {input}/mean_shrunken_residuals.csv.gz -o {output}'
