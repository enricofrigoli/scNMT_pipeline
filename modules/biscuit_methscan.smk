# rule all:
#     input:
#         #join(config['outdir'], f'biscuit/p1_A1/p1_A1_HCG.bed') # test pipeline for one sample
#         #[join(config['outidr'], f'biscuit/{sample}/{sample}_HCG.bed') for sn in metadata.index]
#         join(config['outdir'], 'methscan/VMR_matrices'),
#         multiext(join(config['outdir'], f'qc_plots/{config['dataset']}'), '_CpG_count_vs_meth_frac.png')


rule prepare_biscuit_index:
    input:
        ancient(config['reference']['genome'])
    output:
        directory(join(dirname(config['reference']['genome']), 'biscuit_index'))
    params:
        prefix = join(dirname(config['reference']['genome']), 'biscuit_index', basename(config['reference']['genome']).removesuffix('.gz').removesuffix('.bgz')),
        alg = config['biscuit_index_alg']
    conda: '../envs/biscuit.yaml'
    shell:
        '''
        mkdir -p {output:q}
        biscuit index {input:q} -p {params.prefix:q} {params.alg}
        '''


rule trim_adaptors:
    input:
        ancient(lambda wildcards: multiext(join(fqid_to_dir[wildcards.fqid], '{fqid}/fastq/{fqid}'), '_R1.fastq.gz', '_R2.fastq.gz'))
    output:
        temp(multiext(join(config['outdir'], 'biscuit/{sample}/{fqid}'), '_R1_val_1.fq.gz', '_R2_val_2.fq.gz'))
    params:
        outdir = join(config['outdir'], 'biscuit/{sample}')
    log: join(config['outdir'], 'biscuit/{sample}/{fqid}.trim_galore.log')
    conda: '../envs/biscuit.yaml'
    shell:
        'trim_galore --paired {input} --output_dir {params.outdir} 2> {log}'


rule align_to_ref:
    input:
        trimmed_reads = lambda wildcards: expand(rules.trim_adaptors.output, sample=wildcards.sample, fqid=sample_to_fqid[wildcards.sample]),
        index = ancient(rules.prepare_biscuit_index.output)
    output:
        temp(join(config['outdir'], 'biscuit/{sample}/{sample}.biscuit_aligned.bam'))
    params:
        base = rules.prepare_biscuit_index.params.prefix,
        biscuit_args = config['biscuit_align_args']
    log: join(config['outdir'], 'biscuit/{sample}/{sample}.biscuit_align.log')
    threads: min(4, workflow.cores) 
    conda: '../envs/biscuit.yaml'
    shell:
        'biscuit align -@ {threads} {params.base} {params.biscuit_args} {input.trimmed_reads} 2> {log} | samtools view -b -o {output}'


rule deduplicate_and_sort_bam:
    input:
        bam = rules.align_to_ref.output,
        ref_genome = ancient(config['reference']['genome'])
    output:
        join(config['outdir'], 'biscuit/{sample}/{sample}.dedup_sorted.bam')
    log: join(config['outdir'], 'biscuit/{sample}/{sample}.dupsifter.stat')
    conda: '../envs/biscuit.yaml'
    shell:
        'dupsifter {input.ref_genome} {input.bam} -O {log} | samtools sort -o {output}'


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
        '{base}.bam'
    output:
        '{base}.bam.bai'
    conda: '../envs/biscuit.yaml'
    shell:
        'samtools index {input}'


rule extract_variants:
    input:
        bam = rules.deduplicate_and_sort_bam.output,
        bai = rules.deduplicate_and_sort_bam.output[0] + '.bai',
        ref_genome = ancient(config['reference']['genome'])
    output:
        temp(join(config['outdir'], 'biscuit/{sample}/{sample}_variants.vcf.bgz'))
    log: join(config['outdir'], 'biscuit/{sample}/{sample}.biscuit_pileup.log')
    threads: min(4, workflow.cores)
    conda: '../envs/biscuit.yaml'
    shell:
        'biscuit pileup -@ {threads} -N {input.ref_genome} {input.bam} 2> {log} | bgzip -@ {threads} -o {output}'


rule extract_methylation:
    input:
        rules.extract_variants.output
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
        input_format = 'biscuit_short'
    log: join(config['outdir'], 'methscan/methscan_prepare.log')
    conda: '../envs/methscan.yaml'
    run:
        shell('methscan prepare --input-format {params.input_format} {input} {output} 2> {log}')
        # rename the columns of methscan data to the original sample names
        cell_stats = pd.read_csv(join(str(output), 'cell_stats.csv'))
        cell_stats['cell_name'] = cell_stats['cell_name'].str.removesuffix('_HCG')
        cell_stats['cell_name'].to_csv(join(str(output), 'column_header.txt'), header=False, index=False)
        cell_stats.to_csv(join(str(output), 'cell_stats.csv'), index=False)


rule filter_methscan_data:
    input:
        rules.prepare_methscan_data_and_rename_columns.output
    output:
        directory(join(config['outdir'], 'methscan/filtered_data'))
    log: join(config['outdir'], 'methscan/methscan_filter.log')
    conda: '../envs/methscan.yaml'
    shell:
        'methscan filter --min-sites 10000 --min-meth 10 --max-meth 90 {input} {output} 2> {log}'


rule find_methscan_VMRs:
    input:
        rules.filter_methscan_data.output
    output:
        join(config['outdir'], f'methscan/{config['dataset']}_VMRs.bed')
    log: join(config['outdir'], 'methscan/methscan_scan.log')
    conda: '../envs/methscan.yaml'
    threads: workflow.cores
    shell:
        '''
        methscan smooth {input} 2>> {log}
        methscan scan --threads {threads} {input} {output} 2>> {log}
        '''


rule construct_methscan_matrix:
    input:
        vmrs = rules.find_methscan_VMRs.output,
        data_dir = rules.filter_methscan_data.output
    output:
        directory(join(config['outdir'], 'methscan/VMR_matrices'))
    log: join(config['outdir'], 'methscan/methscan_matrix.log')
    conda: '../envs/methscan.yaml'
    threads: workflow.cores
    shell:
        'methscan matrix --threads {threads} {input.vmrs} {input.data_dir} {output} 2> {log}'


rule build_meth_anndata:
    input:
        rules.construct_methscan_matrix.output
    output:
        join(config['outdir'], f'{config['dataset']}.biscuit_methscan.h5ad')
    conda: '../envs/anndata.yaml'
    shell:
        'python3 scripts/summarize_methscan_matrix.py {input}/mean_shrunken_residuals.csv.gz -o {output}'


# rule generate_methscan_qc_plots:
#     input:
#         methscan_data = rules.prepare_methscan_data_and_rename_columns.output,
#         tss = ancient(config['ref_genome']['tss'])
#     output:
#         multiext(join(config['outdir'], f'qc_plots/{config['dataset']}'), '_CpG_count_vs_meth_frac.png')
#     params:
#         outprefix = join(config['outdir'], f'qc_plots/{config['dataset']}')
#     conda: '../envs/plotting.yaml'
#     shell:
#         'python3 scripts/plot_methscan_qc.py -t {input.tss} {input.methscan_data} -o {params.outprefix}'