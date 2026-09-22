# gDNA QC. Included from modality.smk only when generate_QC_plots is set.

qc_dir = join(config['outdir'], 'qc')
# Names come from reference.profile_regions and are validated as identifiers, so
# they are safe as wildcards, path components and plot titles.
profile_regions = config['reference'].get('profile_regions') or {}
# CpG comes from the HCG beds the core pipeline already prepares; GpC needs its
# own compact dataset, built from the GCH beds that were otherwise unused.
profile_data_dirs = {
    'CpG': join(config['outdir'], 'methscan/compact_data'),
    'GpC': join(config['outdir'], 'methscan/compact_data_GpC'),
}


rule qc_flagstat:
    input:
        bam = rules.deduplicate_and_sort_bam.output.bam,
        bai = rules.index_bam.output
    output:
        temp(join(config['outdir'], 'biscuit/{sample}/{sample}.flagstat'))
    threads: 1
    resources:
        mem_mb=2000,
        walltime=60
    conda: '../envs/biscuit.yaml'
    shell:
        'samtools flagstat {input.bam:q} > {output:q}'


rule qc_flagstat_metrics:
    input:
        reports = expand(rules.qc_flagstat.output, sample=sample_to_fqid.keys()),
        script = workflow.source_path('../scripts/summarize_flagstat.py')
    output:
        join(qc_dir, 'tables/flagstat_metrics.tsv')
    threads: 1
    resources:
        mem_mb=4000,
        walltime=60
    conda: '../envs/plotting.yaml'
    shell:
        'python3 {input.script:q} {input.reports:q} -o {output:q}'


rule qc_dupsifter_metrics:
    input:
        stats = expand(
            rules.deduplicate_and_sort_bam.output.stats, sample=sample_to_fqid.keys()
        ),
        script = workflow.source_path('../scripts/summarize_dupsifter.py')
    output:
        join(qc_dir, 'tables/dupsifter_metrics.tsv')
    threads: 1
    resources:
        mem_mb=4000,
        walltime=60
    conda: '../envs/plotting.yaml'
    shell:
        'python3 {input.script:q} {input.stats:q} -o {output:q}'


# HCH is the only cytosine class that is neither CpG nor GpC, so it is the one
# conversion estimate the NOMe-seq GpC label cannot inflate.
rule qc_conversion_metrics:
    input:
        averages = expand(
            rules.extract_variants.output.meth_average, sample=sample_to_fqid.keys()
        ),
        script = workflow.source_path('../scripts/summarize_hch_conversion.py')
    output:
        join(qc_dir, 'tables/conversion_metrics.tsv')
    threads: 1
    resources:
        mem_mb=4000,
        walltime=60
    conda: '../envs/plotting.yaml'
    shell:
        'python3 {input.script:q} {input.averages:q} -o {output:q}'


rule qc_per_cell_metrics:
    input:
        tables = [
            rules.qc_conversion_metrics.output,
            rules.qc_dupsifter_metrics.output,
            rules.qc_flagstat_metrics.output,
        ],
        script = workflow.source_path('../scripts/collect_per_cell_metrics.py')
    output:
        join(qc_dir, f'tables/{config["dataset"]}.gDNA.per_cell_metrics_mqc.tsv')
    params:
        section_id = 'gdna_per_cell_metrics',
        section_name = 'gDNA per-cell metrics'
    threads: 1
    resources:
        mem_mb=4000,
        walltime=60
    conda: '../envs/plotting.yaml'
    shell:
        r'''
        python3 {input.script:q} {input.tables:q} \
            --section-id {params.section_id:q} \
            --section-name {params.section_name:q} \
            -o {output:q}
        '''


# Conversion, duplication and mapping rate are the metrics whose failures tend
# to follow plate position rather than the cell itself.
rule qc_plate_heatmap:
    input:
        metrics = rules.qc_per_cell_metrics.output,
        script = workflow.source_path('../scripts/plot_plate_heatmap.py')
    output:
        join(qc_dir, 'plate/plate_qc.png')
    params:
        title = f'{config["dataset"]} gDNA plate QC'
    threads: 1
    resources:
        mem_mb=4000,
        walltime=60
    conda: '../envs/plotting.yaml'
    shell:
        r'''
        python3 {input.script:q} {input.metrics:q} \
            -m conversion_rate_pct \
            -m duplicate_pct \
            -m mapped_pct \
            --title {params.title:q} \
            -o {output:q}
        '''


rule qc_methscan_cell_stats:
    input:
        data_dir = rules.prepare_methscan_data_and_rename_columns.output,
        script = workflow.source_path('../scripts/plot_methscan_cell_stats.py')
    output:
        join(qc_dir, 'cell_stats.png')
    threads: 1
    resources:
        mem_mb=4000,
        walltime=60
    conda: '../envs/plotting.yaml'
    shell:
        'python3 {input.script:q} {input.data_dir:q}/cell_stats.csv -o {output:q}'


# The GCH beds are already written per cell but nothing consumed them. Profiling
# needs its own compact dataset; it is not filtered, because the point of these
# plots is to show the cells a filter would have removed.
rule prepare_methscan_data_GpC:
    input:
        expand(
            join(config['outdir'], 'biscuit/{sample}/{sample}_GCH.bed'),
            sample=sample_to_fqid.keys()
        )
    output:
        directory(join(config['outdir'], 'methscan/compact_data_GpC'))
    params:
        methscan_args = config['methscan_prepare_args']
    log: join(config['outdir'], 'methscan/methscan_prepare_GpC.log')
    # As costly as the CpG preparation: one job over every cell's GCH bed.
    threads: 1
    resources:
        mem_mb=32000,
        walltime=720
    conda: '../envs/methscan.yaml'
    run:
        shell('methscan prepare {params.methscan_args} {input} {output} 2> {log}')
        # methscan profile reads cell names from column_header.txt, so strip the
        # mark suffix that the per-cell BED filenames introduced.
        cell_stats = pd.read_csv(join(str(output), 'cell_stats.csv'))
        cell_stats['cell_name'] = cell_stats['cell_name'].str.removesuffix('_GCH')
        cell_stats['cell_name'].to_csv(
            join(str(output), 'column_header.txt'), header=False, index=False
        )
        cell_stats.to_csv(join(str(output), 'cell_stats.csv'), index=False)


rule qc_methscan_profile:
    input:
        data_dir = lambda wildcards: profile_data_dirs[wildcards.mark],
        regions = lambda wildcards: profile_regions[wildcards.region]
    output:
        join(qc_dir, 'profiles/{region}_{mark}.csv')
    params:
        methscan_args = config['methscan_profile_args']
    log: join(qc_dir, 'profiles/{region}_{mark}.methscan_profile.log')
    threads: 1
    resources:
        mem_mb=16000,
        walltime=240
    conda: '../envs/methscan.yaml'
    shell:
        'methscan profile {params.methscan_args} {input.regions:q} {input.data_dir:q} {output:q} 2> {log:q}'


rule qc_plot_region_profiles:
    input:
        cpg = join(qc_dir, 'profiles/{region}_CpG.csv'),
        gpc = join(qc_dir, 'profiles/{region}_GpC.csv'),
        script = workflow.source_path('../scripts/plot_region_profiles.py')
    output:
        join(qc_dir, 'profiles/{region}_profiles.pdf')
    threads: 1
    resources:
        mem_mb=4000,
        walltime=60
    conda: '../envs/plotting.yaml'
    shell:
        r'''
        python3 {input.script:q} \
            --cpg {input.cpg:q} \
            --gpc {input.gpc:q} \
            --region {wildcards.region:q} \
            -o {output:q}
        '''


rule qc_multiqc:
    input:
        trimming_reports = expand(
            [rules.trim_adaptors.output.report1, rules.trim_adaptors.output.report2],
            zip,
            sample=[
                sample for sample, fqids in sample_to_fqid.items() for _ in fqids
            ],
            fqid=[fqid for fqids in sample_to_fqid.values() for fqid in fqids],
        ),
        metrics = rules.qc_per_cell_metrics.output,
        multiqc_config = workflow.source_path('../multiqc_config.yaml')
    output:
        report = join(qc_dir, 'multiqc/multiqc_report.html'),
        data = directory(join(qc_dir, 'multiqc/multiqc_report_data'))
    params:
        outdir = join(qc_dir, 'multiqc'),
        title = f'{config["dataset"]} gDNA'
    log: join(qc_dir, 'multiqc/multiqc.log')
    # Parses two trimming reports per cell of the batch.
    threads: 1
    resources:
        mem_mb=8000,
        walltime=120
    conda: '../envs/multiqc.yaml'
    shell:
        r'''
        multiqc {input.trimming_reports:q} {input.metrics:q} \
            --outdir {params.outdir:q} \
            --filename multiqc_report \
            --config {input.multiqc_config:q} \
            --title {params.title:q} \
            --force --no-ansi \
            > {log:q} 2>&1
        '''
