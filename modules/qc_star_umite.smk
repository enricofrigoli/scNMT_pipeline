# cDNA QC. Included from modality.smk only when generate_QC_plots is set.
# Every rule here reads files the core pipeline already produces, so enabling
# QC adds work but never changes an analysis result.

qc_dir = join(config['outdir'], 'qc')


rule qc_star_metrics:
    input:
        logs = expand(
            rules.align_to_ref.output.final_log, sample=sample_to_fqid.keys()
        ),
        script = workflow.source_path('../scripts/summarize_star_logs.py')
    output:
        join(qc_dir, 'tables/star_metrics.tsv')
    threads: 1
    resources:
        mem_mb=4000,
        walltime=60
    conda: '../envs/plotting.yaml'
    shell:
        'python3 {input.script:q} {input.logs:q} -o {output:q}'


# umicount reports its per-cell funnel into a log rather than a table. The rule
# depends on the matrices, which Snakemake tracks, and reads the log beside them.
rule qc_umicount_metrics:
    input:
        matrices = rules.count_umis.output,
        script = workflow.source_path('../scripts/summarize_umicount_log.py')
    output:
        join(qc_dir, 'tables/umicount_metrics.tsv')
    params:
        umicount_log = join(config['outdir'], 'umicount/umicount.log')
    threads: 1
    resources:
        mem_mb=4000,
        walltime=60
    conda: '../envs/plotting.yaml'
    shell:
        'python3 {input.script:q} {params.umicount_log:q} -o {output:q}'


rule qc_per_cell_metrics:
    input:
        tables = [
            rules.qc_star_metrics.output,
            rules.qc_umicount_metrics.output,
        ],
        script = workflow.source_path('../scripts/collect_per_cell_metrics.py')
    output:
        join(qc_dir, f'tables/{config["dataset"]}.cDNA.per_cell_metrics_mqc.tsv')
    params:
        section_id = 'cdna_per_cell_metrics',
        section_name = 'cDNA per-cell metrics'
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


rule qc_multiqc:
    input:
        star_logs = expand(
            rules.align_to_ref.output.final_log, sample=sample_to_fqid.keys()
        ),
        metrics = rules.qc_per_cell_metrics.output,
        multiqc_config = workflow.source_path('../multiqc_config.yaml')
    output:
        report = join(qc_dir, 'multiqc/multiqc_report.html'),
        data = directory(join(qc_dir, 'multiqc/multiqc_report_data'))
    params:
        outdir = join(qc_dir, 'multiqc'),
        title = f'{config["dataset"]} cDNA'
    log: join(qc_dir, 'multiqc/multiqc.log')
    # Parses one STAR log per cell of the batch.
    threads: 1
    resources:
        mem_mb=8000,
        walltime=120
    conda: '../envs/multiqc.yaml'
    shell:
        r'''
        multiqc {input.star_logs:q} {input.metrics:q} \
            --outdir {params.outdir:q} \
            --filename multiqc_report \
            --config {input.multiqc_config:q} \
            --title {params.title:q} \
            --force --no-ansi \
            > {log:q} 2>&1
        '''
