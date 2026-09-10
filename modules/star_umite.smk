rule fix_gtf_exon_ids:
    input:
        gtf = ancient(config['reference']['genes'])
    output:
        gtf = join(config['outdir'], 'reference/genes_with_exon_id.gtf')
    run:
        import os

        os.makedirs(os.path.dirname(output.gtf), exist_ok=True)

        processed_lines = []

        with open(input.gtf) as f:
            lines = f.readlines()

        for line in lines:

            if line.startswith("#"):
                processed_lines.append(line)
                continue

            fields = line.strip().split("\t")

            if len(fields) < 9:
                processed_lines.append(line)
                continue

            if fields[2] != "exon":
                processed_lines.append(line)
                continue

            attrs = fields[8]
            attr_dict = {}

            for item in attrs.split(";"):
                item = item.strip()
                if not item:
                    continue

                parts = item.split(" ", 1)
                if len(parts) == 2:
                    key = parts[0]
                    value = parts[1].strip('"')
                    attr_dict[key] = value

            if "exon_id" not in attr_dict:
                if "transcript_id" in attr_dict and "exon_number" in attr_dict:
                    attr_dict["exon_id"] = (
                        attr_dict["transcript_id"] + "_" + attr_dict["exon_number"]
                    )

            new_attrs = "; ".join(
                f'{k} "{v}"' for k, v in attr_dict.items()
            ) + ";"

            fields[8] = new_attrs
            processed_lines.append("\t".join(fields) + "\n")

        with open(output.gtf, "w") as f:
            f.writelines(processed_lines)


# Index generation lives in the top-level Snakefile: one index per reference
# genome is shared by every batch, so it must not be declared per batch.


# rule link_fastq_files:
#     # so that umiextract produces the desired filename
#     input:
#         ancient(lambda wildcards: multiext(join(ilse_to_fastqdir[wildcards.seqid], '{fqid}/fastq/{fqid}'), '_R1.fastq.gz', '_R2.fastq.gz'))
#     output:
#         temp(multiext(join(config['outdir'], 'alignments/{sample}/{seqid}_{fqid}'), '_R1.fastq.gz', '_R2.fastq.gz'))
#     shellz:
#         '''
#         ln -s {input[0]} {output[0]}
#         ln -s {input[1]} {output[1]}
#         '''


# rule merge_fastq:
#     input:
#         read1 = ancient(lambda wildcards: sample_to_read1[wildcards.sample]),
#         read2 = ancient(lambda wildcards: sample_to_read2[wildcards.sample])
#     output:
#         read1_merged = temp(join(config['outdir'], 'star_alignments/{sample}/{sample}_merged_R1.fastq.gz')),
#         read2_merged = temp(join(config['outdir'], 'star_alignments/{sample}/{sample}_merged_R2.fastq.gz'))
#     shell:
#         r'''
#         cat {input.read1} > {output.read1_merged}
#         cat {input.read2} > {output.read2_merged}
#         '''


rule extract_umi:
    input:
        read1 = lambda wildcards: fqid_to_reads[wildcards.fqid]['read1'],
        read2 = lambda wildcards: fqid_to_reads[wildcards.fqid]['read2']
    output:
        read1_umi = temp(join(config['outdir'], 'star_alignments/{sample}/{fqid}_R1_umiextract.fastq.gz')),
        read2_umi = temp(join(config['outdir'], 'star_alignments/{sample}/{fqid}_R2_umiextract.fastq.gz'))
        #temp(multiext(join(config['outdir'], 'alignments/{sample}/{seqid}_{fqid}'), '_R1_umiextract.fastq.gz', '_R2_umiextract.fastq.gz'))
    params:
        outdir = join(config['outdir'], 'star_alignments/{sample}'),
        umiextract_args = config['umiextract_args']
    # umiextract parallelises over read pairs, and this rule passes exactly one,
    # so extra cores would be reserved but never used.
    threads: 1
    conda: '../envs/umite.yaml'
    shell:
        r'''
        umiextract \
            {params.umiextract_args} \
            -c {threads} \
            -1 {input.read1} \
            -2 {input.read2} \
            -d {params.outdir}
        '''


rule align_to_ref:
    input:
        read1 = lambda wildcards: expand(rules.extract_umi.output.read1_umi, sample=wildcards.sample, fqid=sample_to_fqid[wildcards.sample]),
        read2 = lambda wildcards: expand(rules.extract_umi.output.read2_umi, sample=wildcards.sample, fqid=sample_to_fqid[wildcards.sample]),
        indices = ancient(config['star_index_inputs'])
    output:
        join(config['outdir'], 'star_alignments/{sample}/{sample}_Aligned.out.bam')
    params:
        read1_comma = lambda wildcards, input: ','.join(input.read1),
        read2_comma = lambda wildcards, input: ','.join(input.read2),
        outprefix = join(config['outdir'], 'star_alignments/{sample}/{sample}_'),
        star_args = config['star_align_args'],
        index_dir = config['star_index_dir']
    log: join(config['outdir'], 'star_alignments/{sample}/STAR_alignment.log')
    threads: min(4, workflow.cores)
    conda: '../envs/star.yaml'
    shell:
        r'''
        STAR \
            {params.star_args} \
            --runThreadN {threads} \
            --genomeDir {params.index_dir:q} \
            --genomeLoad LoadAndKeep \
            --readFilesIn {params.read1_comma} {params.read2_comma} \
            --outFileNamePrefix {params.outprefix} \
            > {log}
        '''


rule sort_bam_by_query_name:
    input:
        rules.align_to_ref.output
    output:
        temp(join(config['outdir'], 'alignments/{sample}/{sample}_Aligned.qn_sorted.bam'))
    conda: '../envs/star.yaml'
    shell:
        'samtools cat {input} | samtools sort -n -o {output}' # umicount requires the BAM file to be sorted by query name (instead of genomic location)


# umicount stores the parsed features per strand only when reads will query
# them by strand, and refuses a dump whose mode differs from --stranded. The
# mode is therefore part of the dump name so the two rules cannot disagree.
stranded = config['umicount_stranded']


rule parse_dump_GTF:
    input:
        temp(rules.fix_gtf_exon_ids.output)
    output:
        join(config['outdir'], f'reference/umicount_GTF_dump.{stranded}.pkl')
    params:
        stranded = stranded
    conda: '../envs/umite.yaml'
    shell:
        'umicount -g {input} --GTF_dump {output} --stranded {params.stranded}'


rule count_umis:
    input:
        bams = expand(rules.sort_bam_by_query_name.output, sample=sample_to_fqid.keys()),
        gtf_dump = rules.parse_dump_GTF.output
    output:
        multiext(join(config['outdir'], 'umicount/umite'), '.D.tsv', '.R.tsv', '.U.tsv')
    params:
        outdir = join(config['outdir'], 'umicount'),
        umicount_args = config['umicount_args']
    log: join(config['outdir'], 'umicount/umicount.log')
    threads: min(4, workflow.cores)
    conda: '../envs/umite.yaml'
    shell:
        r'''
        umicount \
            {params.umicount_args} \
            --combine_unspliced \
            -c {threads} \
            --GTF_skip_parse {input.gtf_dump} \
            --bams {input.bams} \
            -d {params.outdir} \
            -l {log}
        '''


# rule rename_umicount_files:
#     input:
#         rules.count_umis.output
#     output:
#         multiext(join(config['outdir'], f'umicount/{config['dataset']}_umite'), '.D.tsv', '.RE.tsv', '.RI.tsv', '.UE.tsv', '.UI.tsv')
#     params:
#         filename_prefix = config['dataset'],
#         samplename_suffix = '_Aligned.qn_sorted.bam'
#     shell:
#         'python3 scripts/rename_umicount_output.py --filename_prefix {params.filename_prefix} --samplename_suffix {params.samplename_suffix} {input}'


rule build_trsc_anndata:
    input:
        join(config['outdir'], 'umicount/umite.U.tsv'),
        gtf_dump = rules.parse_dump_GTF.output
    output:
        join(config['outdir'], f'{config["dataset"]}.star_umite.h5ad')
    params:
        filename_prefix = config['dataset'],
        samplename_suffix = '_Aligned.qn_sorted.bam'
    conda: '../envs/anndata.yaml'
    shell:
        r'''
        python3 scripts/summarize_umicount_tsv.py \
            --filename_prefix {params.filename_prefix} \
            --samplename_suffix {params.samplename_suffix} \
            -g {input.gtf_dump} \
            {input[0]} \
            -o {output}
        '''