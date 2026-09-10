from datetime import datetime
import hashlib
from pathlib import Path
import re

import yaml
from snakemake.utils import min_version

from workflow_config import (
    STAR_INDEX_FILES,
    normalize_run_config,
    prepare_batch_configs,
    prepare_modality_config,
)


min_version("9.19")

configfile: "snakeconfig.yaml"
with open(workflow.source_path("defaultconfig.yaml")) as handle:
    default_config = yaml.safe_load(handle)
config = normalize_run_config(default_config | config, workflow.basedir)
batch_configs = prepare_batch_configs(config, workflow.basedir)

layer_configs = {}
star_cleanup_groups = {}
star_index_builds = {}
all_outputs = []
for batch, batch_config in batch_configs.items():
    modalities = ("cDNA", "gDNA") if batch_config["modality"] == "both" else (batch_config["modality"],)
    for modality in modalities:
        try:
            layer = prepare_modality_config(batch_config, modality)
        except ValueError as exc:
            raise ValueError(f"{batch or batch_config['dataset']}/{modality}: {exc}") from exc
        layer_configs[(batch, modality)] = layer
        module_name = f"{batch}_{modality}_" if batch is not None else f"{modality}_"
        module:
            name: module_name
            snakefile: "modules/modality.smk"
            config: layer

        use rule * from module_name as module_name*

        all_outputs.append(str(Path(layer["outdir"]) / f"{layer['dataset']}.{layer['pipeline']}.h5ad"))
        if modality == "gDNA":
            all_outputs.append(str(Path(layer["outdir"]) / "qc_plots/cell_stats.png"))
        else:
            index_dir = str(Path(layer["star_index_dir"]).resolve())
            index_key = hashlib.sha256(index_dir.encode()).hexdigest()[:16]
            group = star_cleanup_groups.setdefault(index_key, {
                "index_dir": index_dir, "index_inputs": layer["star_index_inputs"], "bams": [],
            })
            group["bams"].extend(
                str(Path(layer["outdir"]) / f"star_alignments/{sample}/{sample}_Aligned.out.bam")
                for sample in layer["sample_to_fqid"]
            )
            if not layer["star_index_generated"]:
                continue
            # One generated index per genome, so every batch sharing it must agree.
            layer_name = f"{batch or layer['dataset']}/{modality}"
            settings = {
                "genome": layer["reference"]["genome"],
                "genes": layer["reference"]["genes"],
                "index_args": layer["star_index_args"],
            }
            build = star_index_builds.setdefault(
                index_dir, {**settings, "overhang": 0, "layers": []}
            )
            conflicting = sorted(key for key, value in settings.items() if build[key] != value)
            if conflicting:
                raise ValueError(
                    f"{index_dir} would be generated with different "
                    f"{', '.join(conflicting)} by {build['layers'][0]} and {layer_name}; "
                    "supply reference.star_index per batch to keep separate indexes"
                )
            # STAR recommends max(read length) - 1 across everything using the index.
            build["overhang"] = max(build["overhang"], layer["star_sjdb_overhang"])
            build["layers"].append(layer_name)


if star_index_builds:
    rule generate_star_index:
        input:
            ref_genome = lambda wildcards: ancient(star_index_builds[wildcards.index_dir]["genome"]),
            gene_annotation = lambda wildcards: ancient(star_index_builds[wildcards.index_dir]["genes"])
        output:
            multiext("{index_dir}/", *STAR_INDEX_FILES)
        params:
            star_args = lambda wildcards: star_index_builds[wildcards.index_dir]["index_args"],
            overhang = lambda wildcards: star_index_builds[wildcards.index_dir]["overhang"],
            outprefix = "{index_dir}/star_genome_generate_"
        log:
            # Kept beside the index: the log must carry the output's wildcards, and
            # two genomes generating indexes would otherwise share one log path.
            "{index_dir}/star_genome_generate.log"
        wildcard_constraints:
            index_dir = "|".join(re.escape(index_dir) for index_dir in star_index_builds)
        threads: workflow.cores
        conda: "envs/star.yaml"
        shell:
            r'''
            STAR \
                {params.star_args} \
                --runThreadN {threads} \
                --runMode genomeGenerate \
                --genomeDir {wildcards.index_dir:q} \
                --genomeFastaFiles {input.ref_genome:q} \
                --sjdbGTFfile {input.gene_annotation:q} \
                --sjdbOverhang {params.overhang} \
                --outFileNamePrefix {params.outprefix:q} \
                > {log:q} 2>&1
            '''


if star_cleanup_groups:
    # A shared index must remain loaded until every batch using it has aligned.
    # Counts do not depend on this flag: adding a batch must not recount old ones.
    rule unload_star_genome:
        input:
            bams = lambda wildcards: star_cleanup_groups[wildcards.index_key]["bams"],
            index_files = lambda wildcards: star_cleanup_groups[wildcards.index_key]["index_inputs"]
        output:
            touch(str(Path(config["outdir"]) / "reference/star_cleanup/{index_key}.flag"))
        params:
            index_dir = lambda wildcards: star_cleanup_groups[wildcards.index_key]["index_dir"],
            outprefix = str(Path(config["outdir"]) / "reference/star_cleanup/{index_key}_")
        log:
            str(Path(config["outdir"]) / "reference/star_cleanup/{index_key}.log")
        wildcard_constraints:
            index_key = r"[0-9a-f]{16}"
        conda: "envs/star.yaml"
        shell:
            r'''
            star_status=0
            STAR --genomeLoad Remove --genomeDir {params.index_dir:q} \
                --outFileNamePrefix {params.outprefix:q} > {log:q} 2>&1 || star_status=$?
            if [ "$star_status" -ne 0 ]; then
                if [ "$star_status" -eq 105 ] && grep -Fq \
                    'Did not find the genome in memory, did not remove any genomes from shared memory' {log:q}; then
                    :
                else
                    cat {log:q} >&2
                    exit "$star_status"
                fi
            fi
            '''

    all_outputs.extend(
        str(Path(config["outdir"]) / f"reference/star_cleanup/{index_key}.flag")
        for index_key in star_cleanup_groups
    )


rule all:
    default_target: True
    input:
        all_outputs
    run:
        with open(Path(config["outdir"]) / "last_success_config.yaml", "w") as handle:
            timestamp = datetime.today().isoformat(sep=" ", timespec="seconds")
            handle.write(f"# Configuration of the most recent successful run: {timestamp}\n")
            yaml.safe_dump(config, handle, sort_keys=False)
