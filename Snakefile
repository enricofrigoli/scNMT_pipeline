from datetime import datetime
import hashlib
from pathlib import Path

import yaml
from snakemake.utils import min_version

from workflow_config import normalize_run_config, prepare_batch_configs, prepare_modality_config


min_version("9.19")

configfile: "snakeconfig.yaml"
with open(workflow.source_path("defaultconfig.yaml")) as handle:
    default_config = yaml.safe_load(handle)
config = normalize_run_config(default_config | config, workflow.basedir)
batch_configs = prepare_batch_configs(config, workflow.basedir)

layer_configs = {}
star_cleanup_groups = {}
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
