from os.path import join

import pandas as pd


sample_to_fqid = config["sample_to_fqid"]
fqid_to_reads = config["fqid_to_reads"]

wildcard_constraints:
    sample = r"[\w.-]+",
    fqid = r"[\w.-]+",
    region = r"[\w.-]+",
    mark = r"CpG|GpC"

include: f"{config['pipeline']}.smk"

# QC rules are additive and off by default. Rules defined inside a conditional
# include are still re-exported by the top-level "use rule * from <module>".
if config.get("generate_QC_plots"):
    include: f"qc_{config['pipeline']}.smk"
