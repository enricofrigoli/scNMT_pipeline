from os.path import join

import pandas as pd


sample_to_fqid = config["sample_to_fqid"]
fqid_to_reads = config["fqid_to_reads"]

wildcard_constraints:
    sample = r"[\w.-]+",
    fqid = r"[\w.-]+"

include: f"{config['pipeline']}.smk"
