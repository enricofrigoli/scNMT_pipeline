import argparse
import os
import pickle

import pandas as pd
import anndata as ad


# UMITE emits exactly these read-assignment counters alongside the genes; see
# ReadCategory in umite.umicount (_unique is tracked but never written out).
READ_CATEGORIES = ['_unmapped', '_multimapping', '_no_feature', '_ambiguous']

# Must track GTF_DUMP_VERSION in umite.umicount, which bumps whenever the
# layout of the pickled GTF data changes.
GTF_DUMP_VERSION = 2


parser = argparse.ArgumentParser()
parser.add_argument('umicount_tsv')
parser.add_argument('--filename_prefix', type=str, required=True)
parser.add_argument('--samplename_suffix', type=str, required=True)
parser.add_argument('-g', '--dumped_gtf', required=False)
parser.add_argument('-o', '--output', required=True)
args = parser.parse_args()




df_umicount = pd.read_csv(args.umicount_tsv, sep='\t', index_col=0)

# The counters are per-cell QC rather than expression: keep them in .obs and
# leave only genes in the matrix. Selecting them first fails loudly if UMITE
# ever renames a category, instead of silently counting it as a gene.
missing = [column for column in READ_CATEGORIES if column not in df_umicount.columns]
if missing:
    raise ValueError(
        f"{args.umicount_tsv} is missing UMITE read categories: {', '.join(missing)}. "
        f"Expected {', '.join(READ_CATEGORIES)} before the gene columns."
    )
read_categories = df_umicount[READ_CATEGORIES].copy()
df_umicount = df_umicount.drop(columns=READ_CATEGORIES)

# rename rows (samples) and columns (genes)
sample_names = [
    os.path.basename(row_name).removesuffix(args.samplename_suffix)
    for row_name in df_umicount.index
]
df_umicount.index = sample_names
read_categories.index = sample_names

# umicount wraps the parsed GTF tuple in a versioned envelope, so that a dump
# cannot be read back under a strand mode it was not parsed with. Only the gene
# attributes are needed here, and those are identical in either mode.
gene_id_to_name = {}
with open(args.dumped_gtf, 'rb') as f:
    dumped_gtf = pickle.load(f)
if not isinstance(dumped_gtf, dict) or dumped_gtf.get('version') != GTF_DUMP_VERSION:
    raise ValueError(
        f"{args.dumped_gtf} is not a version {GTF_DUMP_VERSION} UMITE GTF dump. "
        "Delete it and rerun so parse_dump_GTF writes it with the installed umite."
    )
for gene_id, gene_names in dumped_gtf['gtf_data'][2].items():
    if len(gene_names) == 0 or gene_names[0] == '':
        gene_id_to_name[gene_id] = gene_id
    else:
        gene_id_to_name[gene_id] = gene_names[0]

df_umicount.rename(gene_id_to_name, axis=1, inplace=True)

# sum up columns with the same gene name
df_umicount = df_umicount.T.groupby(level=0).sum().T

# write AnnData
adata = ad.AnnData(df_umicount, obs=read_categories)
adata.write(args.output, compression='gzip')
