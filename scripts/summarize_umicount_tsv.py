import argparse
import os
import pickle

import pandas as pd
import anndata as ad




parser = argparse.ArgumentParser()
parser.add_argument('umicount_tsv')
parser.add_argument('--filename_prefix', type=str, required=True)
parser.add_argument('--samplename_suffix', type=str, required=True)
parser.add_argument('-g', '--dumped_gtf', required=False)
parser.add_argument('-o', '--output', required=True)
args = parser.parse_args()




df_umicount = pd.read_csv(args.umicount_tsv, sep='\t', index_col=0)

# UMITE writes assignment counters alongside genes; exclude them from expression.
df_umicount.drop(
    columns=['_unmapped', '_multimapping', '_no_feature', '_ambiguous'],
    errors='ignore',
    inplace=True,
)

# rename rows (samples) and columns (genes)
df_umicount.index = [os.path.basename(row_name).removesuffix(args.samplename_suffix) for row_name in df_umicount.index]

gene_id_to_name = {}
with open(args.dumped_gtf, 'rb') as f:
    for gene_id, gene_names in pickle.load(f)[2].items():
        if len(gene_names) == 0 or gene_names[0] == '':
            gene_id_to_name[gene_id] = gene_id
        else:
            gene_id_to_name[gene_id] = gene_names[0]

df_umicount.rename(gene_id_to_name, axis=1, inplace=True)

# sum up columns with the same gene name
df_umicount = df_umicount.T.groupby(level=0).sum().T

# write AnnData
adata = ad.AnnData(df_umicount)
adata.write(args.output, compression='gzip')
