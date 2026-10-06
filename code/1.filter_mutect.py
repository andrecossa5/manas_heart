"""
Mutect2 calls on the merged BAMs (results/ALL_FILTERED_ctx.tsv.gz, from the Nextflow call_variants
workflow) -> hard filters -> results/FILTER.1.tsv.

Also writes the sites to force-call in the single LCM samples (results/forcecall.tsv.gz, input of the
Nextflow force_call workflow, which produces results/ALLELIC_TABLE.tsv.gz): every site seen with >=1 alt
read in a heart LCM sample of the callset (data/Heart_metadata.csv), plus the FILTER.1 sites.
"""

import os
import pandas as pd


##


# Paths
path_main = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))     # repository root
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results')


##


# Merged-BAM Mutect2 calls, with trinucleotide context
df_samples = pd.read_csv(os.path.join(path_data, 'input', 'heart_samples.csv'))
df = pd.read_csv(os.path.join(path_results, 'ALL_FILTERED_ctx.tsv.gz'), sep='\t')
df['mutation_id'] = df['CHROM'] + '_' + df['POS'].astype(str) + '_' + df['REF'] + '_' + df['ALT']
assert df['chunk'].isin(df_samples['chunk']).all()

# FILTER.1: SNVs only (a trinucleotide context is defined), and hard thresholds on the calls
df = df.dropna(subset=['SBS6'])
df['PASS'] = (df['SB_pval']>0.1) & \
             (df['median_BQ']>=30) & \
             (df['MPOS']>5) & \
             (df['AD_placenta']<2) & \
             (df['NLOD']>10) & \
             (df['TLOD']>25) & \
             (df['POPAF']>2) & \
             ((df['orientation_ratio']>0.3) & (df['orientation_ratio']<0.7))
df_filtered = df[df['PASS']].copy()
df_filtered.to_csv(os.path.join(path_results, 'FILTER.1.tsv'), sep='\t', index=False)
print(f'FILTER.1: {df_filtered["mutation_id"].nunique()} SNVs')


##


# Sites to force-call: single-LCM calls with >=1 alt read (in the samples of the merged BAMs) + FILTER.1
df_LCM = pd.read_csv(os.path.join(path_data, 'Heart_metadata.csv'), usecols=['mutation_id', 'Sample_ID', 'NV'])
samples_common = set(df_samples['Sample_ID']) & set(df_LCM['Sample_ID'])
df_LCM = df_LCM.query('NV>0 and Sample_ID in @samples_common')


def split_sites(mutation_id):
    sites = mutation_id.str.split('_', expand=True)
    sites.columns = ['CHROM', 'POS', 'REF', 'ALT']
    return sites.drop_duplicates().reset_index(drop=True)


df_forcecall = pd.concat([split_sites(df_LCM['mutation_id']), split_sites(df_filtered['mutation_id'])])
df_forcecall.reset_index(drop=True).to_csv(os.path.join(path_results, 'forcecall.tsv.gz'), sep='\t', index=False)
print(f'forcecall: {len(df_forcecall)} sites')
