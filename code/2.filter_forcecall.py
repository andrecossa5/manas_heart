"""
Very lenient filtering on the force-called table (results/ALLELIC_TABLE.tsv.gz, from the Nextflow
force_call workflow): a site is kept if it is called in >=2 heart LCM samples, has >=5 alt reads in total,
and the mean AF over its positive samples is >=30% above its mean AF over all heart samples (signal not
spread thinly over every sample). Writes results/ALLELIC_TABLE_FILTERED.tsv.gz.
"""

import os
import pandas as pd


##


# Paths
path_main = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))     # repository root
path_results = os.path.join(path_main, 'results')


##


# Read forcecall results
df = pd.read_csv(os.path.join(path_results, 'ALLELIC_TABLE.tsv.gz'), sep='\t')

# Reshape and wrangle
df['mutation_id'] = df['CHROM'].astype(str) + '_' + df['POS'].astype(str) + '_' + df['REF'] + '_' + df['ALT']
df['region'] = df['chunk'].str.split('.').str[0]

# Annotate forcecalled mutations
placenta = (
    df.query('tissue=="placenta"')
    .groupby('mutation_id')
    .apply(lambda x: pd.Series({
        'mean_AF_placenta': x['AF'].mean(),
        'mean_AF_pos_placenta': x.loc[x['AF']>0, 'AF'].mean(),
    }))
)
heart = (
    df.query('tissue=="heart"')
    .groupby('mutation_id')
    .apply(lambda x: pd.Series({
        'mean_AF_heart': x['AF'].mean(),
        'mean_AF_pos_heart': x.loc[x['AF']>0, 'AF'].mean(),
        'n_heart': x['Sample_ID'].nunique(),
        'sum_AD': x['AD_alt'].sum(),
        'mean_AD_pos': x.loc[x['AD_alt']>0, 'AD_alt'].mean()
    }))
)
annot = placenta.join(heart, how='outer')
annot.loc[annot['mean_AF_pos_heart'].isna()] = 0
annot['AF_ratio'] = (annot['mean_AF_pos_heart'] - annot['mean_AF_heart']) / (annot['mean_AF_heart']+10**(-18))

# Filter annotated muts
annot = annot.loc[
    (annot['n_heart']>=2) & \
    (annot['AF_ratio']>=.3) & \
    (annot['sum_AD']>=5)
]
muts = annot.index.tolist()
df = df.query('mutation_id in @muts').copy()
print(f'ALLELIC_TABLE_FILTERED: {df["mutation_id"].nunique()} SNVs')

# Write
df.to_csv(os.path.join(path_results, 'ALLELIC_TABLE_FILTERED.tsv.gz'), sep='\t', index=False)
