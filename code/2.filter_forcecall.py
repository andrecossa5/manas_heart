"""
Very lenient filtering on forcecall results.
"""

import os
import pandas as pd


##


# Paths
path_main = '/Users/cossa/Desktop/projects/manas_heart'
path_data = os.path.join(path_main, 'data')
path_input = os.path.join(path_main, 'data/input')
path_filtered = os.path.join(path_main, 'results')


##


# Read forcecall results
df = pd.read_csv(os.path.join(path_filtered, 'ALLELIC_TABLE.tsv.gz'), sep='\t')

# Reshape and wrangle
df['mutation_id'] = df['CHROM'].astype(str) + '_' + df['POS'].astype(str) + '_' + df['REF'] + '_' + df['ALT']
df['region'] = df['chunk'].str.split('.').str[0]


##


# Annotate forcecalled mutations
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


##


# Filter annotated muts
annot = annot.loc[
    (annot['n_heart']>=2) & \
    (annot['AF_ratio']>=.3) & \
    (annot['sum_AD']>=5) 
]
muts = annot.index.tolist()
df = df.query('mutation_id in @muts').copy()
df['mutation_id'].nunique()

# Write 
df.to_csv(os.path.join(path_filtered, 'ALLELIC_TABLE_FILTERED.tsv.gz'), sep='\t', index=False)


##