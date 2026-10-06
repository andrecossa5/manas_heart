"""
Final callset: the positive SNVs that are not flagged as artefacts (results/ALLELIC_TABLE_FILTERED_ANNOTATED.tsv.gz),
with their position in the phylogeny of the donor (data/Filteredmutations_14061_Sample_subset_snv_assigned_to_branches.txt)
and a few annotations. Writes results/ALLELIC_TABLE_FINAL.tsv.gz.
"""

import os
import numpy as np
import pandas as pd


##


# Paths
path_main = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))     # repository root
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results')

df = pd.read_csv(os.path.join(path_results, 'ALLELIC_TABLE_FILTERED_ANNOTATED.tsv.gz'), sep='\t')
tree_assignment = pd.read_csv(
    os.path.join(path_data, 'Filteredmutations_14061_Sample_subset_snv_assigned_to_branches.txt'), sep='\t'
)


##


# Drop artefacts
df_filtered = df.query('AF>0 and artifact_flag == "No artefact"')

# Where were the SNVs assigned in the phylogeny? Not assigned: called in the merged BAMs only
tree_assignment['mutation_id'] = (
    tree_assignment['Chr'].astype(str) + '_' +
    tree_assignment['Pos'].astype(str) + '_' +
    tree_assignment['Ref'].astype(str) + '_' +
    tree_assignment['Alt'].astype(str)
)
df_filtered = (
    df_filtered.merge(
        tree_assignment[['mutation_id', 'desc_samples_orgin']],
        on='mutation_id', how='left'
    )
)
df_filtered.loc[df_filtered['desc_samples_orgin'].isna(), 'desc_samples_orgin'] = 'Unassigned'

# Annotate CpGs and calling strategy
CpGs = ['A[C>T]G', 'C[C>T]G', 'G[C>T]G', 'T[C>T]G']
df_filtered['calling_strategy'] = (
    np.where(df_filtered['desc_samples_orgin'] == 'Unassigned', 'Merged', 'LCM')
)
df_filtered['in_CpGs'] = (
    np.where(df_filtered['SBS96'].isin(CpGs), 'CpG', 'Non-CpG')
)

# Sensible: >=3 alt reads in a sample and no AF above 0.25
muts_too_high_AF = (
    df_filtered.query('tissue!="placenta"')
    .groupby('mutation_id')
    ['AF'].max().loc[lambda x: x>0.25].index
)
muts_sensible = (
    df_filtered.query('tissue!="placenta"')
    .groupby('mutation_id')
    ['AD_alt'].max().loc[lambda x: x>=3].index
)
MUTS = set(muts_sensible) - (set(muts_too_high_AF))
df_filtered['in_sensible'] = df_filtered['mutation_id'].isin(MUTS)

# Shared: >=5 alt reads in a chunk, in >=2 regions
MIN_AD_CHUNK = 5
MIN_CHUNKS = 2
muts_chunks = (
    df_filtered.query('tissue!="placenta"')
    .groupby(['region', 'chunk', 'mutation_id'])
    ['AD_alt'].sum().loc[lambda x: x>=MIN_AD_CHUNK]
    .reset_index()
    .groupby('mutation_id')['region'].nunique().loc[lambda x: x>=MIN_CHUNKS].index
)
MUTS = set(muts_chunks) - (set(muts_too_high_AF))
df_filtered['in_shared'] = df_filtered['mutation_id'].isin(MUTS)

print(f'ALLELIC_TABLE_FINAL: {df_filtered["mutation_id"].nunique()} SNVs')
df_filtered.to_csv(os.path.join(path_results, 'ALLELIC_TABLE_FINAL.tsv.gz'), sep='\t', index=False)
