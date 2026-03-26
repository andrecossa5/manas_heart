"""
Create the final table of somatic mutations.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import plotting_utils as plu
import matplotlib.pyplot as plt
matplotlib.use('macOSX')
plu.set_rcParams()


##


# Paths
path_main = '/Users/cossa/Desktop/projects/manas_heart'
path_data = os.path.join(path_main, 'data')
path_input = os.path.join(path_main, 'data/input')
path_filtered = os.path.join(path_main, 'results/filtered.1')


##


# Read all files

# Sample chunk mapping
df_samples = pd.read_csv(os.path.join(path_input, 'heart_samples.csv'))
df_samples = df_samples.rename(columns={'Chunk': 'chunk', 'Sample': 'Sample_ID'})

# Filtered mutations
df = pd.read_csv(os.path.join(path_filtered, 'ALL_FILTERED.tsv'), sep='\t')
df['mutation_id'] = df.apply(lambda x: f"{x['CHROM']}_{x['POS']}_{x['REF']}_{x['ALT']}", axis=1)

# Original muts
df_LCM = pd.read_csv(os.path.join(path_data, 'Heart_metadata.csv'))

# Filter on the same samples
samples_merged = df_samples['Sample_ID'].unique()
samples_LCM = df_LCM['Sample_ID'].unique()
common = list(set(samples_merged) & set(samples_LCM))
df_LCM = df_LCM.query('NV>=3 and Sample_ID in @common')

# Checks
assert df['chunk'].isin(df_samples['chunk']).all()

# Explore 
# df[['AD_placenta', 'AD_heart', 'DP_heart', 'DP_placenta', 'AF_placenta', 'AF_heart']].describe()
# df[['median_MQ', 'median_BQ', 'SB_pval', 'MPOS']].describe()
# df[['AD_placenta', 'AD_heart']].hist(bins=20)
# plt.show()
# ...

# PASS: final test
df['PASS'] = (df['AD_heart']>=3) & \
             (df['DP_heart']>=10) & \
             (df['DP_placenta']>=10) & \
             (df['SB_pval']>=0.01) & \
             (df['median_BQ']>=20) & \
             (df['AD_placenta']<2)
df_filtered = df[df['PASS']].copy()
df_filtered
df['PASS'].value_counts(normalize=True)

df['AF_heart'].describe()
df['AF_placenta'].hist()
df['AF_heart'].hist()
plt.show()

# Single LCM
total_merged_set = set(df['mutation_id'].unique())
filtered_merged_set = set(df_filtered['mutation_id'].unique())
LCM_set = set(df_LCM['mutation_id'].unique())
print(f"Fraction of total mutations in LCM: {len(total_merged_set & LCM_set) / len(total_merged_set):.4%}")
print(f"Fraction of filtered mutations in LCM: {len(filtered_merged_set & LCM_set) / len(filtered_merged_set):.4%}")
print(f"Fraction of LCM mutations in total: {len(total_merged_set & LCM_set) / len(LCM_set):.4%}")
print(f"Fraction of LCM mutations in filtered: {(len(filtered_merged_set & LCM_set) / len(LCM_set)):.4%}")


##


# Group at region level
grouped = df.groupby(['mutation_id','region'])
grouped_sums = (
    grouped[['AD_heart', 'AD_placenta', 'DP_heart', 'DP_placenta']]
    .sum()
)
grouped_means = (
    grouped[['SB_pval', 'median_BQ']]
    .mean()
)
grouped_merged = grouped_sums.join(grouped_means).reset_index()

grouped_merged['PASS'] = (grouped_merged['AD_heart']>=3) & \
                         (grouped_merged['DP_heart']>=10) & \
                         (grouped_merged['DP_placenta']>=10) & \
                         (grouped_merged['SB_pval']>=0.01) & \
                         (grouped_merged['median_BQ']>=30) & \
                         (grouped_merged['AD_placenta']==0)

grouped_merged = grouped_merged[grouped_merged['PASS']].copy()
grouped_merged['AF_heart'] = grouped_merged['AD_heart'] / ( grouped_merged['DP_heart'] + 10**-18 )
grouped_merged['AF_placenta'] = grouped_merged['AD_placenta'] / ( grouped_merged['DP_placenta'] + 10**-18 )

grouped_merged['AF_heart'].describe()
grouped_merged['AD_heart'].describe()
grouped_merged['DP_heart'].describe()

df_LCM.loc[df_LCM['VAF']!=0, 'VAF'].min()

grouped_merged_set = set(grouped_merged['mutation_id'].unique())
print(f"Fraction of LCM mutations in grouped_merged: {(len(grouped_merged_set & LCM_set) / len(LCM_set)):.4%}")


##