"""
Analyze filtered calls.
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
path_filtered = os.path.join(path_main, 'results')
path_figures = os.path.join(path_main, 'figures')


##


# Read all files

# Sample chunk mapping
df_samples = pd.read_csv(os.path.join(path_input, 'heart_samples.csv'))
df_samples = df_samples.rename(columns={'Chunk': 'chunk', 'Sample': 'Sample_ID'})

# ALL FILTERED mutations
df = pd.read_csv(os.path.join(path_filtered, 'ALL_FILTERED.tsv.gz'), sep='\t')
df['mutation_id'] = df.apply(lambda x: f"{x['CHROM']}_{x['POS']}_{x['REF']}_{x['ALT']}", axis=1)

# ALL FILTERED stats
stats = pd.read_csv(os.path.join(path_filtered, 'ALL_FILTERED.stats'), sep='\t')
stats['mutation_id'] = df.apply(lambda x: f"{x['CHROM']}_{x['POS']}_{x['REF']}_{x['ALT']}", axis=1)

# Original muts
df_LCM = pd.read_csv(os.path.join(path_data, 'Heart_metadata.csv'))

# Filter on the same samples
samples_merged = df_samples['Sample_ID'].unique()
samples_LCM = df_LCM['Sample_ID'].unique()
common = list(set(samples_merged) & set(samples_LCM))
df_LCM = df_LCM.query('NV>=4 and Sample_ID in @common')

# Checks
assert df['chunk'].isin(df_samples['chunk']).all()


##


df['AD_placenta'].max()

# Explore 
# df[['AD_placenta', 'AD_heart', 'DP_heart', 'DP_placenta', 'AF_placenta', 'AF_heart']].describe()
# df[['median_MQ', 'median_BQ', 'SB_pval', 'MPOS']].describe()
# df[['AD_placenta', 'AD_heart']].hist(bins=20)
# plt.show()
# ...

df_ = stats.groupby('filter')['n_records'].sum().to_frame('n').reset_index()
order = df_.sort_values('n', ascending=False)['filter'].to_list()

fig, ax = plt.subplots(figsize=(3.5,3.5))
plu.bar(df_, x='filter', y='n', ax=ax, x_order=order)
plu.format_ax(ax, xlabel='', ylabel='n records', rotx=90, reduced_spines=True)
ax.set_yscale('log')
plt.tight_layout()
fig.savefig(os.path.join(path_figures, 'all_stats.pdf'))

##

fig, ax = plt.subplots(figsize=(3.5,3.5))
plu.box(stats, x='filter', y='n_records', ax=ax, x_order=order, color='steelblue')
plu.format_ax(ax, xlabel='', ylabel='n records', rotx=90, reduced_spines=True)
ax.set_yscale('log')
plt.tight_layout()
fig.savefig(os.path.join(path_figures, 'all_stats_box.pdf'))

##

fig, axs = plt.subplots(2,1,figsize=(3.5,4.5))

ax = axs[0]
plu.dist(df, x='AF_heart', ax=ax, color='steelblue')
plu.dist(df, x='AF_placenta', ax=ax, color='orange')
plu.format_ax(ax, xlabel='AF', ylabel='Density', reduced_spines=True)
plu.add_legend(ax=ax, colors={'Heart': 'steelblue', 'Placenta': 'orange'}, loc='upper right', bbox_to_anchor=(1,1))

ax = axs[1]
plu.dist(df, x='AD_heart', ax=ax, color='steelblue')
plu.dist(df, x='AD_placenta', ax=ax, color='orange')
plu.format_ax(ax, xlabel='AD', ylabel='Density', reduced_spines=True)

plt.tight_layout()
plt.show()
fig.savefig(os.path.join(path_figures, 'AF_ALL_FILTERED.pdf'))

df['AF_placenta'].max()


##



# PASS: final test
df['PASS'] = (df['AD_heart']>3) & \
             (df['DP_heart']>=10) & \
             (df['DP_placenta']>=10) & \
             (df['SB_pval']>=0.01) & \
             (df['AD_placenta']<2)
df_filtered = df[df['PASS']].copy()
df_filtered

##



fig, ax = plt.subplots(figsize=(1.5,3.5))
df_ = df['PASS'].value_counts(normalize=True).to_frame('fraction').reset_index()
plu.bar(df_, x='PASS', y='fraction', ax=ax, color='steelblue')
plu.format_ax(ax, xlabel='', ylabel='Fraction of records', reduced_spines=True)
fig.tight_layout()
plt.show()
fig.savefig(os.path.join(path_figures, 'frac_second_filter.pdf'))

##

fig, axs = plt.subplots(1,2,figsize=(6,3.5))

ax = axs[0]
plu.dist(df_filtered, x='AF_heart', ax=ax, color='steelblue')
plu.dist(df_filtered, x='AF_placenta', ax=ax, color='orange')
plu.format_ax(ax, xlabel='AF', ylabel='Density', reduced_spines=True)
plu.add_legend(ax=ax, colors={'Heart': 'steelblue', 'Placenta': 'orange'}, loc='upper right', bbox_to_anchor=(1,1))

ax = axs[1]
plu.dist(df_filtered, x='AD_heart', ax=ax, color='steelblue')
plu.dist(df_filtered, x='AD_placenta', ax=ax, color='orange')
plu.format_ax(ax, xlabel='AD', ylabel='Density', reduced_spines=True)

plt.tight_layout()
plt.show()
fig.savefig(os.path.join(path_figures, 'AF_FILTERED.pdf'))


##


# Single LCM
total_merged_set = set(df['mutation_id'].unique())
filtered_merged_set = set(df_filtered['mutation_id'].unique())
LCM_set = set(df_LCM['mutation_id'].unique())

# Here ...
sets = {'First_filter': total_merged_set, 'Second_filter': filtered_merged_set, 'LCM': LCM_set}
labels = list(sets.keys())
overlap_matrix = pd.DataFrame(
    [[len(sets[a] & sets[b]) for b in labels] for a in labels],
    index=labels, columns=labels
)

fig, ax = plt.subplots(figsize=(3,3))
plu.plot_heatmap(overlap_matrix, ax=ax, palette='Blues', annot=True, fmt='d', annot_size=7)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'overlap_matrix.pdf'))
plt.plot()

common = list(total_merged_set & LCM_set)

x = df.loc[df['mutation_id'].isin(common)].groupby('mutation_id')['AF_heart'].mean().loc[common]
y = df_LCM.loc[df_LCM['mutation_id'].isin(common)].groupby('mutation_id')['VAF'].mean().loc[common]

import seaborn as sns
from scipy.stats import pearsonr

r, pval = pearsonr(x, y)

fig, ax = plt.subplots(figsize=(3.5,3.5))
ax.plot(x, y, 'ko', alpha=0.7)
sns.regplot(x=x, y=y, ax=ax, scatter=False)
plu.format_ax(ax, xlabel='AF (merged)', ylabel='AF (LCM)', 
              title=f'Pearson r = {r:.2f}, p = {pval:.2e}',
              reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'AF_correlation.pdf'))


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
                         (grouped_merged['AD_placenta']<2)

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