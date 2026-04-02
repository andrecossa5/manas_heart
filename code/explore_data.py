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
# stats.query('filter=="passed"')['n_records'].sum()

# Original muts
df_LCM = pd.read_csv(os.path.join(path_data, 'Heart_metadata.csv'))

# Filter on the same samples
samples_merged = df_samples['Sample_ID'].unique()
samples_LCM = df_LCM['Sample_ID'].unique()
common = list(set(samples_merged) & set(samples_LCM))
df_LCM = df_LCM.query('NV>=3 and Sample_ID in @common')

# Checks
assert df['chunk'].isin(df_samples['chunk']).all()

##


df['DP_placenta'].describe()
df['DP_heart'].describe()
df['AD_heart'].describe()
df['DP_placenta'].describe()
df['AF_heart'].describe()
df['AF_placenta'].describe()


# Explore 
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

fig, axs = plt.subplots(1,3,figsize=(9,3))

ax = axs[0]
plu.dist(df, x='AD_heart', ax=ax, color='steelblue')
plu.dist(df, x='AD_placenta', ax=ax, color='darkorange')
plu.format_ax(ax, xlabel='AD', ylabel='Density', reduced_spines=True)
plu.add_legend(ax=ax, colors={'Heart': 'steelblue', 'Placenta': 'darkorange'}, loc='upper right', bbox_to_anchor=(1,1))
ax.axvline(100, color='darkorange', linestyle='--')

ax = axs[1]
plu.dist(df, x='DP_heart', ax=ax, color='steelblue')
plu.dist(df, x='DP_placenta', ax=ax, color='darkorange')
plu.format_ax(ax, xlabel='DP', ylabel='Density', reduced_spines=True)

ax = axs[2]
plu.dist(df, x='AF_heart', ax=ax, color='steelblue')
plu.dist(df, x='AF_placenta', ax=ax, color='darkorange')
plu.format_ax(ax, xlabel='AF', ylabel='Density', reduced_spines=True)

plt.tight_layout()
fig.savefig(os.path.join(path_figures, 'dist_BASE_filter.pdf'))

##

# Overlap
BASE_set = set(df['mutation_id'].unique())
LCM_set = set(df_LCM['mutation_id'].unique())
len(BASE_set)
len(LCM_set)

sets = {'BASE': BASE_set, 'LCM': LCM_set}
labels = list(sets.keys())
overlap_matrix = pd.DataFrame(
    [[len(sets[a] & sets[b]) for b in labels] for a in labels],
    index=labels, columns=labels
)

fig, ax = plt.subplots(figsize=(2.5,2.5))
plu.plot_heatmap(overlap_matrix, ax=ax, palette='Blues', annot=True, fmt='d', annot_size=7)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'overlap_matrix.pdf'))
plt.plot()

# Correlation commons
common = list(BASE_set & LCM_set)
x = (
    df.loc[df['mutation_id'].isin(common)].groupby('mutation_id')
    [['AD_heart', 'DP_heart']].sum()
    .loc[common]
    .assign(AF_heart=lambda x: x['AD_heart'] / (x['DP_heart'] + 10**-18))
)  
y = (
    df_LCM.loc[df_LCM['mutation_id'].isin(common)].groupby('mutation_id')
    [['NV', 'NR']].sum()
    .loc[common]
    .assign(AF_heart=lambda x: x['NV'] / (x['NR'] + 10**-18))
) 

import seaborn as sns
from scipy.stats import pearsonr

fig, axs = plt.subplots(1,2,figsize=(6.5,3))

ax = axs[0]
r, pval = pearsonr(x['AD_heart'], y['NV'])
ax.plot(x['AD_heart'], y['NV'], 'ko', alpha=0.7)
sns.regplot(x=x['AD_heart'], y=y['NV'], ax=ax, scatter=False)
plu.format_ax(ax, xlabel='Summed AD chunks', ylabel='Summed AD LCMs', 
              title=f'Pearson r = {r:.2f}, p = {pval:.2e}',
              reduced_spines=True)

ax = axs[1]
r, pval = pearsonr(x['AF_heart'], y['AF_heart'])
ax.plot(x['AF_heart'], y['AF_heart'], 'ko', alpha=0.7)
sns.regplot(x=x['AF_heart'], y=y['AF_heart'], ax=ax, scatter=False)
plu.format_ax(ax, xlabel='Pseudo-bulk AF (chunks)', ylabel='Pseudo-bulk AF (LCMs)', 
              title=f'Pearson r = {r:.2f}, p = {pval:.2e}',
              reduced_spines=True)

fig.tight_layout()
plt.show()

fig.savefig(os.path.join(path_figures, 'common_muts_correlation.pdf'))


##


# FILTER.1

# Add stats and flags
df['INDEL'] = df['ALT'].map(lambda x: True if len(x)>1 else False)

# t = 100
# np.sum(df.loc[df['AD_placenta']==0, 'NLOD']>=t) / np.sum(df['AD_placenta']==0)

df['PASS'] = (df['SB_pval']>0.1) & \
             (df['AD_placenta']<2) & \
             (df['AF_heart']<0.5) & \
             (df['AF_heart']>.01) & \
             (df['NLOD']>=100) & \
             (~df['INDEL'])


df_filtered = df[df['PASS']].copy()
df_filtered['mutation_id'].nunique()


##


fig, ax = plt.subplots(figsize=(1.5,3.5))
df_ = df['PASS'].value_counts(normalize=True).to_frame('fraction').reset_index()
plu.bar(df_, x='PASS', y='fraction', ax=ax, color='steelblue')
plu.format_ax(ax, xlabel='', ylabel='Fraction of records', reduced_spines=True)
fig.tight_layout()
plt.show()
fig.savefig(os.path.join(path_figures, 'frac_FILTER.1.pdf'))


##

fig, axs = plt.subplots(1,3,figsize=(9,3))

ax = axs[0]
plu.dist(df_filtered, x='AD_heart', ax=ax, color='steelblue')
plu.dist(df_filtered, x='AD_placenta', ax=ax, color='darkorange')
plu.format_ax(ax, xlabel='AD', ylabel='Density', reduced_spines=True)
plu.add_legend(ax=ax, colors={'Heart': 'steelblue', 'Placenta': 'darkorange'}, loc='upper right', bbox_to_anchor=(1,1))
ax.axvline(100, color='darkorange', linestyle='--')

ax = axs[1]
plu.dist(df_filtered, x='DP_heart', ax=ax, color='steelblue')
plu.dist(df_filtered, x='DP_placenta', ax=ax, color='darkorange')
plu.format_ax(ax, xlabel='DP', ylabel='Density', reduced_spines=True)

ax = axs[2]
plu.dist(df_filtered, x='AF_heart', ax=ax, color='steelblue')
plu.dist(df_filtered, x='AF_placenta', ax=ax, color='darkorange')
plu.format_ax(ax, xlabel='AF', ylabel='Density', reduced_spines=True)

plt.tight_layout()
fig.savefig(os.path.join(path_figures, 'dist_FILTER.1_filter.pdf'))


##