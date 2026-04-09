"""
FILTER.1.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import plotting_utils as plu
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
matplotlib.use('macOSX')
plu.set_rcParams()


##


colors = {
    "C>A": "#03BDEF",
    "C>G": "#010101",
    "C>T": "#E42926",
    "T>A": "#CBCACA",
    "T>C": "#A2CF63",
    "T>G": "#ECC7C5",
}
MUT_ORDER = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]


##


def calculate_sbs96(df, context='SBS96'):

    bases = ['A', 'C', 'G', 'T']
    total = len(df)
    groups = ['SBS6', context]
    counts = (
        df.groupby(groups)
        .size()
        .div(total)
        .reset_index(name='fraction')
    )

    # Build complete index of all 96 contexts
    ctx_per_mut = {}
    for mut in MUT_ORDER:
        ref, alt = mut[0], mut[2]
        if context == 'SBS96':
            ctx_per_mut[mut] = sorted([f"{p}[{ref}>{alt}]{n}" for p in bases for n in bases])
        else:
            ctx_per_mut[mut] = sorted([f"{p}{ref}{n}" for p in bases for n in bases])

    full_idx = pd.MultiIndex.from_tuples(
        [(mut, ctx) for mut in MUT_ORDER for ctx in ctx_per_mut[mut]],
        names=['SBS6', context]
    )
    counts = (
        counts.set_index(groups)['fraction']
        .reindex(full_idx, fill_value=0)
        .reset_index()
    )

    return counts


##


def mut_profile(df=None, counts=None, context='SBS96', figsize=(12, 3), legend_kwargs={}) -> matplotlib.figure.Figure:
    """
    Plot raw fraction of MT-SNVs across SBS96 (or 3nt) contexts,
    stratified by mutation type (one axis each).
    """
    if counts is None:
        total = len(df)
        counts = calculate_sbs96(df, context=context)
    else:
        total = None

    fig, axs = plt.subplots(
        1, len(MUT_ORDER), figsize=figsize, sharey=True,
        constrained_layout=True
    )

    for i, mut in enumerate(MUT_ORDER):
        ax = axs[i]
        df_ = counts.query('SBS6 == @mut')
        x_order = sorted(df_[context].unique())
        plu.bar(
            df_, x=context, y='fraction',
            color=colors[mut],
            x_order=x_order,
            width=0.8, alpha=1.0, edgecolor=None,
            with_label=False, ax=ax
        )
        n_mut = int(round(df_['fraction'].sum() * total)) if total is not None else None
        plu.format_ax(
            ax, xlabel=mut, rotx=90,
            title=f'n: {n_mut}' if n_mut is not None else '',
            ylabel='Fraction of total SBSs' if i == 0 else '',
            reduced_spines=True, xticks_size=6
        )

    plu.add_legend(
        ax=axs[-1],
        colors={'H': '#444444', 'L': '#bbbbbb'},
        label='Strand', ncols=1,
        loc='upper left', bbox_to_anchor=(1, 1),
        **legend_kwargs
    )

    return fig


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
df = pd.read_csv(os.path.join(path_filtered, 'ALL_FILTERED_ctx.tsv.gz'), sep='\t')
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
df['AD_placenta'].describe()
df['AF_heart'].describe()
df['AF_placenta'].describe()
df['mutation_id'].nunique()
df.columns

##

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
ax.axvline(1000, color='darkorange', linestyle='--')

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

# Spectrum
df_unique = df.drop_duplicates('mutation_id').dropna(subset=['SBS6'])
counts = calculate_sbs96(df_unique, context='SBS96')
fig = mut_profile(counts=counts, context='SBS96', figsize=(12, 3))
fig.savefig(os.path.join(path_figures, 'spectrum_BASE_filter.pdf'))

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
plu.plot_heatmap(overlap_matrix, ax=ax, palette='Blues', annot=True, fmt='d', annot_size=10)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'overlap_matrix.pdf'))

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

##

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
fig.savefig(os.path.join(path_figures, 'common_muts_correlation.pdf'))

##

# What's missing?
df_missing = df_LCM.loc[~df_LCM['mutation_id'].isin(BASE_set)].copy()
df_present = df_LCM.loc[df_LCM['mutation_id'].isin(BASE_set)].copy()
df_ = pd.concat([
    df_missing.groupby('mutation_id')['Sample_ID'].nunique().to_frame('n_samples').reset_index().assign(status='Missing'),
    df_present.groupby('mutation_id')['Sample_ID'].nunique().to_frame('n_samples').reset_index().assign(status='Present')
])
df_.groupby('status')['n_samples'].describe()

np.sum(df_missing.groupby('mutation_id')['Sample_ID'].nunique()>2)

fig, axs = plt.subplots(1,3,figsize=(9,3))

ax = axs[0]
colors = {'Missing': 'grey', 'Present': 'darkorange'}
plu.dist(df_missing, x='NV', ax=ax, color='grey')
plu.dist(df_present, x='NV', ax=ax, color='darkorange')
plu.format_ax(ax, xlabel='AD', ylabel='Density', reduced_spines=True)
plu.add_legend(ax=ax, colors=colors, 
               loc='upper right', bbox_to_anchor=(1,1))

ax = axs[1]
plu.dist(df_missing, x='VAF', ax=ax, color='grey')
plu.dist(df_present, x='VAF', ax=ax, color='darkorange')
plu.format_ax(ax, xlabel='AF', ylabel='Density', reduced_spines=True)

ax = axs[2]
df_ = pd.concat([
    df_missing.groupby('mutation_id')['Sample_ID'].nunique().to_frame('n_samples').reset_index().assign(status='Missing'),
    df_present.groupby('mutation_id')['Sample_ID'].nunique().to_frame('n_samples').reset_index().assign(status='Present')
])
plu.dist(df_, x='n_samples', ax=ax, by='status', categorical_cmap=colors)
plu.format_ax(ax, xlabel='N samples', ylabel='Density', reduced_spines=True)

plt.tight_layout()
fig.savefig(os.path.join(path_figures, 'missing_from_BASE_filter.pdf'))


##


# FILTER.1

# Add stats and flags
df.dropna(subset=['SBS6'], inplace=True)

# Filter criteria:
df['PASS'] = (df['SB_pval']>0.1) & \
             (df['median_BQ']>=30) & \
             (df['MPOS']>5) & \
             (df['AD_placenta']<2) & \
             (df['NLOD']>10) & \
             (df['TLOD']>25) & \
             (df['POPAF']>3)

# Filter
df.loc[lambda x: x['PASS']]['mutation_id'].nunique()
df_filtered = df[df['PASS']].copy()
df_filtered['mutation_id'].nunique()

df_filtered['DP_placenta'].describe()

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


# Overlap
BASE_set = set(df_filtered['mutation_id'].unique())
LCM_set = set(df_LCM['mutation_id'].unique())
len(BASE_set)
len(LCM_set)

sets = {'FILTER.1': BASE_set, 'LCM': LCM_set}
labels = list(sets.keys())
overlap_matrix = pd.DataFrame(
    [[len(sets[a] & sets[b]) for b in labels] for a in labels],
    index=labels, columns=labels
)

fig, ax = plt.subplots(figsize=(2.5,2.5))
plu.plot_heatmap(overlap_matrix, ax=ax, palette='Blues', annot=True, fmt='d', annot_size=10)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'overlap_matrix.1.pdf'))

# Correlation commons
common = list(BASE_set & LCM_set)
x = (
    df_filtered.loc[df_filtered['mutation_id'].isin(common)].groupby('mutation_id')
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


##


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
fig.savefig(os.path.join(path_figures, 'common_muts_correlation.1.pdf'))


##


# Spectrum
df_unique = df_filtered.drop_duplicates('mutation_id').dropna(subset=['SBS6'])
counts = calculate_sbs96(df_unique, context='SBS96')
fig = mut_profile(counts=counts, context='SBS96', figsize=(12, 3))
fig.savefig(os.path.join(path_figures, 'spectrum_FILTER.1.pdf'))


##


# Analysis of SBS6 characteristics
fig, axs = plt.subplots(1,2,figsize=(7.5, 3))

ax = axs[0]
plu.violin(df_filtered, 'SBS6', 'AF_heart', ax=ax, categorical_cmap=colors, x_order=MUT_ORDER)
plu.strip(df_filtered, 'SBS6', 'AF_heart', ax=ax, categorical_cmap=colors, x_order=MUT_ORDER)
plu.format_ax(ax, xlabel='', ylabel='AF', reduced_spines=True)

ax = axs[1]
df_ = (
    df_filtered.groupby('mutation_id')
    ['region'].nunique()
    .to_frame('n regions')
    .join(
        df_filtered[['mutation_id', 'SBS6']]
        .drop_duplicates()
        .set_index('mutation_id')
    )
)
plu.violin(df_, 'SBS6', 'n regions', ax=ax, categorical_cmap=colors, x_order=MUT_ORDER)
plu.strip(df_, 'SBS6', 'n regions', ax=ax, categorical_cmap=colors, x_order=MUT_ORDER)
plu.format_ax(ax, xlabel='', ylabel='n regions', reduced_spines=True)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'SBS6_AF_nregions.pdf'))

##

fig, axs = plt.subplots(1,2,figsize=(7.5, 3))

ax = axs[0]
plu.violin(df_filtered, 'SBS6', 'AD_heart', ax=ax, categorical_cmap=colors, x_order=MUT_ORDER)
plu.strip(df_filtered, 'SBS6', 'AD_heart', ax=ax, categorical_cmap=colors, x_order=MUT_ORDER)
plu.format_ax(ax, xlabel='', ylabel='AD', reduced_spines=True)

ax = axs[1]
plu.violin(df_filtered, 'SBS6', 'SB_pval', ax=ax, categorical_cmap=colors, x_order=MUT_ORDER)
plu.strip(df_filtered, 'SBS6', 'SB_pval', ax=ax, categorical_cmap=colors, x_order=MUT_ORDER)
plu.format_ax(ax, xlabel='', ylabel='SB p-value', reduced_spines=True)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'SBS6_AD_SB_pval.pdf'))


# Save 
df_filtered.to_csv(os.path.join(path_filtered, 'FILTER.1.tsv'), sep='\t', index=False)


##


# FILTER.2
df = pd.read_csv(os.path.join(path_filtered, 'FILTER.1.tsv'), sep='\t')
df.columns

df['PASS'] = False
df['PASS'] = (df['SB_pval']>0.75)

df_filtered = df[df['PASS']].copy()
df_filtered['mutation_id'].nunique()

df_unique = df_filtered.drop_duplicates('mutation_id').dropna(subset=['SBS6'])
counts = calculate_sbs96(df_unique, context='SBS96')
fig = mut_profile(counts=counts, context='SBS96', figsize=(12, 3))
fig.savefig(os.path.join(path_figures, 'spectrum_FILTER.2.pdf'))



df_unique.groupby('SBS6').size()