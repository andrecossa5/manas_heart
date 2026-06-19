"""
Final checks on final callset.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import plotting_utils as plu
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from sklearn.metrics import pairwise_distances
from scipy.cluster.hierarchy import linkage, leaves_list
matplotlib.use('macOSX')
plu.set_rcParams()


##


sb6_colors = {
    "C>A": "#03BDEF",
    "C>G": "#010101",
    "C>T": "#E42926",
    "T>A": "#CBCACA",
    "T>C": "#A2CF63",
    "T>G": "#ECC7C5",
}
MUT_ORDER = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]


##


def rescale_distances(D):
    """
    Rescale (row-wise) pairwise distances to [0,1].
    """
    min_dist = D[~np.eye(D.shape[0], dtype=bool)].min()
    max_dist = D[~np.eye(D.shape[0], dtype=bool)].max()
    D = (D-min_dist)/(max_dist-min_dist)
    np.fill_diagonal(D, 0)
    return D


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

    df = df.drop_duplicates('mutation_id') if df is not None else None

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
            color=sb6_colors[mut],
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
path_filtered = os.path.join(path_main, 'results')
path_figures = os.path.join(path_main, 'figures')

# Read data
df = pd.read_csv(os.path.join(path_filtered, 'ALLELIC_TABLE_FILTERED_ANNOTATED.tsv.gz'), sep='\t')
tree_assignment = pd.read_csv(
    os.path.join(
        path_data, 
        'Filteredmutations_14061_Sample_subset_snv_assigned_to_branches.txt'
    ), sep='\t'
)
df_LCM = pd.read_csv(os.path.join(path_data, 'Heart_metadata.csv'))

##


# Preliminary: check artifacts spectrum
fig = mut_profile(df.query('artifact_flag == "Artefact"'), context='SBS96', figsize=(12, 3))
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'forcecall_artifact_spectrum.pdf'))

fig = mut_profile(df.query('artifact_flag == "No artefact"'), context='SBS96', figsize=(12, 3))
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'forcecall_no_artifact_spectrum.pdf'))

fig = mut_profile(df, context='SBS96', figsize=(12, 3))
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'forcecall_all_spectrum.pdf'))

##

# Filter artifacts flag, extract final_callset
df_filtered = df.query('AF>0 and artifact_flag == "No artefact"')
final_callset = set(df_filtered['mutation_id'].unique())


##


# 1. How many final muts were present in the original single-LCM calls?

# LCM
common = set(df['Sample_ID'].unique()) & set(df_LCM['Sample_ID'].unique())      
df_LCM = df_LCM.query('NV>0 and Sample_ID in @common')
lcm = set(df_LCM['mutation_id'].unique())

# Merged bams
df_merged = pd.read_csv(os.path.join(path_filtered, 'FILTER.1.tsv'), sep='\t')
merged = set(df_merged['mutation_id'].unique())

# Calculate overlaps (7 disjoint Venn regions: lcm, merged, final_callset)
only_lcm         = len(lcm - merged - final_callset)
only_merged      = len(merged - lcm - final_callset)
only_final       = len(final_callset - lcm - merged)
lcm_merged       = len((lcm & merged) - final_callset)
lcm_final        = len((lcm & final_callset) - merged)
merged_final     = len((merged & final_callset) - lcm)
lcm_merged_final = len(lcm & merged & final_callset)

# Venn diagram
fig, ax = plt.subplots(figsize=(3.5, 3.5))
v = venn3(
    subsets=(
        only_lcm, only_merged, lcm_merged,
        only_final, lcm_final, merged_final, lcm_merged_final
    ),
    set_labels=('LCM', 'Merged', 'Final callset'),
    set_colors=('#03BDEF', '#A2CF63', '#E42926'),
    ax=ax
)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'callset_venn.pdf'))


##


# 2. Where were final mutations assigned to the phylogeny?
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

fig, ax = plt.subplots(figsize=(3.5, 3))
plu.counts_plot(df_filtered.drop_duplicates('mutation_id'), x='desc_samples_orgin', ax=ax)
plu.format_ax(ax, xlabel='', ylabel='n SNVs', rotx=90, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'tree_assignment.pdf'))


##


# 3. Spectra lcm only, merged only, final callset.
fig = mut_profile(df_filtered.query('desc_samples_orgin == "Unassigned"'), context='SBS96', figsize=(10, 2.5))
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'merged_final_spectrum.pdf'))

fig = mut_profile(df_filtered.query('desc_samples_orgin != "Unassigned"'), context='SBS96', figsize=(10, 2.5))
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'lcm_final_spectrum.pdf'))

fig = mut_profile(df_filtered, context='SBS96', figsize=(10, 2.5))
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'final_417_spectrum.pdf'))


##


# Annotate CpGs and calling strategy
CpGs = ['A[C>T]G', 'C[C>T]G', 'G[C>T]G', 'T[C>T]G']
df_filtered['calling_strategy'] = (
    np.where(df_filtered['desc_samples_orgin'] == 'Unassigned', 'Merged', 'LCM')
)
df_filtered['in_CpGs'] = (
    np.where(df_filtered['SBS96'].isin(CpGs), 'CpG', 'Non-CpG')
)


##


## 5. Signal sparsity with increasing resolution
fig, axs = plt.subplots(1,3,figsize=(9, 3))

ax = axs[0]
df_ = (
    df_filtered.groupby(['mutation_id', 'region'])
    ['AD_alt'].sum()
    .reset_index()
    .merge(df_filtered[['mutation_id', 'region', 'calling_strategy']].drop_duplicates(), 
           on=['mutation_id', 'region'], how='left')
)
median = df_['AD_alt'].median()
median_lcm = df_.query('calling_strategy == "LCM"')['AD_alt'].median()
median_merged = df_.query('calling_strategy == "Merged"')['AD_alt'].median()
plu.dist(df_, x='AD_alt', ax=ax)
plu.format_ax(ax, xlabel='AD_alt', ylabel='Density', rotx=0, reduced_spines=True, title='AD summed by region')
ax.axvline(median, color='red', linestyle='--')
ax.text(0.5, 0.95, f'Median: {median}', fontsize=8, transform=ax.transAxes)
ax.axvline(median_lcm, color='blue', linestyle='--')
ax.text(0.5, 0.90, f'Median LCM: {median_lcm}', fontsize=8, transform=ax.transAxes)
ax.axvline(median_merged, color='green', linestyle='--')
ax.text(0.5, 0.85, f'Median Merged: {median_merged}', fontsize=8, transform=ax.transAxes)

ax = axs[1]
df_ = (
    df_filtered.groupby(['mutation_id', 'chunk'])
    ['AD_alt'].sum()
    .reset_index()
    .merge(df_filtered[['mutation_id', 'chunk', 'calling_strategy']].drop_duplicates(), 
           on=['mutation_id', 'chunk'], how='left')
)
median = df_['AD_alt'].median()
median_lcm = df_.query('calling_strategy == "LCM"')['AD_alt'].median()
median_merged = df_.query('calling_strategy == "Merged"')['AD_alt'].median()
plu.dist(df_, x='AD_alt', ax=ax)
plu.format_ax(ax, xlabel='AD_alt', ylabel='Density', rotx=0, reduced_spines=True, title='AD summed by chunk')
ax.axvline(median, color='red', linestyle='--')
ax.text(0.5, 0.95, f'Median: {median}', fontsize=8, transform=ax.transAxes)
ax.axvline(median_lcm, color='blue', linestyle='--')
ax.text(0.5, 0.90, f'Median LCM: {median_lcm}', fontsize=8, transform=ax.transAxes)
ax.axvline(median_merged, color='green', linestyle='--')
ax.text(0.5, 0.85, f'Median Merged: {median_merged}', fontsize=8, transform=ax.transAxes)

ax = axs[2]
df_ = df_filtered.copy()
median = df_['AD_alt'].median()
median_lcm = df_.query('calling_strategy == "LCM"')['AD_alt'].median()
median_merged = df_.query('calling_strategy == "Merged"')['AD_alt'].median()

plu.dist(df_, x='AD_alt', ax=ax)
plu.format_ax(ax, xlabel='AD_alt', ylabel='Density', rotx=0, reduced_spines=True, title='single LCM sample calls')
ax.axvline(median, color='red', linestyle='--')
ax.text(0.5, 0.95, f'Median: {median}', fontsize=8, transform=ax.transAxes)
ax.axvline(median_lcm, color='blue', linestyle='--')
ax.text(0.5, 0.90, f'Median LCM: {median_lcm}', fontsize=8, transform=ax.transAxes)
ax.axvline(median_merged, color='green', linestyle='--')
ax.text(0.5, 0.85, f'Median Merged: {median_merged}', fontsize=8, transform=ax.transAxes)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'signal_sparsity.pdf'))


##


# 6. Sample level differences in LCM/merged, CpG/non-CpG calls
fig, ax = plt.subplots(1,4, figsize=(7, 2.5))

plu.box(df_filtered, x='in_CpGs', y='AF', ax=ax[0])
plu.format_ax(ax[0], xlabel='', reduced_spines=True)
plu.box(df_filtered, x='calling_strategy', y='AF', ax=ax[1])
plu.format_ax(ax[1], xlabel='', reduced_spines=True)

df_ = (
    df_filtered
    .groupby(['in_CpGs', 'mutation_id'])['Sample_ID'].nunique()
    .reset_index(name='n_samples')
)
plu.box(df_, x='in_CpGs', y='n_samples', ax=ax[2])
plu.format_ax(ax[2], xlabel='', reduced_spines=True)

df_ = (
    df_filtered
    .groupby(['calling_strategy', 'mutation_id'])['Sample_ID'].nunique()
    .reset_index(name='n_samples')
)
plu.box(df_, x='calling_strategy', y='n_samples', ax=ax[3])
plu.format_ax(ax[3], xlabel='', reduced_spines=True)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'AF_sharedness.pdf'))


##


# 7. Final callsets: sensible vs shared

# Sensible
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

# Shared
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

# df_filtered.query('in_sensible')['mutation_id'].nunique()
# df_filtered.query('in_shared')['mutation_id'].nunique()


##


# 8. Differences in SBS96 spectra, sensible vs shared
fig = mut_profile(df_filtered.query('in_shared'), context='SBS96', figsize=(10, 2.5))
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'shared_spectrum.pdf'))

fig = mut_profile(df_filtered.query('in_sensible'), context='SBS96', figsize=(10, 2.5))
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'sensible_spectrum.pdf'))


##


# 9. Visualize callset distances and mutations AF
fig, axs = plt.subplots(1,3,figsize=(9,3.5))

ticks_sizes = {'region': 8, 'chunk': 6, 'Sample_ID': 3}
callset = 'complete'
for i, group in enumerate(['region', 'chunk', 'Sample_ID']):

    ax = axs[i]
    X = (
        df_filtered
        #.query('in_sensible')
        .query('tissue!="placenta"')
        .groupby(['mutation_id', group])
        [['AD_alt', 'DP']].sum()
        .reset_index()
        .assign(AF=lambda x: x['AD_alt'] / (x['DP'] + 10**(-18)))
        .pivot(index=group, columns='mutation_id', values='AF').fillna(0)
    )

    # Distances 
    D = pairwise_distances((X.values), metric='cosine')
    D = rescale_distances(D)
    order = leaves_list(linkage(D, method='average'))
    region_order = X.index[order].tolist()
    D = pd.DataFrame(D, index=X.index, columns=X.index)

    ax.imshow(D.iloc[order, order], cmap='Spectral', vmin=0.1, vmax=.9)
    plu.format_ax(
        ax, 
        xticks=D.columns[order], 
        yticks=[], 
        xlabel='', ylabel='',
        rotx=90, 
        xticks_size=ticks_sizes[group],
        title=group
    )
    plu.add_cbar(D.values.flatten(), ax=ax, 
                 label='Distance', palette='Spectral', vmin=0.1, vmax=.9)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, f'{callset}_distance_clustering.pdf'))

# Muts
X = (
    df_filtered
    # .query('in_sensible')
    .query('tissue!="placenta"')
    .groupby(['mutation_id', 'region'])
    [['AD_alt', 'DP']].sum()
    .reset_index()
    .assign(AF=lambda x: x['AD_alt'] / (x['DP'] + 10**(-18)))
    .pivot(index='region', columns='mutation_id', values='AF').fillna(0)
)
D = pairwise_distances((X.values), metric='cosine')
D = rescale_distances(D)
order = leaves_list(linkage(D, method='average'))
region_order = X.index[order].tolist()
D = pd.DataFrame(D, index=X.index, columns=X.index)
D_muts = pairwise_distances((X.values.T), metric='cosine')
order = leaves_list(linkage(D_muts, method='average'))
mut_order = X.columns[order].tolist()

fig, axs = plt.subplots(2,1,figsize=(10,4), sharex=True)

ax = axs[0]
ax.imshow(X.loc[region_order, mut_order], cmap='afmhot_r', vmin=0, vmax=.2, aspect='auto')
plu.format_ax(ax, xticks=mut_order, yticks=region_order, rotx=90, xticks_size=3)
plu.add_cbar(X.values.flatten(), ax=ax, 
             label='AF', palette='afmhot_r', vmin=0, vmax=.2)

ax = axs[1]
ax.imshow(X.loc[region_order, mut_order], cmap='afmhot_r', vmin=0, vmax=.03, aspect='auto')
plu.format_ax(ax, xticks=mut_order, yticks=region_order, rotx=90, xticks_size=3)
plu.add_cbar(X.values.flatten(), ax=ax, 
             label='AF', palette='afmhot_r', vmin=0, vmax=.03)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, f'{callset}_{group}_distance_clustering.pdf'))



##


# Write
(
    df_filtered
    .to_csv(os.path.join(path_filtered, 'ALLELIC_TABLE_FINAL.tsv.gz'), sep='\t', index=False)
)


##