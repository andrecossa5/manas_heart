"""
Spatial analysis.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import plotting_utils as plu
import matplotlib.pyplot as plt
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
COMP = str.maketrans('ACGT', 'TGCA')


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


def DA_muts(df, groupby: str, groups: list[str]|str):
    """
    Differential abundance of mutations between samples in `groups` vs the rest, 
    where groups are defined by `groupby` (e.g. region).
    """

    # Pivot
    AD = df.pivot_table(index='Sample_ID', columns='mutation_id', values='AD_alt').fillna(0)
    DP = df.pivot_table(index='Sample_ID', columns='mutation_id', values='DP').fillna(0)
    AF = df.pivot_table(index='Sample_ID', columns='mutation_id', values='AF').fillna(0)

    # Get groups
    groups = [groups] if isinstance(groups, str) else groups
    group_samples = df.loc[df[groupby].isin(groups), 'Sample_ID'].unique()
    rest_samples = df.loc[~df[groupby].isin(groups), 'Sample_ID'].unique()

    # Stats
    AD_group = AD.loc[group_samples]
    AD_rest = AD.loc[rest_samples]
    DP_group = DP.loc[group_samples]
    DP_rest = DP.loc[rest_samples]
    AF_group = AF.loc[group_samples]
    AF_rest = AF.loc[rest_samples]
    AF_group_mean = AF_group.mean()
    AF_rest_mean = AF_rest.mean()
    AF_group_max = AF_group.max()
    AF_rest_max = AF_rest.max()
    prevalence_group = (AF_group>0).sum() / group_samples.size
    prevalence_rest = (AF_rest>0).sum() / rest_samples.size
    pseudobulk_AF_group = AD_group.sum() / DP_group.sum()
    pseudobulk_AF_rest = AD_rest.sum() / DP_rest.sum()
    FC = (pseudobulk_AF_group - pseudobulk_AF_rest) / (pseudobulk_AF_rest + 10**(-10))
    pvals = [ mannwhitneyu(AF_group[mut], AF_rest[mut])[1] for mut in AF.columns ]

    # Package results
    results = pd.DataFrame({
        'AF_group_mean': AF_group_mean,
        'AF_rest_mean': AF_rest_mean,
        'AF_group_max': AF_group_max,
        'AF_rest_max': AF_rest_max,
        'FC': FC,
        'prevalence_group': prevalence_group,
        'prevalence_rest': prevalence_rest,
        'pseudobulk_AF_group': pseudobulk_AF_group,
        'pseudobulk_AF_rest': pseudobulk_AF_rest,
        'pval': pvals
    })
    results = results.sort_values('pseudobulk_AF_group', ascending=False)
    
    return results


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


def draw_heatmap(
    df: pd.DataFrame,
    muts_heart: list[str],
    muts_regions: list[str],
    region_order: list[str],
    cmap: str = 'mako',
    vmax: float | None = None,
    plot: str = 'AF',
    ax: plt.Axes | None = None,
    ):
    """
    Draw heatmap of mutation values across samples.
    """

    # AF and AD matrices
    muts = list(set(muts_heart) | set(muts_regions))
    AF = df.pivot_table(index='Sample_ID', columns='mutation_id', values='AF').fillna(0)[muts]
    AD = df.pivot_table(index='Sample_ID', columns='mutation_id', values='AD_alt').fillna(0)[muts]

    # Sample -> region map
    sample_region = (
        df[['Sample_ID', 'region']].drop_duplicates()
        .set_index('Sample_ID')['region']
    )
    sample_region = sample_region.reindex(AF.index)

    # --- Row order: samples grouped by region (clustering order) ---
    row_order = (
        pd.DataFrame({'region': sample_region}, index=sample_region.index)
        .assign(region_rank=lambda d: pd.Categorical(d['region'], categories=region_order, ordered=True).codes)
        .sort_values('region_rank')
        .index.tolist()
    )

    # --- Column order: heart block first, then per-region blocks ---
    heart_samples = [s for s in df.loc[df['tissue']=='heart', 'Sample_ID'].unique() if s in AF.index]
    heart_prev = (AF.loc[heart_samples, muts_heart] > 0).mean()
    heart_cols = heart_prev.sort_values(ascending=False).index.tolist()
    region_only = [m for m in muts_regions if m not in set(heart_cols)]
    prev_per_region = {}
    for r in region_order:
        samples_r = [s for s in df.loc[df['region']==r, 'Sample_ID'].unique() if s in AF.index]
        prev_per_region[r] = (AF.loc[samples_r, region_only] > 0).mean()
    prev_df = pd.DataFrame(prev_per_region)
    dom_region = prev_df.idxmax(axis=1)

    region_cols = []
    for r in region_order:
        block = prev_df.loc[dom_region == r, r].sort_values(ascending=False).index.tolist()
        region_cols.extend(block)

    col_order = heart_cols + region_cols

    # Prep plotting matrix and vmin/vmax
    if plot == 'AD':
        X_heat = AD.loc[row_order, col_order]
        vmax = 3 if vmax is None else vmax
    else:
        X_heat = AF.loc[row_order, col_order]
        vmax = np.percentile(X_heat.values[X_heat.values > 0], 95) if vmax is None else vmax

    # Ax
    ax.imshow(X_heat.values, aspect='auto', cmap=cmap, vmin=0, vmax=vmax)

    # Lines to separate heart block and region blocks
    if heart_cols and region_cols:
        ax.axvline(len(heart_cols) - 0.5, color='white', lw=0.5)
    offset = len(heart_cols)
    for r in region_order[:-1]:
        offset += int((dom_region == r).sum())
        ax.axvline(offset - 0.5, color='white', lw=0.5)
    row_region_seq = sample_region.loc[row_order].values
    boundaries = np.where(row_region_seq[:-1] != row_region_seq[1:])[0]
    for b in boundaries:
        ax.axhline(b + 0.5, color='white', lw=0.5)

    # Column block labels (above heatmap)
    block_starts, block_labels = [], []
    if heart_cols:
        block_starts.append(0)
        block_labels.append('Heart')
    offset = len(heart_cols)
    for r in region_order:
        block_size = int((dom_region == r).sum())
        if block_size > 0:
            block_starts.append(offset)
            block_labels.append(r)
            offset += block_size
    block_ends = block_starts[1:] + [offset]
    last_idx = len(block_labels) - 1
    for i, (start, end, label) in enumerate(zip(block_starts, block_ends, block_labels)):
        if i == last_idx:
            x, ha = end - 0.5, 'right'
        else:
            x, ha = start - 0.5, 'left'
        ax.text(
            x, -0.5, label,
            ha=ha, va='bottom', fontsize=8,
            rotation=0, rotation_mode='anchor', clip_on=False,
        )

    # Cosmetic
    plu.format_ax(
        ax, yticks=[], xticks=[],
        ylabel=f'Sample (n={X_heat.shape[0]})',
        xlabel=f'SNV (n={X_heat.shape[1]})'
    )
    plu.add_cbar(
        X_heat.values.flatten(), 
        ax=ax, palette=cmap, vmin=0, vmax=vmax, 
        label='n reads' if plot=='AD' else 'AF'
    )

    return ax


##


# Paths
path_main = '/Users/cossa/Desktop/projects/manas_heart'
path_data = os.path.join(path_main, 'data')
path_input = os.path.join(path_main, 'data/input')
path_filtered = os.path.join(path_main, 'results')
path_figures = os.path.join(path_main, 'figures')


##


# Read forcecall results
df = pd.read_csv(os.path.join(path_filtered, 'ALLELIC_TABLE_NO_ARTIFACTS.tsv.gz'), sep='\t')
df['mutation_id'].nunique()

##

# Whole heart Differential Abundance (DA)
results = DA_muts(df, 'tissue', groups='heart')
muts_heart = (
    results
    .query('pval<=0.01 and prevalence_group>=0.75 and prevalence_rest<=0.1')
    .index.to_list()
)
len(muts_heart)

# Single-regions
L = []
for group in df['region'].unique():
    results = DA_muts(df, 'region', groups=group)
    L.append(results)
results = pd.concat(L)
muts_regions = (
    results 
    .query('pval<=0.01 and prevalence_group>=0.3 and FC>=.1 and prevalence_rest<=0.1')
    .index.to_list()
)
len(muts_regions)

# Combine
muts = list(set(muts_heart) | set(muts_regions)) 
len(muts)

# Spectrum differences
fig = mut_profile(df=df.query('mutation_id in @muts'))
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'forcecall_enriched_spectrum.pdf'))


##


# Cluster regions by their pseudobulk muts profiles
AF = df.pivot_table(index='Sample_ID', columns='mutation_id', values='AF').fillna(0)
AD = df.pivot_table(index='Sample_ID', columns='mutation_id', values='AD_alt').fillna(0)

##

# Cluster regions by their pseudobulk muts profiles
X = (
    df# .query('mutation_id in @muts')
    .groupby(['mutation_id', 'region'])
    [['AD_alt', 'DP']].sum()
    .reset_index()
    .assign(AF=lambda x: x['AD_alt'] / (x['DP'] + 10**(-18)))
    .pivot(index='region', columns='mutation_id', values='AF').fillna(0)
)

D = pairwise_distances(X, metric='cosine')
order = leaves_list(linkage(D, method='average'))
region_order = X.index[order].tolist()
D = rescale_distances(D)
D = pd.DataFrame(D, index=X.index, columns=X.index)

fig, ax = plt.subplots(figsize=(3.5,3.5))
ax.imshow(D.iloc[order, order], cmap='Spectral', vmin=0, vmax=1)
plu.format_ax(ax, xticks=D.columns[order], yticks=D.index[order], rotx=90)
plu.add_cbar(D.values.flatten(), ax=ax, 
             label='Cosine distance', palette='Spectral', vmin=0, vmax=1)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'region_clustering.pdf'))

##

fig, ax = plt.subplots(figsize=(10,4))
draw_heatmap(df, muts_heart, muts_regions, region_order=region_order, plot='AF', ax=ax, vmax=0.1)
fig.tight_layout()
plu.save_best_pdf_quality(
    fig, (10,4), path_figures, 'heatmap_AF.pdf', 1000

)



