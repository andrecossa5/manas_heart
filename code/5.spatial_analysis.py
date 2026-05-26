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
from scipy.cluster.hierarchy import linkage, leaves_list, cophenet
from scipy.cluster.hierarchy import linkage as _link, cophenet
from sklearn.decomposition import NMF
from scipy.spatial.distance import squareform
from scipy.stats import mannwhitneyu, fisher_exact, chi2_contingency
from statsmodels.stats.multitest import multipletests
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


def fisher_asymmetry(df, regions, region_col='region', muts=None):
    """
    Per-mutation Fisher's exact test on pseudobulk read counts across `regions`.

    H0: P(read=alt | region) is equal across all `regions` => clone contributes equally.
    H1: contribution is asymmetric.

    For K=2 regions: 2x2 test; effect size AI = (AF_0 - AF_1)/(AF_0 + AF_1).
    For K=3 regions: 2x3 test; effect sizes AI_LR (first vs last) and f_C (centre fraction).

    Caveat: pseudobulk pools chunks within a region, ignoring chunk-level overdispersion.
    """

    if muts is not None:
        df = df[df['mutation_id'].isin(muts)]
    df = df[df[region_col].isin(regions)]

    AD = (
        df.pivot_table(index='mutation_id', columns=region_col, values='AD_alt', aggfunc='sum')
        .reindex(columns=regions).fillna(0)
    )
    DP = (
        df.pivot_table(index='mutation_id', columns=region_col, values='DP', aggfunc='sum')
        .reindex(columns=regions).fillna(0)
    )
    REF = DP - AD
    AF = AD / DP.replace(0, np.nan)

    pvals, ors = [], []
    K = len(regions)
    for mut in AD.index:
        table = np.vstack([AD.loc[mut].values, REF.loc[mut].values])
        if table.sum() == 0 or (table.sum(axis=1) == 0).any() or (table.sum(axis=0) == 0).any():
            pvals.append(1.0); ors.append(np.nan); continue
        if K == 2:
            res = fisher_exact(table)
            pvals.append(res.pvalue); ors.append(res.statistic)
        else:
            # scipy<1.15 fisher_exact supports only 2x2; use chi-square for 2xK.
            # Pseudobulk read counts are large => asymptotics are fine.
            chi2_res = chi2_contingency(table, correction=False)
            pvals.append(chi2_res.pvalue); ors.append(np.nan)

    out = pd.DataFrame(index=AD.index)
    for r in regions:
        out[f'AF_{r}'] = AF[r]
        out[f'AD_{r}'] = AD[r]
        out[f'DP_{r}'] = DP[r]

    eps = 1e-12
    if len(regions) == 2:
        a, b = regions
        out['AI'] = (AF[a].fillna(0) - AF[b].fillna(0)) / (AF[a].fillna(0) + AF[b].fillna(0) + eps)
        out['odds_ratio'] = ors
    else:
        a, b = regions[0], regions[-1]
        out['AI_LR'] = (AF[a].fillna(0) - AF[b].fillna(0)) / (AF[a].fillna(0) + AF[b].fillna(0) + eps)
        if len(regions) == 3:
            c = regions[1]
            total = AF[a].fillna(0) + AF[b].fillna(0) + AF[c].fillna(0)
            out['f_C'] = AF[c].fillna(0) / (total + eps)

    out['pval'] = pvals
    out['qval'] = multipletests(pvals, method='fdr_bh')[1]
    out = out.sort_values('qval')

    return out


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


# Regional burdens
fig, ax = plt.subplots(figsize=(2.5,3.5))
plu.counts_plot(df.drop_duplicates(['mutation_id', 'region']).query('tissue=="heart"'), 'region', ax=ax)
plu.format_ax(ax=ax, xlabel='', ylabel='Number of SNVs', rotx=90, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'regional_burdens.pdf'))

##

fig, ax = plt.subplots(figsize=(2.5,3.5))
df_ = (
    df
    .query('tissue=="heart"')
    .groupby(['Sample_ID', 'region'])
    ['mutation_id'].nunique().to_frame('n')
    .reset_index()
)
x_order = df_.groupby('region')['n'].mean().sort_values(ascending=False).index

plu.box(df_, x='region', y='n', color='white', ax=ax, x_order=x_order)
plu.strip(df_, x='region', y='n', ax=ax, x_order=x_order)
plu.format_ax(ax=ax, xlabel='', ylabel='Number of SNVs', rotx=90, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'regional_burdens_box.pdf'))


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
    df.query('mutation_id in @muts')
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


## 


# Contribution Ventricles to septum sections

# Stage 1: ventricular origin per mutation (LV vs RV, 2x2 Fisher)
vent = fisher_asymmetry(df, ['Left_Ventricle', 'Right_Ventricle'], muts=None)

# Stage 2: septal asymmetry across the three septum sections (2x3 Fisher)
sept = fisher_asymmetry(
    df, ['Left_septum', 'Centre_septum', 'Right_septum'], muts=None
)
sept = (
    sept[['AF_Left_septum', 'AF_Centre_septum', 'AF_Right_septum', 'pval']]
    .query('pval <= 0.05')
)
sept

# Get ventriculars
(
    vent.loc[vent.index.isin(sept.index)]
    [['AF_Left_Ventricle', 'AF_Right_Ventricle', 'pval']]
    .sort_values('pval')
)



##


# Clonal decomposition via Poisson-NMF (KL-loss) on AD counts.
# X (regions x mutations) of alt-read counts factorized as W (regions x K) @ H (K x mutations).
# Biology:  W[r,k] = abundance of clone k in region r;  H[k,m] = mutation m's loading on clone k.
# KL-loss is the Poisson model on counts: high-coverage entries naturally weigh more.
# Per-region coverage differences end up in W magnitudes and are normalized away by W_frac.
# L1 sparsity on H concentrates each clone's fingerprint on a few mutations.

# Region-level pseudobulk alt-read counts
pb = df.groupby(['Sample_ID', 'mutation_id'])[['AD_alt', 'DP']].sum().reset_index()
AD_mat = pb.pivot(index='Sample_ID', columns='mutation_id', values='AD_alt').fillna(0)
X_counts = AD_mat.values.astype(float)

# Pick K by cophenetic correlation of consensus W-clustering (Brunet et al. 2004).
# Idea: for each K, run NMF n_runs times with different seeds. For each run, assign
# each sample to its dominant clone (argmax over W row) and build a connectivity matrix
# C (1 if two samples share a clone, else 0). Average C across runs -> consensus matrix.
# Stable K -> consensus is near-binary -> cophenetic correlation of (1 - consensus) is high.
# Pick the largest K before the cophenetic drop.

Ks = list(range(2, 30))
n_runs = 100
n_samples = X_counts.shape[0]
coph = []
for K in Ks:
    C = np.zeros((n_samples, n_samples))
    for seed in range(n_runs):
        nmf = NMF(
            n_components=K, init='random', max_iter=2000, random_state=seed,
            beta_loss='kullback-leibler', solver='mu',
            alpha_H=0.1, alpha_W=0.0, l1_ratio=1.0,
        )
        W_run = nmf.fit_transform(X_counts)
        assign = np.argmax(W_run, axis=1)
        C += (assign[:, None] == assign[None, :]).astype(float)
    C /= n_runs
    # Cophenetic correlation of (1 - consensus) as distance
    dist = 1.0 - C
    np.fill_diagonal(dist, 0.0)
    Z = _link(squareform(dist, checks=False), method='average')
    coph_corr, _ = cophenet(Z, squareform(dist, checks=False))
    coph.append(coph_corr)
    print(f'K={K}  cophenetic={coph_corr:.3f}')

# Diagnose
print(np.argmax(coph), Ks[np.argmax(coph)], coph[np.argmax(coph)])

# Choose more stable K
fig, ax = plt.subplots(figsize=(3.5, 2.5))
ax.plot(Ks, coph, marker='o')
plu.format_ax(
    ax=ax, xlabel='K (n clones)', ylabel='Cophenetic correlation'
)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'nmf_K_selection.pdf'))


##


# Refit at chosen K with multiple seeds; keep the best
K = Ks[np.argmax(coph)]
best = None
for seed in range(100):
    nmf = NMF(
        n_components=K, init='random', max_iter=4000, random_state=seed,
        beta_loss='kullback-leibler', solver='mu',
        alpha_H=0.1, alpha_W=0.0, l1_ratio=1.0,
    )
    W = nmf.fit_transform(X_counts)
    if best is None or nmf.reconstruction_err_ < best[0]:
        best = (nmf.reconstruction_err_, W, nmf.components_)

_, W, H = best
W = pd.DataFrame(W, index=AD_mat.index, columns=[f'C{i+1}' for i in range(K)])
H = pd.DataFrame(H, index=W.columns, columns=AD_mat.columns)
W_frac = W.div(W.sum(axis=1), axis=0)


##


# Plot clonal mutational fingerprints
D = pairwise_distances(H.values, metric='cosine')
order = leaves_list(linkage(D, method='average'))
clone_order = H.index[order].tolist()

muts = []
for clone in clone_order:
    top_muts = H.loc[clone].sort_values(ascending=False).head(5).index.tolist()
    muts.extend(top_muts)
muts_ = []
for m in muts:
    if m not in muts_:
        muts_.append(m)

fig, ax = plt.subplots(figsize=(6,4.5))
ax.imshow(
    H.loc[clone_order, muts_], 
    aspect='auto', cmap='afmhot_r',
    vmax=np.percentile(H.loc[:, muts_].values, 98),
    vmin=np.percentile(H.loc[:, muts_].values, 2)
)
plu.format_ax(ax, xlabel='SNVs', ylabel='Clone', 
              xticks=muts_, yticks=clone_order, rotx=90)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'clonal_mutational_fingerprints.pdf'))


##


# Plot clonal composition across regions
region_order = [
    'Right_Ventricle', 'Right_septum', 'Centre_septum', 
    'Left_septum', 'Left_Ventricle'
]
W_mean = (
    W.reset_index()
    .merge(df[['Sample_ID', 'region']].drop_duplicates(), on='Sample_ID')
    .groupby('region')[W.columns].mean()
    .loc[region_order]
)

# Order clones by weighted center-of-mass along region_order (0..n-1).
r_idx = np.arange(len(region_order))
clone_com = (W_mean.values * r_idx[:, None]).sum(axis=0) / \
            W_mean.values.sum(axis=0).clip(min=1e-12)
clone_order = W_mean.columns[np.argsort(clone_com)].tolist()
W_show = W_mean.loc[region_order, clone_order]

fig, ax = plt.subplots(figsize=(5,3))
ax.imshow(
    W_show.values, aspect='auto', cmap='Blues',
    vmax=np.percentile(W_mean.values, 99),
    vmin=np.percentile(W_mean.values, 1)
)
plu.format_ax(ax, xticks=W_show.columns, yticks=W_show.index, xlabel='Clone', ylabel='Region')
plu.add_cbar(
    W_mean.values.flatten(), ax=ax, palette='Blues',
    vmin=np.percentile(W_mean.values, 1),
    vmax=np.percentile(W_mean.values, 99),
    label='Mean W (clone abundance)'
)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'clonal_composition_across_regions.pdf'))


##


# Per-clone left-vs-right imbalance test.
# Left  = Left_Ventricle + Left_septum samples
# Right = Right_Ventricle + Right_septum samples  (Centre_septum excluded)
# Test: Mann-Whitney U on per-sample W[:, k] between Left and Right groups.
# Effect size: AI = (mean_W_right - mean_W_left) / (mean_W_right + mean_W_left)
#   AI > 0 -> right-biased; AI < 0 -> left-biased; |AI| in [0, 1].

sample_to_region = (
    df[['Sample_ID', 'region']].drop_duplicates().set_index('Sample_ID')['region']
)
left_regions  = ['Left_Ventricle',  'Left_septum']
right_regions = ['Right_Ventricle', 'Right_septum']
reg = sample_to_region.reindex(W.index)
left_samples  = W.index[reg.isin(left_regions)]
right_samples = W.index[reg.isin(right_regions)]

rows = []
eps = 1e-12
for k in W.columns:
    wL = W.loc[left_samples, k].values
    wR = W.loc[right_samples, k].values
    mL, mR = wL.mean(), wR.mean()
    AI = (mR - mL) / (mR + mL + eps)
    if wL.size == 0 or wR.size == 0 or (np.all(wL == 0) and np.all(wR == 0)):
        pval = np.nan
    else:
        pval = mannwhitneyu(wR, wL, alternative='two-sided').pvalue
    rows.append(
        {'clone': k, 'mean_W_left': mL, 
         'mean_W_right': mR,
        'AI': AI, 'pval': pval}
    )

# Refactor and plot as volcano plot
clone_LR = pd.DataFrame(rows).set_index('clone')
clone_LR['qval'] = multipletests(clone_LR['pval'].fillna(1.0), method='fdr_bh')[1]
clone_LR = clone_LR.loc[clone_order].sort_values('AI')
clone_LR['-log10(pval)'] = -np.log10(clone_LR['pval'] + 1e-12)

fig, ax = plt.subplots(figsize=(3,3))
plu.volcano(
    clone_LR, x='AI', y='-log10(pval)', xlim=(-1.5, 1.5), ax=ax, fig=fig, labels=clone_LR.index
)
ax.axvline(0, color='red', linestyle='--', lw=0.5)
ax.axhline(-np.log10(0.05), color='red', linestyle='--', lw=0.5)
plu.format_ax(
    ax=ax, reduced_spines=True, 
    xlabel='Asymmetry Index (AI)', ylabel='-log10(p-value)'
)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'clonal_LR_bias.pdf'))


##



