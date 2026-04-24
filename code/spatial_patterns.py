
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


# Paths
path_main = '/Users/cossa/Desktop/projects/manas_heart'
path_data = os.path.join(path_main, 'data')
path_input = os.path.join(path_main, 'data/input')
path_filtered = os.path.join(path_main, 'results')
path_figures = os.path.join(path_main, 'figures')


##


# Read fitlered muts
df = pd.read_csv(os.path.join(path_filtered, 'FILTER.2.tsv'), sep='\t')
df['mutation_id'].nunique()


##


# Burden
X = (
    df.groupby(['mutation_id', 'region'])
    [['AD_heart', 'DP_heart']].sum()
    .reset_index()
    .assign(AF_heart=lambda x: x['AD_heart'] / (x['DP_heart'] + 10**(-18)))
    .pivot(index='region', columns='mutation_id', values='AF_heart').fillna(0)
)

fig, ax = plt.subplots(figsize=(3.5,3.5))

df_ = (
    (X>0).sum(axis=1)
    .to_frame('n')
    .sort_values('n', ascending=False)
    .reset_index()
)
order = df_['region'].to_list()
plu.bar(df_, 'region', 'n', ax=ax, x_order=order)
plu.format_ax(ax, xlabel='', ylabel='n SBSs', reduced_spines=True, rotx=90)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'regions_burden.pdf'))



##


# Clustering
D = pairwise_distances(X, metric='cosine')
order = leaves_list(linkage(D, method='average'))
D = rescale_distances(D)
D = pd.DataFrame(D, index=X.index, columns=X.index)

fig, ax = plt.subplots(figsize=(3.5,3.5))
ax.imshow(D.iloc[order, order], cmap='Spectral', vmin=0, vmax=1)
plu.format_ax(ax, xticks=D.columns[order], yticks=D.index[order], rotx=90)
plu.add_cbar(D.values.flatten(), ax=ax, 
             label='Cosine distance', palette='Spectral', vmin=0, vmax=1)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'cosine_distance_regions.pdf'))


##


# Clonal structure
order = D.index[order]
REGION_ORDER = ['Centre_septum', 'Right_Ventricle', 'Left_Ventricle', 'Left_septum', 'Right_septum']

##


# AF per chunks
X = (
    df
    .pivot(index='chunk', columns='mutation_id', values='AF_heart')
    .fillna(0)
)

# --- Row order: chunks grouped by region, within region by total AF descending ---
chunk_region = (
    df[['chunk', 'region']].drop_duplicates()
    .set_index('chunk')['region']
)
chunk_region = chunk_region[chunk_region.isin(REGION_ORDER)]
chunk_total_af = X.sum(axis=1)
row_order = (
    pd.DataFrame({'region': chunk_region, 'total_af': chunk_total_af}, index=chunk_region.index)
    .assign(region_rank=lambda d: pd.Categorical(d['region'], categories=REGION_ORDER, ordered=True).codes)
    .sort_values(['region_rank', 'total_af'], ascending=[True, False])
    .index
)
row_order = [c for c in row_order if c in X.index]
X = X.loc[row_order]

# --- Column order: block diagonal — dominant region per mutation, then AF descending ---
X_region = X.copy()
X_region.index = chunk_region.loc[X_region.index]
X_region = X_region.groupby(level=0).mean().reindex(REGION_ORDER).fillna(0)

dom_region = X_region.idxmax(axis=0)
dom_af     = X_region.max(axis=0)
col_order = (
    pd.DataFrame({'dom_region': dom_region, 'dom_af': dom_af})
    .assign(region_rank=lambda d: pd.Categorical(d['dom_region'], categories=REGION_ORDER, ordered=True).codes)
    .sort_values(['region_rank', 'dom_af'], ascending=[True, False])
    .index
)
X = X[col_order]

##

vmax = np.percentile(X.values[X.values > 0], 99) if (X.values > 0).any() else 1

fig, ax = plt.subplots(figsize=(8,3))
ax.imshow(X.values, aspect='auto', cmap='afmhot_r', vmin=0, vmax=vmax)
plu.format_ax(ax, yticks=X.index, xticks=[], ylabel='Chunk', xlabel=f'SNV (n={X.shape[1]})', rotx=90)
plu.add_cbar(X.values.flatten(), ax=ax, palette='afmhot_r', vmin=0, vmax=vmax, label='Allelic frequency')
fig.tight_layout()
plu.save_best_pdf_quality(fig, (8,3), path_figures, 'af_chunks.pdf')
plt.show()
