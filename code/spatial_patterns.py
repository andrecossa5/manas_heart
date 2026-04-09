
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
df = pd.read_csv(os.path.join(path_filtered, 'FILTER.1.tsv'), sep='\t')


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
ax.imshow(D.iloc[order, order], cmap='Spectral', vmin=0.1, vmax=.9)
plu.format_ax(ax, xticks=D.columns[order], yticks=D.index[order], rotx=90)
plu.add_cbar(D.values.flatten(), ax=ax, 
             label='Cosine distance', palette='Spectral', vmin=0.1, vmax=0.9)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'cosine_distance_regions.pdf'))




# (df.pivot(index='chunk', columns='mutation_id', values='AF_heart').fillna(0)>0).sum(axis=1)
# (df.pivot(index='chunk', columns='mutation_id', values='AF_heart').fillna(0)>0).sum(axis=0).sort_values(ascending=False).head(20)

df.query('mutation_id=="chr18_37006496_T_A"')[['region', 'mutation_id', 'AF_heart', 'AD_heart', 'DP_heart']]