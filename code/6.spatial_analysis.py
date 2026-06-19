"""
Spatial analysis: searching from asymmetric mutations in the heart.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import plotting_utils as plu
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import pairwise_distances
from scipy.cluster.hierarchy import linkage, leaves_list
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
matplotlib.use('macOSX')
plu.set_rcParams()


##


class Arrow3D(FancyArrowPatch):
    """
    A 3D arrow drawn as a FancyArrowPatch (clean 2D-style arrowhead).
    """
    def __init__(self, x0, y0, z0, x1, y1, z1, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz0 = (x0, y0, z0)
        self._xyz1 = (x1, y1, z1)

    def do_3d_projection(self, renderer=None):
        (x0, y0, z0), (x1, y1, z1) = self._xyz0, self._xyz1
        xs, ys, zs = proj_transform((x0, x1), (y0, y1), (z0, z1), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        return min(zs)


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
path_filtered = os.path.join(path_main, 'results')
path_figures = os.path.join(path_main, 'figures')

# Read data
df = pd.read_csv(os.path.join(path_filtered, 'ALLELIC_TABLE_FINAL.tsv.gz'), sep='\t')
xyz = pd.read_csv(os.path.join(path_data, 'Heart_final_coorindates_135.csv'))
xyz.rename(columns={'name':'Sample_ID'}, inplace=True)
samples = df['Sample_ID'].unique()
xyz.query('Sample_ID in @samples', inplace=True)


##


# Regional burdens
fig, axs = plt.subplots(1,2,figsize=(7,3.5))

df_ = (
    df
    .query('tissue=="heart" and in_sensible')
    .groupby(['Sample_ID', 'region'])
    ['mutation_id'].nunique().to_frame('n')
    .reset_index()
)

ax = axs[0]
x_order = df_.groupby('region')['n'].mean().sort_values(ascending=False).index
plu.box(df_, x='region', y='n', color='white', ax=ax, x_order=x_order)
plu.strip(df_, x='region', y='n', ax=ax, x_order=x_order)
plu.format_ax(ax=ax, xlabel='', ylabel='Number of SNVs', rotx=90, reduced_spines=True)

ax = axs[1]
# For each region, split its detected mutations (sensible callset) into those
# private to the region vs shared with >=1 other region (presence in >=1 sample)
pres = (
    df.query('tissue=="heart" and in_sensible')
    [['region', 'mutation_id']].drop_duplicates()
)
pres['status'] = np.where(
    pres.groupby('mutation_id')['region'].transform('size') > 1, 'Shared', 'Private'
)
pivot = (
    pres.groupby(['region', 'status']).size()
    .unstack(fill_value=0).reindex(columns=['Private', 'Shared'], fill_value=0)
)
order = pivot.sum(axis=1).sort_values(ascending=False).index
sns.barplot(x=order, y=pivot.loc[order].sum(axis=1), color='#bcbddc', label='Shared', ax=ax)
sns.barplot(x=order, y=pivot.loc[order, 'Private'], color='#756bb1', label='Private', ax=ax)
plu.format_ax(ax=ax, xlabel='', ylabel='Number of SNVs', rotx=90, reduced_spines=True)
ax.legend(frameon=False, loc='upper right')

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'regional_burdens.pdf'))


##


# Relationship among septum sections
septum_df = df.loc[lambda x: x['chunk'].str.contains('ept')]
centre = set(septum_df.query('in_sensible and region=="Centre_septum"')['mutation_id'].unique())
left = set(septum_df.query('in_sensible and region=="Left_septum"')['mutation_id'].unique())
right = set(septum_df.query('in_sensible and region=="Right_septum"')['mutation_id'].unique())

I = np.zeros((3,3))
J = np.zeros((3,3))
D = np.zeros((3,3))
for i,x in enumerate([left, centre, right]):
    for j,y in enumerate([left, centre, right]):
        if i<=j:
            I[i,j] = len(x.intersection(y))
            I[j,i] = len(x.intersection(y))
            J[i,j] = len(x.intersection(y))/len(x.union(y))
            J[j,i] = len(x.intersection(y))/len(x.union(y))
            D[i,j] = len(x-y)
            D[j,i] = len(y-x)

I = pd.DataFrame(I.astype(int), index=['Left', 'Centre', 'Right'], columns=['Left', 'Centre', 'Right'])
D = pd.DataFrame(D.astype(int), index=['Left', 'Centre', 'Right'], columns=['Left', 'Centre', 'Right'])
J = pd.DataFrame(J.astype(float), index=['Left', 'Centre', 'Right'], columns=['Left', 'Centre', 'Right'])


fig, ax = plt.subplots(1,3,figsize=(4, 1.6))
plu.plot_heatmap(I, ax=ax[0], cb=False, annot=True, fmt='d', title='Intersection')
plu.plot_heatmap(D, ax=ax[1], cb=False, annot=True, fmt='d', title='Difference')
plu.plot_heatmap(J, ax=ax[2], cb=False, annot=True, fmt='.2f', title='Jaccard Index')
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'in_shared_septum_relationships.pdf'))


##


# Single-lineages heatmaps
right_muts = right - left
left_muts = left - right
centre_muts = centre - (left | right)
muts = right_muts | left_muts | centre_muts

X = (
    df
    .loc[lambda x: x['chunk'].str.contains('ept')]
    .query('mutation_id in @muts')
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

fig, axs = plt.subplots(2,1,figsize=(6,2.75), sharex=True)

ax = axs[0]
ax.imshow(X.loc[region_order, mut_order], cmap='afmhot_r', vmin=0, vmax=.2, aspect='auto')
plu.format_ax(ax, xticks=mut_order, yticks=region_order, rotx=90, xticks_size=6)
plu.add_cbar(X.values.flatten(), ax=ax, 
             label='AF', palette='afmhot_r', vmin=0, vmax=.2)

ax = axs[1]
ax.imshow(X.loc[region_order, mut_order], cmap='afmhot_r', vmin=0, vmax=.03, aspect='auto')
plu.format_ax(ax, xticks=mut_order, yticks=region_order, rotx=90, xticks_size=6)
plu.add_cbar(X.values.flatten(), ax=ax, 
             label='AF', palette='afmhot_r', vmin=0, vmax=.03)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'in_shared_septum_relationships.pdf'))


##

# Single-samples heatmap
X = (
    df
    .loc[lambda x: x['chunk'].str.contains('ept')]
    .query('in_sensible and mutation_id in @muts')
    .query('tissue!="placenta"')
    .pivot_table(index='Sample_ID', columns='mutation_id', values='AF').fillna(0)
)
D = pairwise_distances((X.values), metric='cosine')
D = rescale_distances(D)
order = leaves_list(linkage(D, method='average'))
region_order = X.index[order].tolist()
D = pd.DataFrame(D, index=X.index, columns=X.index)
D_muts = pairwise_distances((X.values.T), metric='cosine')
order = leaves_list(linkage(D_muts, method='average'))
mut_order = X.columns[order].tolist()

# Region annotation for each sample (rows)
sample_region = (
    df.loc[lambda x: x['chunk'].str.contains('ept')]
    .query('tissue!="placenta"')
    .drop_duplicates('Sample_ID')
    .set_index('Sample_ID')['region']
)
region_colors = {
    'Left_septum':   '#4C78A8',
    'Centre_septum': '#F58518',
    'Right_septum':  '#54A24B',
}
# Per-sample region colors, in clustered row order
row_colors = [region_colors[sample_region[s]] for s in region_order]

# Fig
fig, ax = plt.subplots(figsize=(6, 3.5))

# Heatmap
ax.imshow(X.loc[region_order, mut_order], cmap='afmhot_r', vmin=0, vmax=.1, aspect='auto')
plu.format_ax(ax, xticks=mut_order, yticks=[], rotx=90, xticks_size=6)

# Row annotation strip (left)
axins_row = ax.inset_axes((-0.055, 0, 0.05, 1))
cb_row = plt.colorbar(
    matplotlib.cm.ScalarMappable(
    cmap=matplotlib.colors.ListedColormap(row_colors[::-1])),
    cax=axins_row, orientation='vertical'
)
cb_row.ax.set(xticks=[], yticks=[])
cb_row.outline.set_linewidth(0.1)

# Region legend
plu.add_legend(colors=region_colors, label='Region', ax=ax,
               ticks_size=8, artists_size=7, label_size=8,
               loc='upper left', bbox_to_anchor=(0, 1.25), ncols=3)

# AF colorbar
plu.add_cbar(X.values.flatten(), ax=ax,
             label='AF', palette='afmhot_r', vmin=0, vmax=.1)

fig.subplots_adjust(left=0.1, right=0.7, top=0.85, bottom=0.3)
fig.savefig(os.path.join(path_figures, 'single_samples_septum_relationships.pdf'))


##


# Relationship among LS/RS, RV and LV
SEPTUM = "Right_septum"
s = set(df.query('in_sensible and region==@SEPTUM')['mutation_id'].unique())
rv = set(df.query('in_sensible and region=="Right_Ventricle"')['mutation_id'].unique())
lv = set(df.query('in_sensible and region=="Left_Ventricle"')['mutation_id'].unique())

I = np.zeros((3,3))
J = np.zeros((3,3))
D = np.zeros((3,3))
for i,x in enumerate([s, rv, lv]):
    for j,y in enumerate([s, rv, lv]):
        if i<=j:
            I[i,j] = len(x.intersection(y))
            I[j,i] = len(x.intersection(y))
            J[i,j] = len(x.intersection(y))/len(x.union(y))
            J[j,i] = len(x.intersection(y))/len(x.union(y))
            D[i,j] = len(x-y)
            D[j,i] = len(y-x)

s = 'LS' if SEPTUM == 'Left_septum' else 'RS'
I = pd.DataFrame(I.astype(int), index=[s, 'RV', 'LV'], columns=[s, 'RV', 'LV'])
D = pd.DataFrame(D.astype(int), index=[s, 'RV', 'LV'], columns=[s, 'RV', 'LV'])
J = pd.DataFrame(J.astype(float), index=[s, 'RV', 'LV'], columns=[s, 'RV', 'LV'])

fig, ax = plt.subplots(1,3,figsize=(4, 1.6))
plu.plot_heatmap(I, ax=ax[0], cb=False, annot=True, fmt='d', title='Intersection')
plu.plot_heatmap(D, ax=ax[1], cb=False, annot=True, fmt='d', title='Difference')
plu.plot_heatmap(J, ax=ax[2], cb=False, annot=True, fmt='.2f', title='Jaccard Index')
fig.tight_layout()
fig.savefig(os.path.join(path_figures, f'{SEPTUM}_RV_LV_relationships.pdf'))


##


# Ancestries
cs = set(septum_df.query('in_sensible and region=="Centre_septum"')['mutation_id'].unique())
ls = set(septum_df.query('in_sensible and region=="Left_septum"')['mutation_id'].unique())
rs = set(septum_df.query('in_sensible and region=="Right_septum"')['mutation_id'].unique())
ls = set(df.query('in_sensible and region=="Left_septum"')['mutation_id'].unique())
rs = set(df.query('in_sensible and region=="Right_septum"')['mutation_id'].unique())
rv = set(df.query('in_sensible and region=="Right_Ventricle"')['mutation_id'].unique())
lv = set(df.query('in_sensible and region=="Left_Ventricle"')['mutation_id'].unique())
rv_specific = rv - (ls | rs | lv)
lv_specific = lv - (ls | rs | rv)
muts_1 = muts | rv_specific | lv_specific
x = (cs|ls|rs)
len(x & (rv-lv))
print(x & (lv-rv))
len(x & (lv & rv))


##


# Region level
X = (
    df
    .query('mutation_id in @muts_1')
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

fig, axs = plt.subplots(2,1,figsize=(6,2.75), sharex=True)

ax = axs[0]
ax.imshow(X.loc[region_order, mut_order], cmap='afmhot_r', vmin=0, vmax=.2, aspect='auto')
plu.format_ax(ax, xticks=mut_order, yticks=region_order, rotx=90, xticks_size=6)
plu.add_cbar(X.values.flatten(), ax=ax, 
             label='AF', palette='afmhot_r', vmin=0, vmax=.2)

ax = axs[1]
ax.imshow(X.loc[region_order, mut_order], cmap='afmhot_r', vmin=0, vmax=.03, aspect='auto')
plu.format_ax(ax, xticks=mut_order, yticks=region_order, rotx=90, xticks_size=6)
plu.add_cbar(X.values.flatten(), ax=ax, 
             label='AF', palette='afmhot_r', vmin=0, vmax=.03)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'ventricles_septum_lineages.pdf'))


##


# Single-samples heatmap
X = (
    df
    .query('in_sensible and mutation_id in @muts_1')
    .query('tissue!="placenta"')
    .pivot_table(index='Sample_ID', columns='mutation_id', values='AF').fillna(0)
)
D = pairwise_distances((X.values), metric='cosine')
D = rescale_distances(D)
order = leaves_list(linkage(D, method='average'))
region_order = X.index[order].tolist()
D = pd.DataFrame(D, index=X.index, columns=X.index)
D_muts = pairwise_distances((X.values.T), metric='cosine')
order = leaves_list(linkage(D_muts, method='average'))
mut_order = X.columns[order].tolist()

# Region annotation for each sample (rows)
sample_region = (
    df
    .query('tissue!="placenta"')
    .drop_duplicates('Sample_ID')
    .set_index('Sample_ID')['region']
)
region_colors = {
    'Left_septum':   '#4C78A8',
    'Centre_septum': '#F58518',
    'Right_septum':  '#54A24B',
    'Left_Ventricle':  '#E45756',
    'Right_Ventricle': '#B279A2',
}
# Per-sample region colors, in clustered row order
row_colors = [region_colors[sample_region[s]] for s in region_order]

# Fig
fig, ax = plt.subplots(figsize=(6, 3.7))

# Heatmap
ax.imshow(X.loc[region_order, mut_order], cmap='afmhot_r', vmin=0, vmax=.1, aspect='auto')
plu.format_ax(ax, xticks=mut_order, yticks=[], rotx=90, xticks_size=6)

# Row annotation strip (left)
axins_row = ax.inset_axes((-0.055, 0, 0.05, 1))
cb_row = plt.colorbar(
    matplotlib.cm.ScalarMappable(
    cmap=matplotlib.colors.ListedColormap(row_colors[::-1])),
    cax=axins_row, orientation='vertical'
)
cb_row.ax.set(xticks=[], yticks=[])
cb_row.outline.set_linewidth(0.1)

# Region legend
plu.add_legend(colors=region_colors, label='Region', ax=ax,
               ticks_size=8, artists_size=7, label_size=8,
               loc='upper left', bbox_to_anchor=(0, 1.3), ncols=3)

# AF colorbar
plu.add_cbar(X.values.flatten(), ax=ax,
             label='AF', palette='afmhot_r', vmin=0, vmax=.1)

fig.subplots_adjust(left=0.1, right=0.7, top=0.85, bottom=0.3)
fig.savefig(os.path.join(path_figures, 'single_samples_septum_relationships.pdf'))


##


# Single lineages in space
fig = plt.figure(figsize=(3.5, 3.5))

# Map each spatial sample to its region (grey if absent from final table)
xyz_region = xyz['Sample_ID'].map(
    df.drop_duplicates('Sample_ID').set_index('Sample_ID')['region']
)
space_colors = {
    'Left_Ventricle':  '#E45756',
    'Right_Ventricle': '#B279A2',
    'Left_septum':     '#4C78A8',
    'Centre_septum':   '#F58518',
    'Right_septum':    '#54A24B',
}
point_colors = xyz_region.map(space_colors)

ax = fig.add_subplot(projection='3d')
ax.computed_zorder = False  # respect explicit zorder so axes sit behind dots
ax.scatter(
    xyz['x'], xyz['y'], xyz['z'],
    s=20, c=point_colors, edgecolor='white', linewidth=.3, alpha=.9, depthshade=False,
    zorder=5
)
# Lighter grid, transparent panes, hide axis lines, tick marks and labels
for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
    axis._axinfo['grid'].update(color=(0, 0, 0, .12), linewidth=.2)
    axis.set_pane_color((1, 1, 1, 0))
    axis.line.set_linewidth(0)
    axis._axinfo['tick']['inward_factor'] = 0
    axis._axinfo['tick']['outward_factor'] = 0
ax.set(xticklabels=[], yticklabels=[], zticklabels=[])

# View so the origin corner sits at the back-bottom of the box
ax.view_init(elev=20, azim=50)

# Custom x/y/z axes emanating from the back-bottom corner, out along the box edges
xmin, xmax = xyz['x'].min(), xyz['x'].max()
ymin, ymax = xyz['y'].min(), xyz['y'].max()
zmin, zmax = xyz['z'].min(), xyz['z'].max()
x0, y0, z0 = xmin, ymin, zmin
pad = .15
axes_ends = {
    'x': (xmax + (xmax - xmin) * pad, ymin, zmin),
    'y': (xmin, ymax + (ymax - ymin) * pad, zmin),
    'z': (xmin, ymin, zmax + (zmax - zmin) * pad),
}
gap = .07
for label, (xe, ye, ze) in axes_ends.items():
    ax.add_artist(Arrow3D(
        x0, y0, z0, xe, ye, ze,
        mutation_scale=9, lw=.5, arrowstyle='-|>', color='k',
        shrinkA=0, shrinkB=0, zorder=0
    ))
    lx, ly, lz = x0 + (xe - x0) * (1 + gap), y0 + (ye - y0) * (1 + gap), z0 + (ze - z0) * (1 + gap)
    ax.text(lx, ly, lz, label, fontsize=8, ha='center', va='center')

# Region legend
plu.add_legend(
    colors={**space_colors}, label='Region', ax=ax,
    artists_size=6, ticks_size=6, label_size=7,
    loc='upper left', bbox_to_anchor=(.95, .95)
)

fig.subplots_adjust(left=.02, right=.78, top=.98, bottom=.02)
fig.savefig(os.path.join(path_figures, 'samples_3D_space.pdf'))



##


# Single mutations AF in space


"""
RS --> LV only ['chr19_4221919_G_A', 'chr18_34950831_G_A']
RS --> both 'chr20_48472304_C_T'

LS --> LV only ['chr7_50307645_G_A', 'chr4_46740534_C_T', 'chr15_45356903_C_T', 'chr2_78727500_C_T', 'chr6_14589088_G_A']
LS --> both ['chr16_48345763_G_A']

CS --> both ['chr3_140441942_C_T']

CS + LS --> RV only chr4_182179415_T_C
CS + LS --> LS only ['chr4_46740534_C_T', 'chr7_50307645_G_A', 'chr18_76379465_G_A', 'chr15_45356903_C_T', 'chr2_78727500_C_T', 'chr20_47905512_C_T', 'chr6_14589088_G_A']
"""


# Select muts
# mut_order[-5]
mut = 'chr9_8283815_C_T'
df.query('mutation_id == @mut')['Sample_ID'].nunique()

##

# Figure here
af = xyz['Sample_ID'].map(
    df
    .query('mutation_id == @mut')
    .drop_duplicates('Sample_ID')
    .set_index('Sample_ID')['AF']
).fillna(0)

fig = plt.figure(figsize=(5.2, 2.8))
gs = fig.add_gridspec(1, 2, width_ratios=[2, 1], wspace=.5)

# Left: AF in space (3D)
ax = fig.add_subplot(gs[0, 0], projection='3d')
ax.computed_zorder = False  # respect explicit zorder so axes sit behind dots
vmax = .1
ax.scatter(
    xyz['x'], xyz['y'], xyz['z'],
    s=20, c=af, cmap='afmhot_r', vmin=0, vmax=vmax,
    edgecolor='#bbbbbb', linewidth=.3, depthshade=False, zorder=5
)
# Lighter grid, transparent panes, hide axis lines, tick marks and labels
for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
    axis._axinfo['grid'].update(color=(0, 0, 0, .12), linewidth=.2)
    axis.set_pane_color((1, 1, 1, 0))
    axis.line.set_linewidth(0)
    axis._axinfo['tick']['inward_factor'] = 0
    axis._axinfo['tick']['outward_factor'] = 0
ax.set(xticklabels=[], yticklabels=[], zticklabels=[])

# View so the origin corner sits at the back-bottom of the box
ax.view_init(elev=20, azim=50)

# Custom x/y/z axes emanating from the back-bottom corner, out along the box edges
xmin, xmax = xyz['x'].min(), xyz['x'].max()
ymin, ymax = xyz['y'].min(), xyz['y'].max()
zmin, zmax = xyz['z'].min(), xyz['z'].max()
x0, y0, z0 = xmin, ymin, zmin
pad = .15
axes_ends = {
    'x': (xmax + (xmax - xmin) * pad, ymin, zmin),
    'y': (xmin, ymax + (ymax - ymin) * pad, zmin),
    'z': (xmin, ymin, zmax + (zmax - zmin) * pad),
}
gap = .07
for label, (xe, ye, ze) in axes_ends.items():
    ax.add_artist(Arrow3D(
        x0, y0, z0, xe, ye, ze,
        mutation_scale=9, lw=.5, arrowstyle='-|>', color='k',
        shrinkA=0, shrinkB=0, zorder=0
    ))
    lx, ly, lz = x0 + (xe - x0) * (1 + gap), y0 + (ye - y0) * (1 + gap), z0 + (ze - z0) * (1 + gap)
    ax.text(lx, ly, lz, label, fontsize=8, ha='center', va='center')

# AF colorbar
plu.add_cbar(af.values, ax=ax, label='AF', palette='afmhot_r', vmin=0, vmax=vmax)

# Add subplot here
df_ = (
    df.query('mutation_id == @mut')
    .groupby('region')['AD_alt'].sum()
    .to_frame('AD_alt')
    .reset_index()
)
region_abbr = {
    'Left_septum':     'LS',
    'Right_septum':    'RS',
    'Centre_septum':   'CS',
    'Left_Ventricle':  'LV',
    'Right_Ventricle': 'RV',
    'Placenta':        'P',
}
df_['region'] = df_['region'].map(region_abbr)

gs_bar = gs[0, 1].subgridspec(3, 1, height_ratios=[1, 3, 1])
ax_bar = fig.add_subplot(gs_bar[1, 0])
plu.bar(df_, 'region', 'AD_alt', ax=ax_bar)
plu.format_ax(ax=ax_bar, xlabel='', ylabel='AD', rotx=90, reduced_spines=True)

fig.suptitle(mut, x=.5)
fig.subplots_adjust(left=.02, right=.92, top=.9, bottom=.2)
fig.savefig(os.path.join(path_figures, f'{mut}_AF_3D_space.pdf'))


##