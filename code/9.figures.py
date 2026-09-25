"""
Figures for the lineage analysis, following FIGURES.md.

Analysis 1, lineage characterisation
  embryonic_lineage_features.pdf   counts and cell fraction of the two classes
  genotype_call_composition.pdf    present / absent / undetermined per class
  lineage_enrichment.pdf           region VAF over region enrichment
  enrichment_type.pdf              extended patches vs scattered
  example_enriched.pdf             3D examples of the two patterns
  lineage_sweep.pdf                cell fraction of every enriched lineage

Analysis 2, phylogenetic relationships
  phylogenetic_relationships.pdf   region dendrogram beside the region VAF matrix
  lcm_heatmap.pdf                  the same matrix at single-cut resolution
  distances.pdf                    region distances behind the dendrogram

Analysis 3
  septum_ventricle_asymmetry.pdf   ventricular lean, its test, and two lineages
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
from scipy.spatial.distance import squareform
from sklearn.metrics import pairwise_distances
from scipy.stats import binom
import plotting_utils as plu

matplotlib.use('macOSX')
plu.set_rcParams()


##


REGIONS = ['LS', 'CS', 'RS', 'LV', 'RV']
REGION_ABBR = {
    'Left_septum': 'LS', 'Centre_septum': 'CS', 'Right_septum': 'RS',
    'Left_Ventricle': 'LV', 'Right_Ventricle': 'RV',
}
COLORS = {
    'LS': '#4C78A8', 'CS': '#F58518', 'RS': '#54A24B',
    'LV': '#E45756', 'RV': '#B279A2',
}
CLASS_COLORS = {'Pre-dating': '#9e9ac8', 'Heart-specific': '#31a354'}
STATE_COLORS = {'Present': '#31a354', 'Absent': '#bdbdbd', 'Undetermined': '#f0f0f0'}
K_NEIGHBOURS = 3        # neighbourhood for extended vs scattered
N_BOOT = 1000
CHAR_FRACTION = 0.8
N_PERM_JOINT = 20000
POINT_SIZE = 22         # every 3D point the same size


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


def style_3d(ax, xyz):
    """
    Transparent panes, light grid, and x/y/z arrows from the back-bottom corner.
    """
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis._axinfo['grid'].update(color=(0, 0, 0, .12), linewidth=.2)
        axis.set_pane_color((1, 1, 1, 0))
        axis.line.set_linewidth(0)
        axis._axinfo['tick']['inward_factor'] = 0
        axis._axinfo['tick']['outward_factor'] = 0
    ax.set(xticklabels=[], yticklabels=[], zticklabels=[])
    ax.view_init(elev=20, azim=50)
    lo, hi = xyz.min(), xyz.max()
    ends = {
        'x': (hi['x'] + (hi['x'] - lo['x']) * .15, lo['y'], lo['z']),
        'y': (lo['x'], hi['y'] + (hi['y'] - lo['y']) * .15, lo['z']),
        'z': (lo['x'], lo['y'], hi['z'] + (hi['z'] - lo['z']) * .15),
    }
    for label, (xe, ye, ze) in ends.items():
        ax.add_artist(Arrow3D(lo['x'], lo['y'], lo['z'], xe, ye, ze, mutation_scale=7,
                              lw=.4, arrowstyle='-|>', color='k', shrinkA=0, shrinkB=0))
        ax.text(lo['x'] + (xe - lo['x']) * 1.07, lo['y'] + (ye - lo['y']) * 1.07,
                lo['z'] + (ze - lo['z']) * 1.07, label, fontsize=6, ha='center', va='center')


def draw_mutation_3d(ax, j, edge_by_region=False, ring=None):
    """
    One mutation in space: present cuts filled by AF, powered absences open,
    undetermined light grey. Every point is the same size.
    """
    state = STATE.values[:, j]
    pres = state == 'present'
    absent = state == 'absent'
    other = ~pres & ~absent
    ax.scatter(coords['x'][other], coords['y'][other], coords['z'][other],
               s=POINT_SIZE, c='#dcdcdc', edgecolor='none', depthshade=False, zorder=3)
    ax.scatter(coords['x'][absent], coords['y'][absent], coords['z'][absent],
               s=POINT_SIZE, facecolor='white', edgecolor='#9b9b9b', linewidth=.4,
               depthshade=False, zorder=4)
    ax.scatter(coords['x'][pres], coords['y'][pres], coords['z'][pres],
               s=POINT_SIZE, c=2 * AF[pres, j], cmap='afmhot_r', vmin=0, vmax=.6,
               edgecolor=[COLORS[r] for r in reg[pres]] if edge_by_region else 'k',
               linewidth=.9 if edge_by_region else .3, depthshade=False, zorder=5)
    if ring is not None and ring.any():
        ax.scatter(coords['x'][ring], coords['y'][ring], coords['z'][ring],
                   s=POINT_SIZE, facecolor='none', edgecolor='k', linewidth=1.6,
                   depthshade=False, zorder=6)
    style_3d(ax, coords)


##


path_main = '/Users/cossa/Desktop/projects/manas_heart'
path_results = os.path.join(path_main, 'results')
path_figures = os.path.join(path_main, 'figures')

geno = pd.read_csv(os.path.join(path_results, 'GENOTYPES_TRUE.tsv.gz'), sep='\t')
summary = pd.read_csv(os.path.join(path_results, 'LINEAGE_SUMMARY.tsv'), sep='\t').set_index('mutation_id')
enrich = pd.read_csv(os.path.join(path_results, 'REGION_ENRICHMENT.tsv'), sep='\t')
lean_tab = pd.read_csv(os.path.join(path_results, 'VENTRICULAR_LEAN.tsv'), sep='\t')
loo = pd.read_csv(os.path.join(path_results, 'VENTRICULAR_LEAN_LOO.tsv'), sep='\t', index_col=0)
xyz = pd.read_csv(os.path.join(path_main, 'data', 'Heart_final_coorindates_135.csv')).set_index('name')

geno['reg'] = geno['region'].map(REGION_ABBR)
AF_df = geno.pivot(index='Sample_ID', columns='mutation_id', values='AF')
AD_df = geno.pivot(index='Sample_ID', columns='mutation_id', values='AD_alt')
DP_df = geno.pivot(index='Sample_ID', columns='mutation_id', values='DP')
STATE = geno.pivot(index='Sample_ID', columns='mutation_id', values='state_af10')
cuts, muts = AF_df.index, AF_df.columns
AF = AF_df.values
meta = geno.drop_duplicates('Sample_ID').set_index('Sample_ID').loc[cuts]
reg, chunk = meta['reg'].values, meta['chunk'].values
present = (STATE == 'present').values
coords = xyz.loc[cuts, ['x', 'y', 'z']]
row_idx = {r: np.where(reg == r)[0] for r in REGIONS}
rng = np.random.default_rng(0)

klass = np.where(summary['pre_existing'].reindex(muts).values, 'Pre-dating', 'Heart-specific')
class_order = ['Pre-dating', 'Heart-specific']

# Region x mutation tables, and the enrichment ordering used by several figures
fisher_p = enrich.pivot_table(index='region', columns='mutation_id', values='p_fisher').reindex(REGIONS)[muts]
region_af = enrich.pivot_table(index='region', columns='mutation_id', values='AF').reindex(REGIONS)[muts]
logp = -np.log10(fisher_p.clip(lower=1e-10))
max_logp = np.nan_to_num(logp.values).max(0)
is_enriched = (fisher_p < 0.05).sum().reindex(muts).fillna(0).values > 0
enrich_order = np.concatenate([
    np.where(~is_enriched)[0][np.argsort(max_logp[~is_enriched])],
    np.where(is_enriched)[0][np.argsort(max_logp[is_enriched])],
])
n_present = present.sum(0)
mean_cf = np.array([2 * AF[present[:, j], j].mean() if present[:, j].any() else np.nan
                    for j in range(len(muts))])

print(f'{len(muts)} mutations x {len(cuts)} cuts | enriched {int(is_enriched.sum())} | '
      f'Pre-dating {(klass == "Pre-dating").sum()}')


##


# 1a. Counts and cell fraction of the two classes
fig = plt.figure(figsize=(7.5, 3))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 2], wspace=.35)

ax = fig.add_subplot(gs[0, 0])
counts = pd.Series(klass).value_counts().reindex(class_order)
ax.bar(range(2), counts.values, color=[CLASS_COLORS[k] for k in class_order], width=.6)
for i, v in enumerate(counts.values):
    ax.text(i, v + 1, str(v), ha='center', fontsize=8)
plu.format_ax(ax=ax, xticks=class_order, ylabel='n SNVs', rotx=20, reduced_spines=True)

ax = fig.add_subplot(gs[0, 1])
for k in class_order:
    sel = klass == k
    ax.scatter(n_present[sel], mean_cf[sel], s=16, color=CLASS_COLORS[k],
               alpha=.85, edgecolor='none', label=k)
for k, style in zip(class_order, ['-', '--']):
    sel = klass == k
    mx, my = np.nanmean(n_present[sel]), np.nanmean(mean_cf[sel])
    ax.axvline(mx, color=CLASS_COLORS[k], linestyle=style, linewidth=.8, alpha=.9)
    ax.axhline(my, color=CLASS_COLORS[k], linestyle=style, linewidth=.8, alpha=.9)
    ax.text(mx, ax.get_ylim()[1], f' {mx:.0f}', color=CLASS_COLORS[k], fontsize=6,
            va='top', ha='left')
    ax.text(ax.get_xlim()[0], my, f'{my:.2f} ', color=CLASS_COLORS[k], fontsize=6,
            va='bottom' if k == class_order[0] else 'top', ha='left')
plu.format_ax(ax=ax, xlabel='n LCM samples', ylabel='Mean cell fraction',
              reduced_spines=True)
ax.legend(frameon=False, fontsize=7, loc='lower right')

fig.tight_layout()
fig.subplots_adjust(bottom=.22)
fig.savefig(os.path.join(path_figures, 'embryonic_lineage_features.pdf'))
plt.show()


##


# 1a, companion. Composition of the genotype calls in each class
cells = pd.DataFrame({
    'Class': np.repeat(klass, len(cuts)),
    'Call': pd.Series(STATE.values.T.ravel()).map(
        {'present': 'Present', 'absent': 'Absent',
         'weak': 'Undetermined', 'undetermined': 'Undetermined'}).values,
})
cells['Class'] = pd.Categorical(cells['Class'], categories=class_order)
cells['Call'] = pd.Categorical(cells['Call'], categories=list(STATE_COLORS))

fig, ax = plt.subplots(figsize=(6.2, 1.9))
plu.bb_plot(cells, cov1='Class', cov2='Call', categorical_cmap=STATE_COLORS, ax=ax)
plu.format_ax(ax=ax, xlabel='Fraction of genotype calls')
plu.add_legend(colors=STATE_COLORS, label='Genotype call', ax=ax, ticks_size=6,
               artists_size=5, label_size=7, loc='upper left',
               bbox_to_anchor=(1.01, 1), ncols=1)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'genotype_call_composition.pdf'))
plt.show()


##


# 1b. Region VAF over region enrichment, mutations ordered by enrichment
cols = enrich_order
split_at = int((~is_enriched).sum()) - .5
vmax_af = np.nanpercentile(region_af.values, 98)
vmax_p = np.nanpercentile(logp.values, 99)
cmap_af = matplotlib.colormaps['afmhot_r'].with_extremes(bad='#dddddd')
cmap_p = matplotlib.colormaps['viridis'].with_extremes(bad='#dddddd')

fig, axs = plt.subplots(2, 1, figsize=(13, 4.6), sharex=True)

ax = axs[0]
im = ax.imshow(np.ma.masked_invalid(region_af.values[:, cols]), cmap=cmap_af,
               aspect='auto', vmin=0, vmax=vmax_af)
ax.axvline(split_at, color='k', linewidth=1.8)
plu.format_ax(ax=ax, yticks=REGIONS, xticks=[], title='Region VAF')
ax.text(split_at / 2, -.9, f'non-enriched (n={int((~is_enriched).sum())})',
        ha='center', va='bottom', fontsize=7)
ax.text((split_at + len(cols)) / 2, -.9, f'enriched (n={int(is_enriched.sum())})',
        ha='center', va='bottom', fontsize=7)
cb = fig.colorbar(im, ax=ax, pad=.005, fraction=.015)
cb.set_label('VAF', fontsize=7)
cb.ax.tick_params(labelsize=6)

ax = axs[1]
im = ax.imshow(np.ma.masked_invalid(logp.values[:, cols]), cmap=cmap_p,
               aspect='auto', vmin=0, vmax=vmax_p)
ax.axvline(split_at, color='k', linewidth=1.8)
for yi in range(len(REGIONS)):
    for xi, j in enumerate(cols):
        if fisher_p.values[yi, j] < 0.05:
            ax.text(xi, yi, '*', ha='center', va='center', fontsize=7, color='white')
plu.format_ax(ax=ax, yticks=REGIONS, xticks=[],
              title='Region enrichment (-log10(pvalue), Fisher\'s exact test)')
ax.set_xticks(range(len(cols)))
ax.set_xticklabels([m.rsplit('_', 2)[0] for m in muts[cols]], rotation=90, fontsize=3.2)
cb = fig.colorbar(im, ax=ax, pad=.005, fraction=.015)
cb.set_label('-log10 p', fontsize=7)
cb.ax.tick_params(labelsize=6)

fig.subplots_adjust(left=.04, right=.96, top=.88, bottom=.18, hspace=.22)
fig.savefig(os.path.join(path_figures, 'lineage_enrichment.pdf'))
plt.show()


##


# 1c. Extended patches vs scattered, among the enriched mutations only
D_space = pairwise_distances(coords.values)
np.fill_diagonal(D_space, np.inf)
neighbours = np.argsort(D_space, axis=1)


def spatial_class(j, k):
    """
    Extended when at least half of the present cuts have another present cut
    among their k nearest neighbours in space.
    """
    idx = np.where(present[:, j])[0]
    if len(idx) < 2:
        return 'Scattered'
    share = np.mean([len(set(neighbours[i, :k]) & set(idx)) > 0 for i in idx])
    return 'Extended patches' if share >= .5 else 'Scattered'


enriched_idx = np.where(is_enriched)[0]
patch_class = pd.Series({muts[j]: spatial_class(j, K_NEIGHBOURS) for j in enriched_idx})
alt_counts = {
    k: pd.Series({muts[j]: spatial_class(j, k) for j in enriched_idx}).value_counts().to_dict()
    for k in (4, 5)
}
print('extended vs scattered, k=3:', patch_class.value_counts().to_dict(),
      '| k=4:', alt_counts[4], '| k=5:', alt_counts[5])

fig, ax = plt.subplots(figsize=(2.6, 3))
order = ['Extended patches', 'Scattered']
vals = patch_class.value_counts().reindex(order).fillna(0).astype(int)
ax.bar(range(2), vals.values, color=['#756bb1', '#bcbddc'], width=.6)
for i, v in enumerate(vals.values):
    ax.text(i, v + .3, str(v), ha='center', fontsize=8)
plu.format_ax(ax=ax, xticks=order, ylabel='n SNVs', rotx=20, reduced_spines=True,
              title=f'k={K_NEIGHBOURS} neighbours')
ax.text(.5, -.32, f'k=4: {alt_counts[4].get("Extended patches", 0)}/'
        f'{alt_counts[4].get("Scattered", 0)}   k=5: '
        f'{alt_counts[5].get("Extended patches", 0)}/{alt_counts[5].get("Scattered", 0)}',
        transform=ax.transAxes, ha='center', fontsize=6, color='#555555')
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'enrichment_type.pdf'))
plt.show()


##


# 1c, companion. One example of each pattern, in space
def pick_example(label):
    candidates = [j for j in enriched_idx if patch_class[muts[j]] == label]
    return max(candidates, key=lambda j: max_logp[j]) if candidates else None


fig = plt.figure(figsize=(8.5, 3))

ax = fig.add_subplot(1, 3, 1, projection='3d')
ax.computed_zorder = False
ax.scatter(coords['x'], coords['y'], coords['z'], s=POINT_SIZE,
           c=[COLORS[r] for r in reg], edgecolor='white', linewidth=.25,
           depthshade=False, zorder=5)
style_3d(ax, coords)
ax.set_title('Sampling', fontsize=8)
plu.add_legend(colors=COLORS, label='Region', ax=ax, ticks_size=5, artists_size=4,
               label_size=6, loc='upper center', bbox_to_anchor=(.5, .08), ncols=5)

for i, label in enumerate(order):
    j = pick_example(label)
    if j is None:
        continue
    ax = fig.add_subplot(1, 3, i + 2, projection='3d')
    ax.computed_zorder = False
    draw_mutation_3d(ax, j)
    ax.set_title(f'{label}\n{muts[j].rsplit("_", 2)[0]}', fontsize=8)

plu.add_cbar(np.array([0, .6]), ax=ax, label='Cell fraction', palette='afmhot_r',
             vmin=0, vmax=.6)
fig.text(.5, .02, 'grey = undetermined, open = absent, filled = present',
         ha='center', fontsize=6, color='#555555')
fig.subplots_adjust(left=.01, right=.9, top=.88, bottom=.08, wspace=.02)
fig.savefig(os.path.join(path_figures, 'example_enriched.pdf'))
plt.show()


##


# 1d. How far each enriched lineage swept, cut by cut
sweep = []
for j in enriched_idx:
    for i in np.where(present[:, j])[0]:
        sweep.append(dict(mutation=muts[j].rsplit('_', 2)[0], cell_fraction=2 * AF[i, j],
                          region=reg[i]))
sweep = pd.DataFrame(sweep)
sweep_order = (
    sweep.groupby('mutation')['cell_fraction'].mean().sort_values().index.tolist()
)

fig, ax = plt.subplots(figsize=(9, 3.2))
plu.strip(sweep, x='mutation', y='cell_fraction', x_order=sweep_order, size=4,
          color='#756bb1', ax=ax)
for i, m in enumerate(sweep_order):
    mu = sweep.loc[sweep['mutation'] == m, 'cell_fraction'].mean()
    ax.plot([i - .3, i + .3], [mu, mu], color='k', linewidth=1.1, zorder=5)
plu.format_ax(ax=ax, xlabel='Enriched lineages', ylabel='Cell fraction', rotx=90,
              xticks_size=4, reduced_spines=True,
              title='Sweep of enriched lineages (black line: mean)')
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'lineage_sweep.pdf'))
plt.show()


##


# 2a. Region dendrogram beside the region VAF matrix
profiles = np.vstack([AF[row_idx[r]].mean(0) for r in REGIONS])
D_reg = pairwise_distances(profiles, metric='cosine')
Z_reg = linkage(squareform(D_reg, checks=False), method='average')


def clade_set(Z, labels):
    """
    Groups of the tree as (leaf set, merge height), read off the linkage matrix.
    """
    n = len(labels)
    members = {i: {labels[i]} for i in range(n)}
    out = []
    for idx, (a, b, height, _) in enumerate(Z):
        joined = members[int(a)] | members[int(b)]
        members[n + idx] = joined
        if len(joined) < n:
            out.append((frozenset(joined), float(height)))
    return out


n_char = int(CHAR_FRACTION * len(muts))
support = {}
for _ in range(N_BOOT):
    cols_b = rng.choice(len(muts), n_char, replace=False)
    Pb = np.vstack([AF[np.ix_(row_idx[r], cols_b)].mean(0) for r in REGIONS])
    Zb = linkage(squareform(pairwise_distances(Pb, metric='cosine'), checks=False),
                 method='average')
    for grp, _h in clade_set(Zb, REGIONS):
        support[grp] = support.get(grp, 0) + 1

leaf_order = [REGIONS[i] for i in leaves_list(Z_reg)][::-1]   # top to bottom
dominant = region_af.idxmax(0).reindex(muts)
mut_order = np.concatenate([
    np.array([j for j in np.argsort(-max_logp) if dominant.iloc[j] == r])
    for r in leaf_order
]).astype(int)
block_edges = np.cumsum([sum(dominant.values == r) for r in leaf_order])[:-1] - .5

fig = plt.figure(figsize=(10, 2.8))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 3], wspace=.02)

ax = fig.add_subplot(gs[0, 0])
dn = dendrogram(Z_reg, labels=REGIONS, orientation='left', color_threshold=0,
                above_threshold_color='k', ax=ax)
for grp, height in clade_set(Z_reg, REGIONS):
    ys = [dn['ivl'].index(m) for m in grp]
    ax.text(height, 10 * np.mean(ys) + 5, f'{100 * support.get(grp, 0) / N_BOOT:.0f}% ',
            fontsize=6, va='center', ha='right', color='#c0392b')
plu.format_ax(ax=ax, xlabel='Cosine distance', reduced_spines=True)
ax.set_yticks([])
cax = ax.inset_axes((.05, -.28, .5, .06))
cb = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0, vmax_af), cmap='afmhot_r'),
                  cax=cax, orientation='horizontal')
cb.set_label('VAF', fontsize=6, labelpad=1)
cb.ax.tick_params(labelsize=5, length=2, pad=1)

ax = fig.add_subplot(gs[0, 1])
heat_rows = [REGIONS.index(r) for r in leaf_order]
ax.imshow(np.ma.masked_invalid(region_af.values[np.ix_(heat_rows, mut_order)]),
          cmap=cmap_af, aspect='auto', vmin=0, vmax=vmax_af)
for e in block_edges:
    ax.axvline(e, color='k', linewidth=.6)
ax.set_yticks(range(len(leaf_order)))
ax.set_yticklabels(leaf_order, fontsize=8)
ax.yaxis.tick_right()
plu.format_ax(ax=ax, xticks=[], xlabel=f'{len(muts)} mutations, grouped by the region '
                                       'carrying their highest VAF')

fig.subplots_adjust(left=.02, right=.95, top=.95, bottom=.3)
fig.savefig(os.path.join(path_figures, 'phylogenetic_relationships.pdf'))
plt.show()


##


# 2b. The same matrix at single-cut resolution
cut_order = []
for r in leaf_order:
    for c in sorted(set(chunk[row_idx[r]])):
        idx = np.array([i for i in row_idx[r] if chunk[i] == c])
        if len(idx) > 2:
            sub = AF[np.ix_(idx, mut_order)]
            idx = idx[leaves_list(linkage(squareform(
                pairwise_distances(sub, metric='cosine'), checks=False), method='average'))]
        cut_order.extend(idx.tolist())
cut_order = np.array(cut_order)
region_edges = np.cumsum([len(row_idx[r]) for r in leaf_order])[:-1] - .5
chunk_edges = np.where(chunk[cut_order][:-1] != chunk[cut_order][1:])[0] + .5

fig, ax = plt.subplots(figsize=(10, 5))
im = ax.imshow(AF[np.ix_(cut_order, mut_order)], cmap='afmhot_r', aspect='auto',
               vmin=0, vmax=np.nanpercentile(AF, 99))
for e in chunk_edges:
    ax.axhline(e, color='#999999', linewidth=.3)
for e in region_edges:
    ax.axhline(e, color='k', linewidth=1.2)
for e in block_edges:
    ax.axvline(e, color='k', linewidth=.6)
for k, i in enumerate(cut_order):
    ax.add_patch(plt.Rectangle((-2.5, k - .5), 2, 1, color=COLORS[reg[i]],
                               clip_on=False, linewidth=0))
plu.format_ax(ax=ax, xticks=[], yticks=[], xlabel=f'{len(muts)} mutations',
              ylabel='62 LCM cuts, by region then chunk')
plu.add_legend(colors=COLORS, label='Region', ax=ax, ticks_size=6, artists_size=5,
               label_size=7, loc='upper left', bbox_to_anchor=(0, 1.12), ncols=5)
cb = fig.colorbar(im, ax=ax, pad=.01, fraction=.02)
cb.set_label('VAF', fontsize=7)
cb.ax.tick_params(labelsize=6)

fig.subplots_adjust(left=.06, right=.94, top=.9, bottom=.08)
fig.savefig(os.path.join(path_figures, 'lcm_heatmap.pdf'))
plt.show()


##


# 2c. The region distances behind the dendrogram
fig, ax = plt.subplots(figsize=(3.6, 3))
D_show = pd.DataFrame(D_reg, index=REGIONS, columns=REGIONS).loc[leaf_order, leaf_order]
plu.plot_heatmap(D_show, ax=ax, annot=True, fmt='.3f', cb=True,
                 title='Region cosine distance')
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'distances.pdf'))
plt.show()


##


# 3. Septum/ventricle asymmetry: the statistic, its test, and two lineages
lean = lean_tab.set_index('Sample_ID')['lean_LV_minus_RV'].loc[cuts].values
septal = np.isin(reg, ['LS', 'CS', 'RS'])
observed = [lean[reg == r].mean() for r in ['LS', 'RS', 'CS']]
separation = min(observed[0], observed[1]) - observed[2]

null = np.empty(N_PERM_JOINT)
pattern = np.zeros(N_PERM_JOINT, dtype=bool)
septal_idx = np.where(septal)[0]
for t in range(N_PERM_JOINT):
    labels = reg.copy()
    labels[septal_idx] = rng.permutation(reg[septal_idx])
    means = [lean[labels == r].mean() for r in ['LS', 'RS', 'CS']]
    null[t] = min(means[0], means[1]) - means[2]
    pattern[t] = means[0] > 0 and means[1] > 0 and means[2] < 0
# The test counts a replicate only when it reproduces the sign pattern as well
p_joint = np.mean(pattern & (null >= separation))
p_sep = np.mean(null >= separation)

fig = plt.figure(figsize=(10, 5.6))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.15], hspace=.45, wspace=.3)

ax = fig.add_subplot(gs[0, 0])
for i, r in enumerate(['LS', 'RS', 'CS']):
    v = lean[reg == r]
    ax.scatter(i + np.random.uniform(-.14, .14, len(v)), v, s=18, color=COLORS[r],
               edgecolor='k', linewidth=.2, zorder=3)
    ax.plot([i - .28, i + .28], [v.mean()] * 2, color='k', linewidth=1.3, zorder=4)
ax.axhline(0, color='#c0392b', linestyle='--', linewidth=.8)
plu.format_ax(ax=ax, xticks=['LS', 'RS', 'CS'], ylabel='Similarity to LV minus RV',
              title='Ventricular lean, per cut', reduced_spines=True)

ax = fig.add_subplot(gs[0, 1])
bins = np.histogram_bin_edges(null, bins=60)
ax.hist([null[~pattern], null[pattern]], bins=bins, stacked=True,
        color=['#e0e0e0', '#7f7f7f'], edgecolor='white', linewidth=.2,
        label=['sign pattern fails', 'LS>0, RS>0, CS<0'])
ax.axvline(separation, color='#c0392b', linewidth=1.2)
ax.text(.97, .72, f'observed {separation:.3f}\njoint p = {p_joint:.4f}\n'
        f'separation alone {p_sep:.3f}', transform=ax.transAxes,
        color='#c0392b', fontsize=6, va='top', ha='right')
ax.legend(frameon=False, fontsize=5, loc='upper left')
plu.format_ax(ax=ax, xlabel='min(LS, RS) minus CS', ylabel='Permutations',
              title=f'Septal labels shuffled ({N_PERM_JOINT:,})', reduced_spines=True)

ax = fig.add_subplot(gs[0, 2])
loo_show = loo.loc[[i for i in loo.index if i != 'none']][['LS', 'RS', 'CS']]
im = ax.imshow(loo_show.values, cmap='RdBu', vmin=0, vmax=1, aspect='auto')
for yi in range(loo_show.shape[0]):
    for xi in range(loo_show.shape[1]):
        v = loo_show.values[yi, xi]
        if np.isfinite(v):
            ax.text(xi, yi, f'{v:.2f}', ha='center', va='center', fontsize=5)
ax.set_yticks(range(loo_show.shape[0]))
ax.set_yticklabels([i.replace('_Ventricle', 'V').replace('_septum', 'S')
                    for i in loo_show.index], fontsize=5)
plu.format_ax(ax=ax, xticks=list(loo_show.columns),
              title='Leave-one-chunk-out\n(fraction of draws leaning LV)')
cb = fig.colorbar(im, ax=ax, pad=.02, fraction=.03)
cb.ax.tick_params(labelsize=5)

examples = [('chr6_166168803_G_A', 'LS / LV lineage'),
            ('chr1_22366052_G_A', 'CS / RV lineage')]
for i, (mut, title) in enumerate(examples):
    if mut not in set(muts):
        continue
    ax = fig.add_subplot(gs[1, i], projection='3d')
    ax.computed_zorder = False
    j = list(muts).index(mut)
    g_af = summary['global_AF'][mut]
    patch = (binom.sf(AD_df.values[:, j] - 1, DP_df.values[:, j], g_af) < .01) & \
            (AD_df.values[:, j] >= 3) & present[:, j]
    draw_mutation_3d(ax, j, edge_by_region=True, ring=patch)
    ax.set_title(f'{title}\n{mut.rsplit("_", 2)[0]}', fontsize=8)

ax = fig.add_subplot(gs[1, 2])
ax.axis('off')
plu.add_legend(colors=COLORS, label='Region (point outline)', ax=ax, ticks_size=6,
               artists_size=5, label_size=7, loc='upper left', bbox_to_anchor=(0, .9))
plu.add_cbar(np.array([0, .6]), ax=ax, label='Cell fraction', palette='afmhot_r',
             vmin=0, vmax=.6)

fig.savefig(os.path.join(path_figures, 'septum_ventricle_asymmetry.pdf'))
plt.show()

print('written: embryonic_lineage_features, genotype_call_composition, '
      'lineage_enrichment, enrichment_type, example_enriched, lineage_sweep, '
      'phylogenetic_relationships, lcm_heatmap, distances, septum_ventricle_asymmetry')
