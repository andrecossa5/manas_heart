"""
Septum vs ventricles: topology robustness, 3D view of the marker SNVs, and matching of septal cuts
to ventricle chunk / region pseudobulks (NEW_ANALYSIS.md).

1. Topology of Septum (LS+CS+RS) / LV / RV from pooled VAF profiles (cosine, average linkage),
   resampling SNVs, cuts within group, or both.
2. 3D read-count view of every SNV in markers_exclusive.pdf.
3. Chunks in 3D; each septal cut matched to the 6 ventricle chunks and the 2 ventricles.

Outputs: results/SEPTUM_VENTRICLE_TOPOLOGY.tsv, SEPTUM_MATCHING.tsv, SEPTUM_MATCHING_TESTS.tsv;
figures/septum_ventricle_topology.pdf, markers_exclusive_3d.pdf, chunks_3d.pdf, septum_matching.pdf.
"""

import os
import itertools
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr, binom, wilcoxon
from sklearn.metrics import pairwise_distances
import plotting_utils as plu

matplotlib.use('Agg')
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
CLASS_COLORS = {'Pre-gastrulation': '#9e9ac8', 'Other shared': '#d4a72c', 'Heart-specific': '#31a354'}
GROUP_COLORS = {'Septum': '#54A24B', 'LV': COLORS['LV'], 'RV': COLORS['RV']}
GROUPS = ['Septum', 'LV', 'RV']
N_BOOT = 1000
N_PERM_CUT = 10000
MARKER_FRAC = .5
POINT_SIZE = 22
VAF_MAX_3D = .2


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


def clade_set(Z, labels):
    """
    Groups of an average-linkage tree, as (leaf set, merge height) pairs.
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


def cosine_dist(P):
    return pairwise_distances(P, metric='cosine')


def tree_of(P):
    return linkage(squareform(cosine_dist(P), checks=False), method='average')


##


path_main = '/Users/cossa/Desktop/projects/manas_heart'
path_results = os.path.join(path_main, 'results')
path_figures = os.path.join(path_main, 'figures')

geno = pd.read_csv(os.path.join(path_results, 'GENOTYPES_TRUE.tsv.gz'), sep='\t')
summary = pd.read_csv(os.path.join(path_results, 'LINEAGE_SUMMARY.tsv'), sep='\t').set_index('mutation_id')
enrich = pd.read_csv(os.path.join(path_results, 'REGION_ENRICHMENT.tsv'), sep='\t')
xyz = pd.read_csv(os.path.join(path_main, 'data', 'Heart_final_coorindates_135.csv')).set_index('name')

geno['reg'] = geno['region'].map(REGION_ABBR)
AD_df = geno.pivot(index='Sample_ID', columns='mutation_id', values='AD_alt')
DP_df = geno.pivot(index='Sample_ID', columns='mutation_id', values='DP')
cuts, muts = AD_df.index, AD_df.columns
A, D = AD_df.values.astype(float), DP_df.values.astype(float)
AF = A / D
meta = geno.drop_duplicates('Sample_ID').set_index('Sample_ID').loc[cuts]
reg, chunk = meta['reg'].values, meta['chunk'].values
coords = xyz.loc[cuts, ['x', 'y', 'z']]
klass = summary['lineage_class'].reindex(muts).values
group = np.where(np.isin(reg, ['LS', 'CS', 'RS']), 'Septum', reg)
group_idx = {g: np.where(group == g)[0] for g in GROUPS}
print(f'{len(muts)} SNVs x {len(cuts)} cuts | ' +
      ' | '.join(f'{g} {len(i)} cuts, {len(set(chunk[i]))} blocks' for g, i in group_idx.items()))


def pooled(rows, cols):
    """
    Pooled VAF, sum AD / sum DP over the given cuts, for the given SNVs.
    """
    return A[np.ix_(rows, cols)].sum(0) / D[np.ix_(rows, cols)].sum(0)


##


# 1. Topology of Septum / LV / RV. With three leaves the topology is the pair joined first.
all_snvs = np.arange(len(muts))
PAIRS = [('Septum', 'LV'), ('Septum', 'RV'), ('LV', 'RV')]
PAIR_IDX = [(0, 1), (0, 2), (1, 2)]


def profiles(rows_by_group, cols):
    return np.vstack([pooled(rows_by_group[g], cols) for g in GROUPS])


def topology(P):
    d = cosine_dist(P)
    dist = np.array([d[a, b] for a, b in PAIR_IDX])
    return int(dist.argmin()), dist


# Sanity: pooled Septum profile equals the sum of the regional pooled counts
reg_sum_A = sum(A[reg == r].sum(0) for r in ['LS', 'CS', 'RS'])
reg_sum_D = sum(D[reg == r].sum(0) for r in ['LS', 'CS', 'RS'])
assert np.allclose(pooled(group_idx['Septum'], all_snvs), reg_sum_A / reg_sum_D)

obs_top, obs_dist = topology(profiles(group_idx, all_snvs))
print('\n1 Septum / LV / RV topology (pooled VAF, cosine, average linkage)')
print('  observed distances: ' + ', '.join(f'{a}-{b} {d:.3f}' for (a, b), d in zip(PAIRS, obs_dist)))
print(f'  observed first merge: {"-".join(PAIRS[obs_top])}')

rng = np.random.default_rng(10)
schemes = ['SNVs', 'cuts', 'SNVs + cuts']
counts = {s: np.zeros(3) for s in schemes}
dists = {s: [] for s in schemes}
for _ in range(N_BOOT):
    cols_b = rng.choice(all_snvs, len(all_snvs), replace=True)
    rows_b = {g: rng.choice(i, len(i), replace=True) for g, i in group_idx.items()}
    for s, rows, cols in [('SNVs', group_idx, cols_b), ('cuts', rows_b, all_snvs),
                          ('SNVs + cuts', rows_b, cols_b)]:
        t, dist = topology(profiles(rows, cols))
        counts[s][t] += 1
        dists[s].append(dist)
topo = pd.DataFrame([dict(scheme=s, first_merge='-'.join(PAIRS[k]), frequency=counts[s][k] / N_BOOT,
                          observed=k == obs_top) for s in schemes for k in range(3)])
topo.to_csv(os.path.join(path_results, 'SEPTUM_VENTRICLE_TOPOLOGY.tsv'), sep='\t', index=False)
print(f'  first-merge frequency over {N_BOOT} resamplings (cuts are resampled as independent, so '
      f'cut support overstates block-level support):')
print(topo.pivot(index='first_merge', columns='scheme', values='frequency')[schemes].round(3).to_string())

fig = plt.figure(figsize=(10, 2.9))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.3], wspace=.7)
ax = fig.add_subplot(gs[0, 0])
from scipy.cluster.hierarchy import dendrogram
Z_obs = tree_of(profiles(group_idx, all_snvs))
dendrogram(Z_obs, labels=GROUPS, orientation='left', color_threshold=0, above_threshold_color='k', ax=ax)
for grp, height in clade_set(Z_obs, GROUPS):
    ax.text(height, 20, f'{height:.3f}', fontsize=6, va='bottom', ha='center', color='#555555')
plu.format_ax(ax=ax, xlabel='Cosine distance', title='Observed tree', reduced_spines=True)
ax = fig.add_subplot(gs[0, 1])
w = .26
for i, s in enumerate(schemes):
    ax.bar(np.arange(3) + (i - 1) * w, counts[s] / N_BOOT, width=w, label=s,
           color=['#bdbdbd', '#737373', '#252525'][i])
ax.set_xticks(range(3))
ax.set_xticklabels(['-'.join(p) for p in PAIRS], rotation=20, fontsize=6)
ax.legend(frameon=False, fontsize=5.5)
plu.format_ax(ax=ax, ylabel='Frequency first merge', title='Topology support', reduced_spines=True)
ax = fig.add_subplot(gs[0, 2])
arr = np.array(dists['SNVs + cuts'])
bp = ax.boxplot([arr[:, k] for k in range(3)], showfliers=False, patch_artist=True,
                medianprops=dict(color='k'))
for p, c in zip(bp['boxes'], ['#9ecae1', '#fdae6b', '#c7e9c0']):
    p.set(facecolor=c)
ax.scatter(range(1, 4), obs_dist, color='#c0392b', zorder=5, s=14, label='observed')
ax.set_xticklabels(['-'.join(p) for p in PAIRS], fontsize=6)
ax.legend(frameon=False, fontsize=5.5)
plu.format_ax(ax=ax, ylabel='Cosine distance', title='SNVs + cuts resampled', reduced_spines=True)
fig.subplots_adjust(left=.07, right=.97, top=.82, bottom=.22)
fig.savefig(os.path.join(path_figures, 'septum_ventricle_topology.pdf'))


##


# 2. The SNVs of markers_exclusive.pdf in 3D. Selection as in 9.figures.py: region tree on the
# pooled regional VAF of all SNVs; a SNV is a marker if its carrier regions (VAF >= MARKER_FRAC
# of its maximum) are exactly one clade of that tree, exclusive if exactly one region.
region_af = enrich.pivot_table(index='region', columns='mutation_id', values='AF').reindex(REGIONS)[muts]
Z_reg = tree_of(region_af.values)
tree_clades = [c for c, _h in clade_set(Z_reg, REGIONS)]
carriers = [frozenset(np.array(REGIONS)[region_af.values[:, j] >= MARKER_FRAC * np.nanmax(region_af.values[:, j])])
            for j in range(len(muts))]
marker_groups = []
for c in sorted(tree_clades, key=len):
    sel = [j for j in range(len(muts)) if carriers[j] == c]
    if sel:
        marker_groups.append(('+'.join(r for r in REGIONS if r in c), sel))
exclusive_groups = [(r, [j for j in range(len(muts)) if carriers[j] == frozenset({r})]) for r in REGIONS]
selected = [(lab, j) for lab, g in marker_groups + [x for x in exclusive_groups if x[1]] for j in g]
print(f'\n2 markers_exclusive SNVs: {len(selected)} | ' + ', '.join(
    f'{lab} {len(g)}' for lab, g in marker_groups + exclusive_groups if g))
assert len(selected) == 33, 'marker selection no longer matches markers_exclusive.pdf'


def draw_counts_3d(ax, j):
    """
    Cuts with >=1 alternate read filled by VAF, others open.
    """
    carries = AD_df.values[:, j] > 0
    ax.scatter(coords['x'][~carries], coords['y'][~carries], coords['z'][~carries], s=POINT_SIZE,
               facecolor='white', edgecolor='k', linewidth=.3, depthshade=False, zorder=4)
    ax.scatter(coords['x'][carries], coords['y'][carries], coords['z'][carries], s=POINT_SIZE,
               c=AF[carries, j], cmap='afmhot_r', vmin=0, vmax=VAF_MAX_3D,
               edgecolor='k', linewidth=.3, depthshade=False, zorder=5)
    style_3d(ax, coords)


ncol = 6
nrow = int(np.ceil((len(selected) + 1) / ncol))
fig = plt.figure(figsize=(2.6 * ncol, 2.6 * nrow))
ax = fig.add_subplot(nrow, ncol, 1, projection='3d')
ax.computed_zorder = False
ax.scatter(coords['x'], coords['y'], coords['z'], s=POINT_SIZE, c=[COLORS[r] for r in reg],
           edgecolor='white', linewidth=.25, depthshade=False, zorder=5)
style_3d(ax, coords)
ax.set_title('Sampling', fontsize=7)
plu.add_legend(colors=COLORS, label='Region', ax=ax, ticks_size=5, artists_size=4,
               label_size=6, loc='upper center', bbox_to_anchor=(.5, .08), ncols=5)
for k, (lab, j) in enumerate(selected):
    ax = fig.add_subplot(nrow, ncol, k + 2, projection='3d')
    ax.computed_zorder = False
    draw_counts_3d(ax, j)
    ax.set_title(f'{lab}\n{muts[j].rsplit("_", 2)[0].replace("_", ":")}\n{klass[j]}', fontsize=6)
plu.add_cbar(np.array([0, VAF_MAX_3D]), ax=ax, label='VAF', palette='afmhot_r', vmin=0, vmax=VAF_MAX_3D)
fig.text(.5, .005, 'filled = cut with >=1 alternate read (colour: VAF), open = none', ha='center',
         fontsize=6, color='#555555')
fig.subplots_adjust(left=.01, right=.94, top=.95, bottom=.02, wspace=.02, hspace=.2)
fig.savefig(os.path.join(path_figures, 'markers_exclusive_3d.pdf'))


##


# 3. Chunks in space, and septal cuts matched to ventricle chunk / region pseudobulks
chunks = sorted(set(chunk))
chunk_region = {c: reg[chunk == c][0] for c in chunks}
chunk_rows = {c: np.where(chunk == c)[0] for c in chunks}
centroid = pd.DataFrame({c: coords.values[chunk_rows[c]].mean(0) for c in chunks}, index=['x', 'y', 'z']).T

# 14 chunks: each region's colour, lighter for later chunks
chunk_color = {}
for r in REGIONS:
    cs = [c for c in chunks if chunk_region[c] == r]
    for i, c in enumerate(cs):
        chunk_color[c] = matplotlib.colors.to_hex(
            np.array(matplotlib.colors.to_rgb(COLORS[r])) * (1 - .55 * i / max(len(cs), 1)) + .55 * i / max(len(cs), 1))

fig = plt.figure(figsize=(5.5, 4.5))
ax = fig.add_subplot(projection='3d')
ax.computed_zorder = False
ax.scatter(coords['x'], coords['y'], coords['z'], s=POINT_SIZE, c=[chunk_color[c] for c in chunk],
           edgecolor='k', linewidth=.2, depthshade=False, zorder=5)
for c in chunks:
    ax.text(*centroid.loc[c].values, c.replace('_Ventricle', ' V').replace('_septum', ' S'), fontsize=5,
            ha='center', zorder=10)
style_3d(ax, coords)
ax.set_title('62 LCM cuts, coloured by chunk', fontsize=8)
fig.subplots_adjust(left=0, right=1, top=.93, bottom=0)
fig.savefig(os.path.join(path_figures, 'chunks_3d.pdf'))

vent_chunks = [c for c in chunks if chunk_region[c] in ('LV', 'RV')]
septal = np.where(group == 'Septum')[0]
chunk_P = np.vstack([pooled(chunk_rows[c], all_snvs) for c in vent_chunks])
region_P = np.vstack([pooled(group_idx[g], all_snvs) for g in ['LV', 'RV']])
d_chunk = pairwise_distances(AF[septal], chunk_P, metric='cosine')            # cuts x 6 chunks
d_region = pairwise_distances(AF[septal], region_P, metric='cosine')          # cuts x (LV, RV)
sp_dist = pairwise_distances(coords.values[septal], centroid.loc[vent_chunks].values)
best_chunk = d_chunk.argmin(1)
nearest = sp_dist.argmin(1)
best_lv = d_region[:, 0] < d_region[:, 1]            # best region is LV
chunk_is_lv = np.array([chunk_region[c] == 'LV' for c in vent_chunks])
sreg, schunk = reg[septal], chunk[septal]

match = pd.DataFrame({
    'Sample_ID': cuts[septal], 'region': sreg, 'chunk': schunk,
    **{f'd_{c}': d_chunk[:, k] for k, c in enumerate(vent_chunks)},
    'd_LV': d_region[:, 0], 'd_RV': d_region[:, 1],
    'best_chunk': [vent_chunks[k] for k in best_chunk],
    'best_region': np.where(best_lv, 'LV', 'RV'),
    'nearest_chunk': [vent_chunks[k] for k in nearest],
    'best_is_nearest': best_chunk == nearest,
})
match.to_csv(os.path.join(path_results, 'SEPTUM_MATCHING.tsv'), sep='\t', index=False)

print('\n3 Septal cuts matched to ventricle pseudobulks (cosine on 123 SNVs)')
print('  best region, by septal region (cuts):')
print(pd.crosstab(match['region'], match['best_region']).to_string().replace('\n', '\n    '))
print('  best chunk, by septal region:')
print(pd.crosstab(match['region'], match['best_chunk']).to_string().replace('\n', '\n    '))
print(f'  best chunk is the spatially nearest ventricle chunk: {int(match.best_is_nearest.sum())} of {len(match)}')

# (a) Side correspondence: LS -> LV and RS -> RV among lateral cuts. Null: LS / RS labels
# shuffled over cuts, and exactly over the lateral blocks (blocks are the sampling units).
lat = np.where(np.isin(sreg, ['LS', 'RS']))[0]
lat_blocks = sorted(set(schunk[lat]))


def side_stats(is_ls):
    """
    Share of lateral cuts whose best ventricle is on their own side, and the region-level lean
    (LS minus RS of d(RV) - d(LV)) from pooled profiles of the labelled cuts.
    """
    own = np.where(is_ls, best_lv[lat], ~best_lv[lat]).mean()
    lean = []
    for sel in (is_ls, ~is_ls):
        p = pooled(septal[lat][sel], all_snvs)[None, :]
        d_ = pairwise_distances(p, region_P, metric='cosine')[0]
        lean.append(d_[1] - d_[0])
    return own, lean[0] - lean[1]


obs_is_ls = sreg[lat] == 'LS'
obs_own, obs_lean = side_stats(obs_is_ls)
rng = np.random.default_rng(3)
null_cut = np.array([side_stats(rng.permutation(obs_is_ls)) for _ in range(N_PERM_CUT)])
n_ls_blocks = len({b for b in lat_blocks if sreg[lat][schunk[lat] == b][0] == 'LS'})
null_blk = []
for pick in itertools.combinations(lat_blocks, n_ls_blocks):
    null_blk.append(side_stats(np.isin(schunk[lat], pick)))
null_blk = np.array(null_blk)


def upper_p(null, obs):
    return (np.sum(null >= obs - 1e-12) + 1) / (len(null) + 1)


print(f'\n  (a) side correspondence, {len(lat)} lateral cuts in {len(lat_blocks)} blocks')
print(f'    share on own side {obs_own:.3f} | cut-label permutation p = {upper_p(null_cut[:, 0], obs_own):.4f} '
      f'| block-label exact ({len(null_blk)} labellings) p = {np.mean(null_blk[:, 0] >= obs_own - 1e-12):.3f}')
print(f'    region-level lean LS - RS (d_RV - d_LV) {obs_lean:+.4f} | cut-label p = '
      f'{upper_p(null_cut[:, 1], obs_lean):.4f} | block-label exact p = '
      f'{np.mean(null_blk[:, 1] >= obs_lean - 1e-12):.3f}')
print('    pooled region lean (d_RV - d_LV): ' + ', '.join(
    f'{r} {(lambda d_: d_[1] - d_[0])(pairwise_distances(pooled(np.where(reg == r)[0], all_snvs)[None, :], region_P, metric="cosine")[0]):+.3f}'
    for r in ['LS', 'CS', 'RS']))

# (b) Proximity: does profile distance track spatial distance? Exact null: the 6 chunk profiles
# permuted over the 6 chunk positions (720).
def proximity_stats(perm):
    dc = d_chunk[:, perm]
    return (dc.argmin(1) == nearest).mean(), spearmanr(dc.ravel(), sp_dist.ravel())[0]


obs_prox = proximity_stats(np.arange(len(vent_chunks)))
null_prox = np.array([proximity_stats(np.array(p)) for p in itertools.permutations(range(len(vent_chunks)))])
p_frac = np.mean(null_prox[:, 0] >= obs_prox[0] - 1e-12)
p_rho = np.mean(null_prox[:, 1] <= obs_prox[1] + 1e-12)
print(f'\n  (b) proximity: best chunk = nearest chunk in {obs_prox[0]:.3f} of cuts '
      f'(null mean {null_prox[:, 0].mean():.3f}, exact p = {p_frac:.3f}) | Spearman '
      f'(profile distance, spatial distance) {obs_prox[1]:+.3f} (null mean {null_prox[:, 1].mean():+.3f}, '
      f'exact one-sided p = {p_rho:.3f}; {len(null_prox)} permutations)')
print('  note: single-cut profiles are sparse (median DP 25) and LV has 4 chunks vs 2 for RV, so '
      'best-chunk counts favour LV by chance; the permutation nulls and the region-level lean account for it')

pd.DataFrame([
    dict(test='side: share of lateral cuts on own side', observed=obs_own,
         p_cut_permutation=upper_p(null_cut[:, 0], obs_own), p_block_exact=np.mean(null_blk[:, 0] >= obs_own - 1e-12)),
    dict(test='side: region lean LS minus RS', observed=obs_lean,
         p_cut_permutation=upper_p(null_cut[:, 1], obs_lean), p_block_exact=np.mean(null_blk[:, 1] >= obs_lean - 1e-12)),
    dict(test='proximity: best chunk is nearest chunk', observed=obs_prox[0],
         p_cut_permutation=np.nan, p_block_exact=p_frac),
    dict(test='proximity: Spearman profile vs spatial distance', observed=obs_prox[1],
         p_cut_permutation=np.nan, p_block_exact=p_rho),
]).to_csv(os.path.join(path_results, 'SEPTUM_MATCHING_TESTS.tsv'), sep='\t', index=False)
# p_block_exact: exact over block labellings (side tests) or over the 720 chunk-profile permutations (proximity)

# Figure
vcol = dict(zip(vent_chunks, plt.cm.tab10(np.arange(len(vent_chunks)))))
order = np.lexsort((schunk, np.array([REGIONS.index(r) for r in sreg])))
short = lambda c: c.replace('_Ventricle', ' V').replace('_septum', ' S')

fig = plt.figure(figsize=(11, 4.2))
gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.1, 1], wspace=.3)
ax = fig.add_subplot(gs[0, 0])
im = ax.imshow(d_chunk[order], cmap='viridis_r', aspect='auto')
for y, i in enumerate(order):
    ax.plot(best_chunk[i], y, 'o', color='white', markersize=3, markeredgecolor='k', markeredgewidth=.4)
    ax.add_patch(plt.Rectangle((nearest[i] - .5, y - .5), 1, 1, fill=False, edgecolor='#c0392b', linewidth=.8))
    ax.add_patch(plt.Rectangle((-1.2, y - .5), .5, 1, color=COLORS[sreg[i]], clip_on=False, linewidth=0))
ax.set_yticks([])
ax.set_xticks(range(len(vent_chunks)))
ax.set_xticklabels([short(c) for c in vent_chunks], rotation=45, ha='right', rotation_mode='anchor', fontsize=6)
ax.set_title('Septal cuts vs ventricle chunks\n(dot: best profile, red box: nearest in space)', fontsize=7)
cb = fig.colorbar(im, ax=ax, pad=.02, fraction=.04)
cb.set_label('Cosine distance', fontsize=6)
cb.ax.tick_params(labelsize=5)

ax = fig.add_subplot(gs[0, 1], projection='3d')
ax.computed_zorder = False
vi = np.where(np.isin(group, ['LV', 'RV']))[0]
ax.scatter(coords['x'].values[vi], coords['y'].values[vi], coords['z'].values[vi], s=POINT_SIZE,
           c=[vcol[c] if c in vcol else 'grey' for c in chunk[vi]], edgecolor='k', linewidth=.2,
           depthshade=False, zorder=4)
ax.scatter(coords['x'].values[septal], coords['y'].values[septal], coords['z'].values[septal],
           s=POINT_SIZE * 1.6, c=[vcol[vent_chunks[k]] for k in best_chunk], marker='^', edgecolor='k',
           linewidth=.4, depthshade=False, zorder=5)
style_3d(ax, coords)
ax.set_title('Circles: ventricle cuts; triangles: septal cuts,\ncoloured by best-matching chunk', fontsize=7)
plu.add_legend(colors={short(c): vcol[c] for c in vent_chunks}, label='Ventricle chunk', ax=ax,
               ticks_size=5, artists_size=4, label_size=6, loc='upper center', bbox_to_anchor=(.5, .06), ncols=3)

ax = fig.add_subplot(gs[0, 2])
ax.hist(null_cut[:, 0], bins=30, color='#bdbdbd', edgecolor='white', linewidth=.2)
ax.axvline(obs_own, color='#c0392b')
ax.text(.03, .95, f'side match {obs_own:.2f}\ncut perm p = {upper_p(null_cut[:, 0], obs_own):.3f}\n'
        f'block exact p = {np.mean(null_blk[:, 0] >= obs_own - 1e-12):.3f}\n'
        f'proximity {obs_prox[0]:.2f}, exact p = {p_frac:.3f}\n'
        f'Spearman {obs_prox[1]:+.2f}, exact p = {p_rho:.3f}',
        transform=ax.transAxes, fontsize=6, va='top', color='#c0392b')
plu.format_ax(ax=ax, xlabel='Lateral cuts on own side\n(LS/RS labels shuffled)', ylabel='Permutations',
              reduced_spines=True)
fig.subplots_adjust(left=.06, right=.97, top=.85, bottom=.2)
fig.savefig(os.path.join(path_figures, 'septum_matching.pdf'))


##


# 3b. Per-cut test, each septal cut against LV and RV. Two scores per cut:
#   lean = cosine d(RV) - d(LV) of the cut's VAF profile to the pooled ventricle profiles;
#   LLR  = binomial log-likelihood of the cut's AD/DP under the LV VAF profile minus under RV.
# Positive = LV. Uncertainty: SNVs resampled with replacement and ventricle cuts resampled within
# LV and within RV (profiles rebuilt each time); a cut is called when its sign holds in >= CALL_FRAC
# of N_BOOT_CUT replicates. Cut counts themselves are fixed, so per-cut sampling noise enters only
# through the SNV resampling.
N_BOOT_CUT = 1000
CALL_FRAC = .95
LL_FLOOR = 1e-4
A_s, D_s, AF_s = A[septal], D[septal], AF[septal]


def lean_scores(lv_rows, rv_rows, cols):
    P = np.vstack([pooled(lv_rows, cols), pooled(rv_rows, cols)])
    d_ = pairwise_distances(AF_s[:, cols], P, metric='cosine')
    p_ = np.clip(P, LL_FLOOR, 1 - LL_FLOOR)
    ll = [binom.logpmf(A_s[:, cols], D_s[:, cols], p_[k][None, :]).sum(1) for k in range(2)]
    return d_[:, 1] - d_[:, 0], ll[0] - ll[1]


lean_obs, llr_obs = lean_scores(group_idx['LV'], group_idx['RV'], all_snvs)
rng = np.random.default_rng(5)
boot = [lean_scores(rng.choice(group_idx['LV'], len(group_idx['LV']), replace=True),
                    rng.choice(group_idx['RV'], len(group_idx['RV']), replace=True),
                    rng.choice(all_snvs, len(all_snvs), replace=True)) for _ in range(N_BOOT_CUT)]
lean_b, llr_b = np.array([b[0] for b in boot]), np.array([b[1] for b in boot])
frac_lv = (lean_b > 0).mean(0)
frac_lv_llr = (llr_b > 0).mean(0)
call = np.where(frac_lv >= CALL_FRAC, 'LV', np.where(frac_lv <= 1 - CALL_FRAC, 'RV', 'undecided'))
call_llr = np.where(frac_lv_llr >= CALL_FRAC, 'LV', np.where(frac_lv_llr <= 1 - CALL_FRAC, 'RV', 'undecided'))

cut_lean = pd.DataFrame({
    'Sample_ID': cuts[septal], 'region': sreg, 'chunk': schunk, 'DP_median': np.median(D_s, 1),
    'lean_LV_positive': lean_obs, 'lean_lo': np.percentile(lean_b, 2.5, 0), 'lean_hi': np.percentile(lean_b, 97.5, 0),
    'frac_boot_LV': frac_lv, 'call': call,
    'LLR_LV_minus_RV': llr_obs, 'frac_boot_LV_LLR': frac_lv_llr, 'call_LLR': call_llr,
})
cut_lean.to_csv(os.path.join(path_results, 'SEPTUM_CUT_LEAN.tsv'), sep='\t', index=False)

print(f'\n3b Each septal cut vs LV and RV (call: sign held in >= {CALL_FRAC:.0%} of {N_BOOT_CUT} replicates; '
      f'positive = LV)')
for name, col in [('cosine lean', 'call'), ('binomial LLR', 'call_LLR')]:
    print(f'  {name}:')
    print(pd.crosstab(cut_lean['region'], cut_lean[col]).reindex(['LS', 'CS', 'RS']).reindex(
        columns=['LV', 'RV', 'undecided']).fillna(0).astype(int).to_string().replace('\n', '\n    '))
print(f'  cosine lean and LLR agree in sign for {int((np.sign(lean_obs) == np.sign(llr_obs)).sum())} of {len(septal)} cuts '
      f'| Spearman {spearmanr(lean_obs, llr_obs)[0]:+.2f}')
dec = cut_lean.query('region in ["LS", "RS"] and call != "undecided"')
print(f'  decided lateral cuts: {len(dec)} of {len(lat)} | on own side '
      f'{int(((dec.region == "LS") & (dec.call == "LV")).sum() + ((dec.region == "RS") & (dec.call == "RV")).sum())}')
print('  per block (mean lean, cuts called LV/RV): ' + '; '.join(
    f'{b.replace("_septum", "").replace("_", " ")} {g.lean_LV_positive.mean():+.3f} ({(g.call == "LV").sum()}/{(g.call == "RV").sum()})'
    for b, g in cut_lean.groupby('chunk')))

order_c = np.lexsort((schunk, np.array([REGIONS.index(r) for r in sreg])))
fig, axs = plt.subplots(1, 2, figsize=(10, 5.2), sharey=True)
y = np.arange(len(septal))
c_call = {'LV': COLORS['LV'], 'RV': COLORS['RV'], 'undecided': '#bdbdbd'}
ax = axs[0]
for k, i in enumerate(order_c):
    ax.plot([cut_lean.lean_lo[i], cut_lean.lean_hi[i]], [k, k], color=c_call[call[i]], linewidth=1)
    ax.scatter(lean_obs[i], k, s=14, color=c_call[call[i]], edgecolor='k', linewidth=.3, zorder=3)
    ax.add_patch(plt.Rectangle((-.34, k - .5), .015, 1, color=COLORS[sreg[i]], clip_on=False, linewidth=0))
ax.axvline(0, color='k', linewidth=.6, linestyle='--')
for e in np.where(sreg[order_c][:-1] != sreg[order_c][1:])[0] + .5:
    ax.axhline(e, color='k', linewidth=.8)
ax.set_xlim(-.33, .33)
ax.set_yticks(y)
ax.set_yticklabels([f'{sreg[i]} {schunk[i].split(".")[-1] if "." in schunk[i] else "1"}' for i in order_c], fontsize=5)
plu.format_ax(ax=ax, xlabel='Cosine lean, d(RV) - d(LV) (95% CI)', title='Per cut, <- RV   LV ->', reduced_spines=True)
ax = axs[1]
for k, i in enumerate(order_c):
    ax.barh(k, llr_obs[i], color=c_call[call_llr[i]], height=.8)
ax.axvline(0, color='k', linewidth=.6, linestyle='--')
plu.format_ax(ax=ax, xlabel='Binomial log-likelihood ratio, LV vs RV', title='Same cuts, read-count likelihood', reduced_spines=True)
ax.invert_yaxis()
plu.add_legend(colors=c_call, label=f'Sign held in >= {CALL_FRAC:.0%} of resamplings', ax=ax, ticks_size=6,
               artists_size=5, label_size=6, loc='lower right', bbox_to_anchor=(1, .02), ncols=1)
fig.subplots_adjust(left=.07, right=.97, top=.92, bottom=.1, wspace=.08)
fig.savefig(os.path.join(path_figures, 'septum_cut_lean.pdf'))


##


# 3c. Is the LV lean a size / depth effect? LV (18 cuts, 4 blocks, deep) is compared with RV (8 cuts,
# 2 blocks, shallow), so its pooled profile is less noisy and closer to the heart-wide average.
# Lean = d(RV) - d(LV) for a septal group's pooled profile, positive = closer to LV, with LV matched
# to RV in one of two ways:
#   blocks: each of the 6 pairs of LV blocks (RV has 2 blocks);
#   depth:  8 random LV cuts (RV has 8), alt/total counts thinned per SNV to RV's pooled depth.
# In every replicate the septal blocks are also resampled with replacement within the group.
# Run on all SNVs, without the SNVs exclusive to LV or RV (their carriers are one ventricle only),
# and on each lineage class.
N_DRAW = 2000
N_BLK_BOOT = 300
rv_rows, lv_rows = group_idx['RV'], group_idx['LV']
lv_blocks = sorted(set(chunk[lv_rows]))
excl_ventricle = {j for r, g in exclusive_groups if r in ('LV', 'RV') for j in g}
snv_sets = {
    'All SNVs': all_snvs,
    'No LV/RV-exclusive': np.array([j for j in all_snvs if j not in excl_ventricle]),
    'Pre-gastrulation': np.where(klass == 'Pre-gastrulation')[0],
    'Heart-specific': np.where(klass == 'Heart-specific')[0],
}
sept_groups = {'Septum': septal, 'LS': np.where(reg == 'LS')[0], 'RS': np.where(reg == 'RS')[0],
               'CS': np.where(reg == 'CS')[0]}


def resample_blocks(rows, rng_):
    bl = sorted(set(chunk[rows]))
    pick = rng_.choice(bl, len(bl), replace=True)
    return np.concatenate([rows[chunk[rows] == b] for b in pick])


def lean_vs(rows, lv_prof, rv_prof, cols):
    p = pooled(rows, cols)[None, :]
    d_ = pairwise_distances(p, np.vstack([lv_prof, rv_prof]), metric='cosine')[0]
    return d_[1] - d_[0]


rng = np.random.default_rng(7)
dp_rv_all = D[rv_rows].sum(0)
rows_out, draws_out = [], {}
for set_name, cols in snv_sets.items():
    rv_prof = pooled(rv_rows, cols)
    for gname, rows in sept_groups.items():
        raw = lean_vs(rows, pooled(lv_rows, cols), rv_prof, cols)
        by_block, by_depth = [], []
        for pair in itertools.combinations(lv_blocks, 2):
            lv_prof = pooled(lv_rows[np.isin(chunk[lv_rows], pair)], cols)
            by_block += [lean_vs(resample_blocks(rows, rng), lv_prof, rv_prof, cols) for _ in range(N_BLK_BOOT)]
        for _ in range(N_DRAW):
            r = rng.choice(lv_rows, len(rv_rows), replace=False)
            a_, d_ = A[r].sum(0), D[r].sum(0)
            d2 = np.minimum(dp_rv_all, d_).astype(int)
            a2 = rng.hypergeometric(a_.astype(int), (d_ - a_).astype(int), d2)
            lv_prof = (a2 / np.maximum(d2, 1))[cols]
            by_depth.append(lean_vs(resample_blocks(rows, rng), lv_prof, rv_prof, cols))
        for scheme, v in [('blocks', by_block), ('depth', by_depth)]:
            v = np.array(v)
            draws_out[(set_name, gname, scheme)] = (v, raw)
            rows_out.append(dict(snv_set=set_name, n_snvs=len(cols), septal_group=gname, scheme=scheme,
                                 raw_lean=raw, matched_median=np.median(v), lo95=np.percentile(v, 2.5),
                                 hi95=np.percentile(v, 97.5), frac_LV_closer=(v > 0).mean()))
matched = pd.DataFrame(rows_out)
matched.to_csv(os.path.join(path_results, 'SEPTUM_LV_MATCHED.tsv'), sep='\t', index=False)
print('\n3c LV lean with LV matched to RV (lean = d(RV) - d(LV), positive = LV; septal blocks resampled)')
print(matched.round(3).to_string(index=False))

fig, axs = plt.subplots(1, 2, figsize=(11, 3.4), sharey=True)
set_cols = dict(zip(snv_sets, ['#252525', '#e08214', '#9e9ac8', '#31a354']))
for ax, scheme in zip(axs, ['blocks', 'depth']):
    for gi, gname in enumerate(sept_groups):
        for si, set_name in enumerate(snv_sets):
            v, raw = draws_out[(set_name, gname, scheme)]
            x = gi + (si - 1.5) * .18
            ax.boxplot(v, positions=[x], widths=.15, showfliers=False, patch_artist=True,
                       boxprops=dict(facecolor=set_cols[set_name], alpha=.45, linewidth=.5),
                       medianprops=dict(color='k', linewidth=.8), whiskerprops=dict(linewidth=.5),
                       capprops=dict(linewidth=.5))
            ax.scatter(x, raw, marker='_', s=40, color='#c0392b', zorder=5, linewidth=1.2)
    ax.axhline(0, color='k', linestyle='--', linewidth=.6)
    ax.set_xticks(range(len(sept_groups)))
    ax.set_xticklabels(list(sept_groups))
    plu.format_ax(ax=ax, ylabel='Lean, d(RV) - d(LV)' if scheme == 'blocks' else None,
                  title={'blocks': 'LV reduced to 2 blocks (6 pairs)',
                         'depth': 'LV: 8 cuts, thinned to RV depth'}[scheme], reduced_spines=True)
plu.add_legend(colors=set_cols, label='SNV set (red tick: raw lean, full LV)', ax=axs[1], ticks_size=6,
               artists_size=5, label_size=6, loc='upper right', bbox_to_anchor=(1, 1), ncols=1)
fig.subplots_adjust(left=.07, right=.98, top=.88, bottom=.12, wspace=.06)
fig.savefig(os.path.join(path_figures, 'septum_lv_matched.pdf'))



##


# 3d. Distribution of the septal cuts' distances to pooled LV and pooled RV, compared as paired
# values (same cut, two distances): Wilcoxon signed-rank on d(RV) - d(LV) per cut, and the same on
# block means (cuts of a block are not independent; blocks are the sampling units). Also with LV
# matched to RV (8 random LV cuts, counts thinned to RV depth), over N_MATCH draws.
N_MATCH = 500


def paired_summary(d_lv, d_rv, blocks_):
    diff = d_rv - d_lv
    bl = pd.Series(diff).groupby(blocks_).mean()
    return dict(n_cuts=len(diff), median_d_LV=np.median(d_lv), median_d_RV=np.median(d_rv),
                median_diff=np.median(diff), cuts_closer_LV=int((diff > 0).sum()),
                p_cuts=wilcoxon(diff)[1] if np.any(diff != 0) else np.nan,
                n_blocks=len(bl), blocks_closer_LV=int((bl > 0).sum()),
                p_blocks=wilcoxon(bl)[1] if len(bl) >= 5 else np.nan)


grp_cuts = {'Septum': np.arange(len(septal)), 'LS': np.where(sreg == 'LS')[0],
            'CS': np.where(sreg == 'CS')[0], 'RS': np.where(sreg == 'RS')[0],
            'LS + RS': np.where(np.isin(sreg, ['LS', 'RS']))[0]}
rows3d = []
for gname, k in grp_cuts.items():
    rows3d.append(dict(septal_group=gname, lv_profile='full LV',
                       **paired_summary(d_region[k, 0], d_region[k, 1], schunk[k])))
rv_prof_all = region_P[1]
match_stats = {g: [] for g in grp_cuts}
for _ in range(N_MATCH):
    r = rng.choice(lv_rows, len(rv_rows), replace=False)
    a_, d_ = A[r].sum(0), D[r].sum(0)
    d2 = np.minimum(dp_rv_all, d_).astype(int)
    lv_prof = rng.hypergeometric(a_.astype(int), (d_ - a_).astype(int), d2) / np.maximum(d2, 1)
    dd = pairwise_distances(AF[septal], np.vstack([lv_prof, rv_prof_all]), metric='cosine')
    for gname, k in grp_cuts.items():
        match_stats[gname].append(paired_summary(dd[k, 0], dd[k, 1], schunk[k]))
for gname in grp_cuts:
    st = pd.DataFrame(match_stats[gname])
    rows3d.append(dict(septal_group=gname, lv_profile=f'matched LV ({N_MATCH} draws, median)',
                       **st.median().to_dict(), frac_draws_p_cuts_lt_05=(st.p_cuts < .05).mean(),
                       frac_draws_p_blocks_lt_05=(st.p_blocks < .05).mean()))
paired = pd.DataFrame(rows3d)
paired.to_csv(os.path.join(path_results, 'SEPTUM_PAIRED_DISTANCES.tsv'), sep='\t', index=False)
print('\n3d Per-cut distance to pooled LV vs pooled RV (paired Wilcoxon; diff = d_RV - d_LV, positive = LV)')
print(paired.round(4).to_string(index=False))

fig, axs = plt.subplots(1, 3, figsize=(10, 3.2), gridspec_kw=dict(width_ratios=[1.4, 1, 1], wspace=.35))
ax = axs[0]
for i, gname in enumerate(['LS', 'CS', 'RS']):
    k = grp_cuts[gname]
    for j, (lab, col) in enumerate([(0, COLORS['LV']), (1, COLORS['RV'])]):
        v = d_region[k, lab]
        ax.boxplot(v, positions=[i + (j - .5) * .35], widths=.28, showfliers=False, patch_artist=True,
                   boxprops=dict(facecolor=col, alpha=.5, linewidth=.6), medianprops=dict(color='k'))
        ax.scatter(i + (j - .5) * .35 + rng.uniform(-.08, .08, len(v)), v, s=6, color=col, edgecolor='none')
    row = paired[(paired.septal_group == gname) & (paired.lv_profile == 'full LV')].iloc[0]
    ax.text(i, ax.get_ylim()[1] if False else .74, f'p={row.p_cuts:.3f}', ha='center', fontsize=6)
ax.set_xticks(range(3)); ax.set_xticklabels(['LS', 'CS', 'RS'])
ax.set_ylim(.25, .78)
plu.format_ax(ax=ax, ylabel='Cosine distance of cut to pooled profile', title='Per cut: LV (red) vs RV (purple)',
              reduced_spines=True)
ax = axs[1]
for gname, col in [('LS', COLORS['LS']), ('RS', COLORS['RS'])]:
    dfb = pd.DataFrame({'b': schunk[grp_cuts[gname]], 'diff': d_region[grp_cuts[gname], 1] - d_region[grp_cuts[gname], 0]}).groupby('b')['diff'].mean()
    ax.scatter(np.full(len(dfb), 0 if gname == 'LS' else 1) + rng.uniform(-.1, .1, len(dfb)), dfb.values, s=20,
               color=col, edgecolor='k', linewidth=.3)
dfb = pd.DataFrame({'b': schunk[grp_cuts['CS']], 'diff': d_region[grp_cuts['CS'], 1] - d_region[grp_cuts['CS'], 0]}).groupby('b')['diff'].mean()
ax.scatter([2], dfb.values, s=20, color=COLORS['CS'], edgecolor='k', linewidth=.3)
ax.axhline(0, color='k', linestyle='--', linewidth=.6)
ax.set_xticks(range(3)); ax.set_xticklabels(['LS', 'RS', 'CS'])
plu.format_ax(ax=ax, ylabel='Block mean d(RV) - d(LV)', title='Per block', reduced_spines=True)
ax = axs[2]
st = pd.DataFrame(match_stats['Septum'])
ax.hist(st.median_diff, bins=30, color='#bdbdbd', edgecolor='white', linewidth=.2)
ax.axvline(0, color='k', linestyle='--', linewidth=.6)
ax.axvline(paired.query('septal_group == "Septum" and lv_profile == "full LV"').median_diff.iloc[0], color='#c0392b')
plu.format_ax(ax=ax, xlabel='Median d(RV) - d(LV), septum cuts', ylabel='Matched-LV draws',
              title='LV matched to RV\n(red: full LV)', reduced_spines=True)
fig.subplots_adjust(left=.07, right=.98, top=.82, bottom=.17)
fig.savefig(os.path.join(path_figures, 'septum_paired_distances.pdf'))



##


# 3e. The paired per-cut test of 3d, repeated per lineage class. Same statistics, on the SNVs of
# one class at a time (cosine distance of each septal cut to the pooled LV and RV profiles over
# those SNVs); LV matched to RV as in 3d.
class_sets = {'All SNVs': all_snvs, **{k: np.where(klass == k)[0]
                                        for k in ['Pre-gastrulation', 'Heart-specific', 'Other shared']}}
rows3e, diff_cut, matched_med = [], {}, {}
for set_name, cols in class_sets.items():
    P_ = np.vstack([pooled(lv_rows, cols), pooled(rv_rows, cols)])
    dd = pairwise_distances(AF[septal][:, cols], P_, metric='cosine')
    for gname, k in grp_cuts.items():
        diff_cut[(set_name, gname)] = dd[k, 1] - dd[k, 0]
        rows3e.append(dict(snv_set=set_name, n_snvs=len(cols), septal_group=gname, lv_profile='full LV',
                           **paired_summary(dd[k, 0], dd[k, 1], schunk[k])))
    stats_m = {g: [] for g in grp_cuts}
    for _ in range(N_MATCH):
        r = rng.choice(lv_rows, len(rv_rows), replace=False)
        a_, d_ = A[r].sum(0), D[r].sum(0)
        d2 = np.minimum(dp_rv_all, d_).astype(int)
        lv_prof = (rng.hypergeometric(a_.astype(int), (d_ - a_).astype(int), d2) / np.maximum(d2, 1))[cols]
        dm = pairwise_distances(AF[septal][:, cols], np.vstack([lv_prof, P_[1]]), metric='cosine')
        for gname, k in grp_cuts.items():
            stats_m[gname].append(paired_summary(dm[k, 0], dm[k, 1], schunk[k]))
    for gname in grp_cuts:
        st = pd.DataFrame(stats_m[gname])
        matched_med[(set_name, gname)] = st.median_diff.median()
        rows3e.append(dict(snv_set=set_name, n_snvs=len(cols), septal_group=gname,
                           lv_profile=f'matched LV ({N_MATCH} draws, median)', **st.median().to_dict(),
                           frac_draws_p_cuts_lt_05=(st.p_cuts < .05).mean(),
                           frac_draws_p_blocks_lt_05=(st.p_blocks < .05).mean()))
by_class = pd.DataFrame(rows3e)
by_class.to_csv(os.path.join(path_results, 'SEPTUM_PAIRED_BY_CLASS.tsv'), sep='\t', index=False)
print('\n3e Paired per-cut test by lineage class (diff = d_RV - d_LV, positive = LV)')
print(by_class.query('septal_group in ["Septum", "LS", "RS", "CS"]')[
    ['snv_set', 'n_snvs', 'septal_group', 'lv_profile', 'median_d_LV', 'median_d_RV', 'median_diff',
     'cuts_closer_LV', 'p_cuts', 'blocks_closer_LV', 'p_blocks']].round(4).to_string(index=False))

fig, axs = plt.subplots(1, 4, figsize=(12, 3.2), sharey=True)
for ax, (set_name, cols) in zip(axs, class_sets.items()):
    for i, gname in enumerate(['LS', 'CS', 'RS', 'Septum']):
        v = diff_cut[(set_name, gname)]
        col = COLORS.get(gname, '#252525')
        ax.boxplot(v, positions=[i], widths=.5, showfliers=False, patch_artist=True,
                   boxprops=dict(facecolor=col, alpha=.4, linewidth=.6), medianprops=dict(color='k'))
        ax.scatter(i + rng.uniform(-.15, .15, len(v)), v, s=6, color=col, edgecolor='none', zorder=3)
        ax.scatter(i, matched_med[(set_name, gname)], marker='D', s=22, color='white', edgecolor='#c0392b',
                   linewidth=1, zorder=5)
        row = by_class[(by_class.snv_set == set_name) & (by_class.septal_group == gname) &
                       (by_class.lv_profile == 'full LV')].iloc[0]
        ax.text(i, .56, f'p={row.p_cuts:.3f}\n{int(row.cuts_closer_LV)}/{int(row.n_cuts)} cuts', ha='center', fontsize=5.5)
    ax.axhline(0, color='k', linestyle='--', linewidth=.6)
    ax.set_xticks(range(4)); ax.set_xticklabels(['LS', 'CS', 'RS', 'Septum'])
    ax.set_ylim(-.45, .68)
    plu.format_ax(ax=ax, ylabel='Per-cut d(RV) - d(LV)' if ax is axs[0] else None,
                  title=f'{set_name} (n={len(cols)})', reduced_spines=True)
fig.text(.5, .01, 'positive = closer to LV; boxes: full LV; white diamonds: median with LV matched to RV size and depth; '
         'p: paired Wilcoxon over cuts (not block-aware)', ha='center', fontsize=6, color='#555555')
fig.subplots_adjust(left=.06, right=.98, top=.88, bottom=.15, wspace=.08)
fig.savefig(os.path.join(path_figures, 'septum_paired_by_class.pdf'))



##


# 3f. Aggregated shared class, and distance rescaling. Raw cosine on VAF is dominated by the
# SNVs with the highest VAF. Four profile scalings are compared, all fixed in advance and all
# reported; none is tuned to the result:
#   raw:         VAF, cosine distance (as in 3d / 3e)
#   relative:    VAF divided by the SNV's pooled heart VAF (every SNV on the same scale), cosine
#   sqrt:        sqrt(VAF), cosine (variance-stabilising for counts)
#   rel. + corr: relative VAF, correlation distance (also removes each profile's mean level)
# Applied to cut profiles and pooled LV / RV profiles alike; LV matched to RV as in 3d.
heart_vaf = np.clip(A.sum(0) / D.sum(0), LL_FLOOR, None)
TRANSFORMS = {
    'raw': (lambda x: x, 'cosine'),
    'relative': (lambda x: x / heart_vaf, 'cosine'),
    'sqrt': (lambda x: np.sqrt(x), 'cosine'),
    'relative + correlation': (lambda x: x / heart_vaf, 'correlation'),
}
agg_sets = {'All SNVs': all_snvs,
            'Pre-gastrulation + Other shared': np.where(klass != 'Heart-specific')[0],
            'Pre-gastrulation': np.where(klass == 'Pre-gastrulation')[0],
            'Heart-specific': np.where(klass == 'Heart-specific')[0]}
agg_groups = {g: grp_cuts[g] for g in ['Septum', 'LS', 'RS', 'CS']}
rv_vaf = pooled(rv_rows, all_snvs)
lv_vaf = pooled(lv_rows, all_snvs)


def paired_dist(lv_prof, cols, f, metric):
    X = f(AF[septal])[:, cols]
    P_ = f(np.vstack([lv_prof, rv_vaf]))[:, cols]
    d_ = pairwise_distances(X, P_, metric=metric)
    return np.nan_to_num(d_, nan=1.0)


rows3f = []
draws = []
for _ in range(N_MATCH):
    r = rng.choice(lv_rows, len(rv_rows), replace=False)
    a_, d_ = A[r].sum(0), D[r].sum(0)
    d2 = np.minimum(dp_rv_all, d_).astype(int)
    draws.append(rng.hypergeometric(a_.astype(int), (d_ - a_).astype(int), d2) / np.maximum(d2, 1))
for set_name, cols in agg_sets.items():
    for t_name, (f, metric) in TRANSFORMS.items():
        full = paired_dist(lv_vaf, cols, f, metric)
        mm = {g: [] for g in agg_groups}
        for lv_prof in draws:
            dm = paired_dist(lv_prof, cols, f, metric)
            for g, k in agg_groups.items():
                mm[g].append(paired_summary(dm[k, 0], dm[k, 1], schunk[k]))
        for g, k in agg_groups.items():
            rows3f.append(dict(snv_set=set_name, n_snvs=len(cols), transform=t_name, septal_group=g,
                               lv_profile='full LV', **paired_summary(full[k, 0], full[k, 1], schunk[k])))
            st = pd.DataFrame(mm[g])
            rows3f.append(dict(snv_set=set_name, n_snvs=len(cols), transform=t_name, septal_group=g,
                               lv_profile='matched LV (median)', **st.median().to_dict(),
                               frac_draws_p_cuts_lt_05=(st.p_cuts < .05).mean()))
rescaled = pd.DataFrame(rows3f)
rescaled.to_csv(os.path.join(path_results, 'SEPTUM_PAIRED_RESCALED.tsv'), sep='\t', index=False)
print('\n3f Aggregated shared class and distance rescaling (septum cuts, positive = closer to LV)')
print(rescaled.query('septal_group == "Septum"')[
    ['snv_set', 'transform', 'lv_profile', 'median_diff', 'cuts_closer_LV', 'p_cuts', 'blocks_closer_LV',
     'p_blocks']].round(4).to_string(index=False))
print('  LS / RS, matched LV:')
print(rescaled.query('septal_group in ["LS", "RS"] and lv_profile != "full LV"')[
    ['snv_set', 'transform', 'septal_group', 'median_diff', 'cuts_closer_LV', 'p_cuts']].round(4).to_string(index=False))

fig, axs = plt.subplots(1, 2, figsize=(9, 3), sharey=True)
for ax, prof in zip(axs, ['full LV', 'matched LV (median)']):
    sub = rescaled.query('septal_group == "Septum" and lv_profile == @prof')
    frac = sub.pivot(index='snv_set', columns='transform', values='cuts_closer_LV').reindex(
        index=list(agg_sets), columns=list(TRANSFORMS)) / 36
    pv = sub.pivot(index='snv_set', columns='transform', values='p_cuts').reindex(
        index=list(agg_sets), columns=list(TRANSFORMS))
    im = ax.imshow(frac.values, cmap='RdBu_r', vmin=.2, vmax=.8, aspect='auto')
    for yi in range(frac.shape[0]):
        for xi in range(frac.shape[1]):
            ax.text(xi, yi, f'{frac.values[yi, xi]:.2f}\np={pv.values[yi, xi]:.3f}', ha='center', va='center', fontsize=6)
    ax.set_xticks(range(len(TRANSFORMS))); ax.set_xticklabels(list(TRANSFORMS), rotation=25, ha='right', fontsize=6.5)
    ax.set_yticks(range(len(agg_sets)))
    ax.set_yticklabels([f'{k} ({len(v)})' for k, v in agg_sets.items()], fontsize=6.5)
    ax.set_title('Full LV' if prof == 'full LV' else 'LV matched to RV (median of draws)', fontsize=8)
cb = fig.colorbar(im, ax=axs, pad=.02, fraction=.03)
cb.set_label('Share of septal cuts closer to LV', fontsize=6)
cb.ax.tick_params(labelsize=5)
fig.subplots_adjust(left=.22, right=.88, top=.88, bottom=.27, wspace=.08)
fig.savefig(os.path.join(path_figures, 'septum_paired_rescaled.pdf'))



##


# 3g. Final cut-level figures. Raw VAF, cosine distance, full LV profile; statistic = septal cuts
# closer to LV than to RV. p: two-sided paired Wilcoxon on d(RV) - d(LV) over cuts, unadjusted;
# cuts are treated as independent. The shared set was defined after seeing the data.
FIG_SETS = {'All SNVs': all_snvs,
            'Pre-gastrulation + Other shared': agg_sets['Pre-gastrulation + Other shared'],
            'Heart-specific': agg_sets['Heart-specific']}
SHOW_GROUPS = ['Septum', 'LS', 'CS', 'RS']
identity = lambda x: x
cut_tab, test_tab, fig_data = [], [], {}
for set_name, cols in FIG_SETS.items():
    full = paired_dist(lv_vaf, cols, identity, 'cosine')
    fig_data[set_name] = dict(d_lv=full[:, 0], d_rv=full[:, 1], closer=full[:, 0] < full[:, 1])
    cut_tab.append(pd.DataFrame({'snv_set': set_name, 'n_snvs': len(cols), 'Sample_ID': cuts[septal],
                                 'region': sreg, 'd_LV': full[:, 0], 'd_RV': full[:, 1],
                                 'closer_LV': full[:, 0] < full[:, 1]}))
    for g in SHOW_GROUPS:
        k = grp_cuts[g]
        diff = full[k, 1] - full[k, 0]
        test_tab.append(dict(snv_set=set_name, septal_group=g, n_cuts=len(k), cuts_closer_LV=int((diff > 0).sum()),
                             p_wilcoxon=wilcoxon(diff)[1]))
pd.concat(cut_tab).to_csv(os.path.join(path_results, 'SEPTUM_LV_CUTS.tsv'), sep='\t', index=False)
tests = pd.DataFrame(test_tab)
tests.to_csv(os.path.join(path_results, 'SEPTUM_LV_CUTS_TESTS.tsv'), sep='\t', index=False)

print('\n3g Septal cuts closer to LV (raw cosine, full LV)')
print(tests.round(4).to_string(index=False))
expected = {'All SNVs': 25, 'Pre-gastrulation + Other shared': 26, 'Heart-specific': 24}
for set_name, d in fig_data.items():
    assert int(d['closer'].sum()) == expected[set_name], 'does not reproduce the 3e counts'

# Physical space: each septal sample's mean Euclidean distance (coordinate units) to the LV samples
# and to the RV samples; same paired test. Also how well the spatial lean tracks the molecular lean.
SPACE = 'Physical space'
xyz_s = coords.values[septal]
d_lv_sp = pairwise_distances(xyz_s, coords.values[lv_rows]).mean(1)
d_rv_sp = pairwise_distances(xyz_s, coords.values[rv_rows]).mean(1)
fig_data[SPACE] = dict(d_lv=d_lv_sp, d_rv=d_rv_sp, closer=d_lv_sp < d_rv_sp)
cut_tab.append(pd.DataFrame({'snv_set': SPACE, 'n_snvs': np.nan, 'Sample_ID': cuts[septal], 'region': sreg,
                             'd_LV': d_lv_sp, 'd_RV': d_rv_sp, 'closer_LV': d_lv_sp < d_rv_sp}))
for g in SHOW_GROUPS:
    k = grp_cuts[g]
    diff = d_rv_sp[k] - d_lv_sp[k]
    test_tab.append(dict(snv_set=SPACE, septal_group=g, n_cuts=len(k), cuts_closer_LV=int((diff > 0).sum()),
                         p_wilcoxon=wilcoxon(diff)[1]))
pd.concat(cut_tab).to_csv(os.path.join(path_results, 'SEPTUM_LV_CUTS.tsv'), sep='\t', index=False)
tests = pd.DataFrame(test_tab)
tests.to_csv(os.path.join(path_results, 'SEPTUM_LV_CUTS_TESTS.tsv'), sep='\t', index=False)
print('\n  physical space (mean Euclidean distance to LV / RV samples)')
print(tests.query('snv_set == @SPACE').round(4).to_string(index=False))
lean_space = d_rv_sp - d_lv_sp
for set_name in FIG_SETS:
    d = fig_data[set_name]
    rho, p_rho = spearmanr(d['d_rv'] - d['d_lv'], lean_space)
    print(f'  molecular lean ({set_name}) vs spatial lean across the 36 samples: Spearman {rho:+.2f} (p = {p_rho:.3f})')

# Scatter: distance to RV against distance to LV, one panel per SNV set plus physical space
fig, axs = plt.subplots(1, len(fig_data), figsize=(12.5, 3.5))
lim_mol = (np.floor(min(min(d['d_rv'].min(), d['d_lv'].min()) for k_, d in fig_data.items() if k_ != SPACE) * 20) / 20,
           np.ceil(max(max(d['d_rv'].max(), d['d_lv'].max()) for k_, d in fig_data.items() if k_ != SPACE) * 20) / 20)
lim_sp = (np.floor(min(d_rv_sp.min(), d_lv_sp.min()) / 100) * 100, np.ceil(max(d_rv_sp.max(), d_lv_sp.max()) / 100) * 100)
for c, (set_name, d) in enumerate(fig_data.items()):
    ax = axs[c]
    lim = lim_sp if set_name == SPACE else lim_mol
    ax.plot(lim, lim, color='k', linewidth=.6, linestyle='--', zorder=1)
    ax.scatter(d['d_rv'], d['d_lv'], s=22, c=[COLORS[r] for r in sreg], edgecolor='k', linewidth=.3, zorder=3)
    for i, g in enumerate(SHOW_GROUPS):
        row = tests.query('snv_set == @set_name and septal_group == @g').iloc[0]
        ptxt = 'p < 0.001' if row.p_wilcoxon < .001 else f'p = {row.p_wilcoxon:.3f}'
        ax.text(.04, .96 - .055 * i, f'{g}: {int(row.cuts_closer_LV)}/{int(row.n_cuts)} closer to LV, {ptxt}',
                transform=ax.transAxes, fontsize=6, va='top', color=COLORS.get(g, 'k'))
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect('equal')
    unit = 'Euclidean distance' if set_name == SPACE else 'Distance'
    plu.format_ax(ax=ax, xlabel=f'{unit} to RV', ylabel=f'{unit} to LV' if c in (0, len(fig_data) - 1) else None,
                  title=(f'{set_name} (mean distance to samples)' if set_name == SPACE
                         else f'{set_name} (n={len(FIG_SETS[set_name])})'), reduced_spines=True)
plu.add_legend(colors={r: COLORS[r] for r in ['LS', 'CS', 'RS']}, label='Septal section', ax=axs[0],
               ticks_size=6, artists_size=5, label_size=6, loc='lower right', bbox_to_anchor=(1, .02), ncols=1)
fig.text(.5, .01, 'SNV panels: cosine distance on raw VAF. Physical panel: mean Euclidean distance (coordinate units) to the LV / RV '
         'samples. Below the diagonal = closer to LV. p: paired Wilcoxon, two-sided, unadjusted; samples treated as independent. '
         'The shared set was defined after seeing the data.', ha='center', fontsize=5.5, color='#555555')
fig.subplots_adjust(left=.06, right=.99, top=.9, bottom=.14, wspace=.22)
fig.savefig(os.path.join(path_figures, 'septum_lv_cuts.pdf'))

# Stacked bars: cuts closer to LV / to RV per septal section
fig, axs = plt.subplots(1, len(FIG_SETS), figsize=(9.5, 3), sharey=True)
for c, set_name in enumerate(FIG_SETS):
    d = fig_data[set_name]
    ax = axs[c]
    for i, g in enumerate(SHOW_GROUPS):
        k = grp_cuts[g]
        n_lv = int(d['closer'][k].sum())
        n_rv = len(k) - n_lv
        for bottom, n_, col, lab in [(0, n_lv, COLORS['LV'], 'Closer to LV'), (n_lv / len(k), n_rv, COLORS['RV'], 'Closer to RV')]:
            ax.bar(i, n_ / len(k), bottom=bottom, color=col, width=.7, edgecolor='white', linewidth=.6,
                   label=lab if (c == 0 and i == 0) else None)
            if n_:
                ax.text(i, bottom + n_ / len(k) / 2, str(n_), ha='center', va='center', fontsize=7, color='white')
    ax.axhline(.5, color='k', linestyle='--', linewidth=.6)
    ax.set_xticks(range(len(SHOW_GROUPS)))
    ax.set_xticklabels([f'{g}\n(n={len(grp_cuts[g])})' for g in SHOW_GROUPS], fontsize=7)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, .25, .5, .75, 1])
    ax.set_yticklabels(['0', '25', '50', '75', '100'])
    plu.format_ax(ax=ax, ylabel='% LCM samples' if c == 0 else None, title=f'{set_name} (n={len(FIG_SETS[set_name])})',
                  reduced_spines=True)
axs[0].legend(frameon=False, fontsize=6, loc='lower left', bbox_to_anchor=(0, 1.12), ncols=2)
fig.text(.5, .01, 'numbers are LCM samples; CS has 4 (descriptive)', ha='center', fontsize=5.5, color='#555555')
fig.subplots_adjust(left=.07, right=.98, top=.8, bottom=.2, wspace=.1)
fig.savefig(os.path.join(path_figures, 'septum_lv_cuts_share.pdf'))


# Physical space on its own: scatter with the paired tests, and the stacked share of samples closer to LV / RV
d = fig_data[SPACE]
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.4), gridspec_kw=dict(width_ratios=[1.15, 1]))
ax = axs[0]
ax.plot(lim_sp, lim_sp, color='k', linewidth=.6, linestyle='--', zorder=1)
ax.scatter(d['d_rv'], d['d_lv'], s=22, c=[COLORS[r] for r in sreg], edgecolor='k', linewidth=.3, zorder=3)
for i, g in enumerate(SHOW_GROUPS):
    row = tests.query('snv_set == @SPACE and septal_group == @g').iloc[0]
    ptxt = 'p < 0.001' if row.p_wilcoxon < .001 else f'p = {row.p_wilcoxon:.3f}'
    ax.text(.04, .96 - .06 * i, f'{g}: {int(row.cuts_closer_LV)}/{int(row.n_cuts)} closer to LV, {ptxt}',
            transform=ax.transAxes, fontsize=6.5, va='top', color=COLORS.get(g, 'k'))
ax.set_xlim(lim_sp); ax.set_ylim(lim_sp); ax.set_aspect('equal')
plu.format_ax(ax=ax, xlabel='Mean Euclidean distance to RV samples', ylabel='Mean Euclidean distance to LV samples',
              title='Physical space', reduced_spines=True)
plu.add_legend(colors={r: COLORS[r] for r in ['LS', 'CS', 'RS']}, label='Septal section', ax=ax,
               ticks_size=6, artists_size=5, label_size=6, loc='lower right', bbox_to_anchor=(1, .02), ncols=1)
ax = axs[1]
for i, g in enumerate(SHOW_GROUPS):
    k = grp_cuts[g]
    n_lv = int(d['closer'][k].sum())
    n_rv = len(k) - n_lv
    for bottom, n_, col, lab in [(0, n_lv, COLORS['LV'], 'Closer to LV'), (n_lv / len(k), n_rv, COLORS['RV'], 'Closer to RV')]:
        ax.bar(i, n_ / len(k), bottom=bottom, color=col, width=.7, edgecolor='white', linewidth=.6,
               label=lab if i == 0 else None)
        if n_:
            ax.text(i, bottom + n_ / len(k) / 2, str(n_), ha='center', va='center', fontsize=7, color='white')
ax.axhline(.5, color='k', linestyle='--', linewidth=.6)
ax.set_xticks(range(len(SHOW_GROUPS)))
ax.set_xticklabels([f'{g}\n(n={len(grp_cuts[g])})' for g in SHOW_GROUPS], fontsize=7)
ax.set_ylim(0, 1)
ax.set_yticks([0, .25, .5, .75, 1])
ax.set_yticklabels(['0', '25', '50', '75', '100'])
ax.legend(frameon=False, fontsize=6, loc='lower left', bbox_to_anchor=(0, 1.0), ncols=2)
plu.format_ax(ax=ax, ylabel='% LCM samples', reduced_spines=True)
fig.text(.5, .01, 'distance of each septal sample to a ventricle = mean Euclidean distance (coordinate units) to its samples; '
         'p: paired Wilcoxon, two-sided, unadjusted; CS has 4 samples (descriptive)', ha='center', fontsize=5.5, color='#555555')
fig.subplots_adjust(left=.1, right=.98, top=.86, bottom=.2, wspace=.3)
fig.savefig(os.path.join(path_figures, 'septum_lv_space.pdf'))

print('\nwritten: septum_lv_cuts, septum_lv_cuts_share, septum_lv_space, septum_paired_rescaled, septum_paired_by_class, septum_paired_distances, septum_lv_matched, septum_cut_lean, septum_ventricle_topology, markers_exclusive_3d, chunks_3d, septum_matching')
