"""
Figures for the lineage analysis, following FIGURES.md.

Analysis 1, lineage characterisation
  embryonic_lineage_features.pdf   size, total depth and pooled cell fraction of the three classes
  genotype_call_composition.pdf    present / absent / undetermined per class
  embryonic_lineage_carriers.pdf   carrier samples and cell fraction in them
  lineage_enrichment.pdf           region VAF relative to heart over region-vs-rest test
  lineage_enrichment_af.pdf        raw pooled region VAF, carrier-cut share and the same test
  enrichment_type.pdf              extended patches vs scattered
  example_enriched.pdf             3D examples of the two patterns
  example_spread.pdf / example_local.pdf   3D read counts: enriched and spread / enriched but local
  lineage_sweep.pdf                cell fraction of every enriched lineage

Analysis 2, phylogenetic relationships
  phylogenetic_relationships.pdf   region dendrogram beside the region VAF matrix
  phylogenetic_relationships_enriched.pdf   the same on enriched SNVs only
  phylogenetic_relationships_markers.pdf    the same on the SNVs that mark one clade each
  markers_exclusive.pdf            markers plus region-exclusive SNVs: region tree and heatmap, per-sample heatmap
  snv_read_support.pdf             density of mean AD, mean DP and samples with >=1 alt read per SNV
  block_hierarchy.pdf              the same construction on the 14 tissue blocks
  lcm_heatmap.pdf                  the same matrix at single-cut resolution
  distances.pdf                    region distances behind the dendrogram
  hierarchy_support.pdf            clade support across profile constructions and lineage classes

Analysis 3
  septum_ventricle_asymmetry.pdf   ventricular lean, its test, leave-one-chunk-out
"""

import os
import itertools
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
from scipy.spatial.distance import squareform
from sklearn.metrics import pairwise_distances
from scipy.stats import binom, mannwhitneyu
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
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
STATE_COLORS = {'Present': '#252525', 'Absent': '#969696', 'Undetermined': '#e4e4e4'}
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

klass = summary['lineage_class'].reindex(muts).values
class_order = pd.Series(klass).value_counts().index.tolist()      # largest class first

# Region x mutation tables, and the enrichment ordering used by several figures. Enrichment
# is the read-level test of each region against the rest of the heart (one-sided binomial,
# rest-of-heart VAF floored at the site background); the hit rule (BH q < 0.1 and alt reads in
# >=30% of the region's cuts) and its cut-label permutation are computed in 8.lineage_analysis.py
ENRICH_Q = 0.1
MIN_FRAC_CUTS_ALT = 0.3
enrich_perm = pd.read_csv(os.path.join(path_results, 'ENRICHMENT_PERMUTATION.tsv'), sep='\t',
                          header=None, index_col=0)[1]
region_tab = lambda col: enrich.pivot_table(index='region', columns='mutation_id', values=col).reindex(REGIONS)[muts]
region_p = region_tab('p_binom')
region_af = region_tab('AF')
heart_af = ((enrich.groupby('mutation_id')['AD'].sum() / enrich.groupby('mutation_id')['DP'].sum())
            .reindex(muts).values)                         # pooled over all heart cuts
region_rel = (region_af - heart_af[None, :]) / heart_af[None, :]
logp = -np.log10(region_p.clip(lower=1e-10))
max_logp = np.nan_to_num(logp.values).max(0)
region_call = region_tab('enriched').fillna(0).astype(bool)
is_enriched = region_call.any().reindex(muts).fillna(False).values
enrich_order = np.concatenate([
    np.where(~is_enriched)[0][np.argsort(max_logp[~is_enriched])],
    np.where(is_enriched)[0][np.argsort(max_logp[is_enriched])],
])
n_present = present.sum(0)
mean_cf = np.array([2 * AF[present[:, j], j].mean() if present[:, j].any() else np.nan
                    for j in range(len(muts))])

print(f'{len(muts)} mutations x {len(cuts)} cuts | enriched {int(is_enriched.sum())} | '
      f'Pre-gastrulation {(klass == "Pre-gastrulation").sum()}')


##


# 1a. The three lineage classes. Cell fraction is 2 x sum(AD) / sum(DP) pooled over all
# heart cuts; "carrier" samples are the cuts with a present call
pooled_cf = (2 * AD_df.sum() / DP_df.sum()).reindex(muts).values
total_depth = DP_df.sum().reindex(muts).values.astype(float)
LABEL_POOLED = 'Cell fraction, pooled'
LABEL_DEPTH = 'Total depth, all heart\nsamples (log10)'
LABEL_N = 'n carrier samples'
LABEL_MEAN = 'Cell fraction\nin carrier samples'


def counts_panel(ax):
    counts = pd.Series(klass).value_counts().reindex(class_order)
    ax.bar(range(len(class_order)), counts.values, color=[CLASS_COLORS[k] for k in class_order], width=.6)
    for i, v in enumerate(counts.values):
        ax.text(i, v + 1, str(v), ha='center', fontsize=8)
    plu.format_ax(ax=ax, xticks=class_order, ylabel='n SNVs', rotx=20, reduced_spines=True)


def box_dot_panel(ax, values, ylabel):
    """
    Per-SNV values of the classes as box plot and dots. Brackets give all pairwise
    two-sided Mann-Whitney p-values across SNVs, Holm-adjusted over the pairs.
    """
    jitter = np.random.default_rng(1)
    data = [values[klass == k] for k in class_order]
    bp = ax.boxplot(data, positions=range(len(data)), widths=.5, showfliers=False,
                    patch_artist=True, medianprops=dict(color='k', linewidth=1.2),
                    whiskerprops=dict(linewidth=.8), capprops=dict(linewidth=.8),
                    boxprops=dict(linewidth=.8))
    for patch, k in zip(bp['boxes'], class_order):
        patch.set(facecolor=CLASS_COLORS[k], alpha=.35, edgecolor=CLASS_COLORS[k])
    for i, (k, v) in enumerate(zip(class_order, data)):
        ax.scatter(i + jitter.uniform(-.16, .16, len(v)), v, s=9, color=CLASS_COLORS[k],
                   edgecolor='none', alpha=.9, zorder=3)
    pairs = list(itertools.combinations(range(len(data)), 2))
    p_adj = multipletests([mannwhitneyu(data[a], data[b])[1] for a, b in pairs],
                          method='holm')[1]
    lo, hi = ax.get_ylim()
    span = hi - lo
    for (a, b), p_val in zip(pairs, p_adj):
        y = hi + (.04 + .13 * (b - a - 1)) * span
        ax.plot([a, a, b, b], [y - .015 * span, y, y, y - .015 * span], color='k', linewidth=.8)
        ax.text((a + b) / 2, y + .01 * span, 'p < 0.001' if p_val < .001 else f'p = {p_val:.3f}',
                ha='center', va='bottom', fontsize=6)
    ax.set_ylim(lo, hi + .3 * span)
    plu.format_ax(ax=ax, xticks=class_order, ylabel=ylabel, rotx=20, reduced_spines=True)


def class_ratio(values, adjust):
    """
    Heart-specific over pre-gastrulation ratio of a measure (log-linear model, HC3
    errors), optionally at equal total depth: log(value) ~ class [+ log10 depth].
    """
    sel = np.isin(klass, ['Heart-specific', 'Pre-gastrulation'])
    d = pd.DataFrame({'y': np.log(values[sel]), 'H': (klass[sel] == 'Heart-specific').astype(float),
                      'ld': np.log10(total_depth[sel])})
    fit = smf.ols('y ~ H + ld' if adjust else 'y ~ H', d).fit(cov_type='HC3')
    lo, hi = fit.conf_int().loc['H']
    return np.exp(fit.params['H']), np.exp(lo), np.exp(hi), fit.pvalues['H']


def ratio_panel(ax, measures):
    """
    Class ratios per measure: open dot unadjusted, filled dot at equal total depth.
    """
    for i, (values, label) in enumerate(measures):
        for adjust, y, fill in [(False, i + .14, 'white'), (True, i - .14, 'k')]:
            r, lo, hi, p_val = class_ratio(values, adjust)
            ax.plot([lo, hi], [y, y], color='k', linewidth=1, zorder=2)
            ax.scatter(r, y, s=22, facecolor=fill, edgecolor='k', linewidth=1, zorder=3)
            ax.text(2.6, y, 'p < 0.001' if p_val < .001 else f'p = {p_val:.3f}',
                    va='center', fontsize=6)
    ax.axvline(1, color='#c0392b', linestyle='--', linewidth=.8)
    ax.set_xscale('log')
    ax.set_xticks([.8, 1, 1.5, 2])
    ax.set_xticklabels(['0.8', '1', '1.5', '2'])
    ax.minorticks_off()
    ax.set_xlim(.75, 2.5)
    ax.set_ylim(len(measures) - .5, -.5)
    ax.scatter([], [], s=22, facecolor='white', edgecolor='k', label='unadjusted')
    ax.scatter([], [], s=22, facecolor='k', edgecolor='k', label='equal total depth')
    ax.legend(frameon=False, fontsize=6, loc='lower left', bbox_to_anchor=(0, 1.0), ncols=2)
    plu.format_ax(ax=ax, yticks=[m[1] for m in measures],
                  xlabel='Heart-specific / pre-gastrulation', reduced_spines=True)


# Figure 1: size of the classes, coverage, pooled cell fraction and the depth-adjusted ratios
fig = plt.figure(figsize=(13, 3.2))
gs = fig.add_gridspec(1, 4, width_ratios=[.8, 1, 1, 1.15], wspace=.6)
counts_panel(fig.add_subplot(gs[0, 0]))
box_dot_panel(fig.add_subplot(gs[0, 1]), np.log10(total_depth), LABEL_DEPTH)
box_dot_panel(fig.add_subplot(gs[0, 2]), pooled_cf, LABEL_POOLED)
ratio_panel(fig.add_subplot(gs[0, 3]), [(pooled_cf, LABEL_POOLED), (n_present.astype(float), LABEL_N),
                                        (mean_cf, LABEL_MEAN)])
fig.tight_layout()
fig.subplots_adjust(bottom=.28)
fig.savefig(os.path.join(path_figures, 'embryonic_lineage_features.pdf'))

# Figure 3: carrier samples and the cell fraction in them
fig = plt.figure(figsize=(6.5, 3.2))
gs = fig.add_gridspec(1, 2, wspace=.55)
box_dot_panel(fig.add_subplot(gs[0, 0]), n_present.astype(float), LABEL_N)
box_dot_panel(fig.add_subplot(gs[0, 1]), mean_cf, LABEL_MEAN)
fig.tight_layout()
fig.subplots_adjust(bottom=.28)
fig.savefig(os.path.join(path_figures, 'embryonic_lineage_carriers.pdf'))


##


# Figure 2: how the cells of each class were called
cells = pd.DataFrame({
    'Class': np.repeat(klass, len(cuts)),
    'Call': pd.Series(STATE.values.T.ravel()).map(
        {'present': 'Present', 'absent': 'Absent',
         'weak': 'Undetermined', 'undetermined': 'Undetermined'}).values,
})
cells['Class'] = pd.Categorical(cells['Class'], categories=class_order[::-1])   # drawn bottom to top
cells['Call'] = pd.Categorical(cells['Call'], categories=list(STATE_COLORS))

fig, ax = plt.subplots(figsize=(6.2, 2.3))
plu.bb_plot(cells, cov1='Class', cov2='Call', categorical_cmap=STATE_COLORS, ax=ax)
plu.format_ax(ax=ax, xlabel='Fraction of genotype calls')
plu.add_legend(colors=STATE_COLORS, label='Genotype call', ax=ax, ticks_size=6,
               artists_size=5, label_size=7, loc='upper left',
               bbox_to_anchor=(1.01, 1), ncols=1)
fig.tight_layout()
fig.subplots_adjust(right=.78)
fig.savefig(os.path.join(path_figures, 'genotype_call_composition.pdf'))


##


# 1b. Region VAF over region enrichment, mutations ordered by enrichment
cols = enrich_order
split_at = int((~is_enriched).sum()) - .5
vmax_af = np.nanpercentile(region_af.values, 98)
vmax_rel = np.nanpercentile(region_rel.values, 98)
vmax_p = np.nanpercentile(logp.values, 99)
cmap_af = matplotlib.colormaps['afmhot_r'].with_extremes(bad='#dddddd')
cmap_rel = matplotlib.colormaps['RdBu_r'].with_extremes(bad='#dddddd')
cmap_p = matplotlib.colormaps['viridis'].with_extremes(bad='#dddddd')

fig, axs = plt.subplots(2, 1, figsize=(13, 5), sharex=True)

ax = axs[0]
im = ax.imshow(np.ma.masked_invalid(region_rel.values[:, cols]), cmap=cmap_rel, aspect='auto',
               norm=matplotlib.colors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=vmax_rel))
ax.axvline(split_at, color='white', linewidth=3)
plu.format_ax(ax=ax, yticks=REGIONS, xticks=[])
ax.set_title('Region VAF relative to the whole heart', pad=24)
ax.text(1, 1.32, f'{int(enrich_perm.n_enriched_snvs)} enriched SNVs vs '
        f'{enrich_perm.expected_by_chance:.1f} expected by chance '
        f'({int(enrich_perm.n_permutations)} cut-label permutations; global p '
        f'{"< " if enrich_perm.p_global <= 1 / enrich_perm.n_permutations + 1e-9 else "= "}'
        f'{max(enrich_perm.p_global, 1 / enrich_perm.n_permutations):.3f}, '
        f'empirical FDR {enrich_perm.empirical_FDR:.2f})',
        transform=ax.transAxes, ha='right', va='bottom', fontsize=7)
n_non = int((~is_enriched).sum())
for x0, x1, lab in [(-.5, split_at - .3, f'non-enriched (n={n_non})'),
                    (split_at + .3, len(cols) - .5, f'enriched (n={int(is_enriched.sum())})')]:
    ax.plot([x0, x1], [-.9, -.9], color='k', linewidth=.8, clip_on=False)
    ax.text((x0 + x1) / 2, -1.05, lab, ha='center', va='bottom', fontsize=7)
ax.set_ylim(len(REGIONS) - .5, -.5)
cb = fig.colorbar(im, ax=ax, pad=.005, fraction=.015)
cb.set_label('(AF region - AF heart)\n/ AF heart', fontsize=6)
cb.ax.tick_params(labelsize=6)

ax = axs[1]
im = ax.imshow(np.ma.masked_invalid(logp.values[:, cols]), cmap=cmap_p,
               aspect='auto', vmin=0, vmax=vmax_p)
ax.axvline(split_at, color='white', linewidth=3)
for yi in range(len(REGIONS)):
    for xi, j in enumerate(cols):
        if region_call.values[yi, j]:
            ax.text(xi, yi, '*', ha='center', va='center', fontsize=7, color='white')
plu.format_ax(ax=ax, yticks=REGIONS, xticks=[],
              title=f'Region vs rest of heart (-log10 p, one-sided binomial on reads; '
                    f'* hit: BH q < {ENRICH_Q} and alt reads in >={MIN_FRAC_CUTS_ALT:.0%} of region cuts)')
ax.set_xticks(range(len(cols)))
ax.set_xticklabels([m.rsplit('_', 2)[0].replace('_', ':') for m in muts[cols]],
                   rotation=90, fontsize=4.8, ha='center', va='top')
for lab, j in zip(ax.get_xticklabels(), cols):       # enriched SNVs in bold
    lab.set_fontweight('bold' if is_enriched[j] else 'normal')
ax.tick_params(axis='x', length=2, width=.4, pad=1.5)
cb = fig.colorbar(im, ax=ax, pad=.005, fraction=.015)
cb.set_label('-log10 p', fontsize=7)
cb.ax.tick_params(labelsize=6)

fig.subplots_adjust(left=.04, right=.96, top=.88, bottom=.2, hspace=.22)
fig.savefig(os.path.join(path_figures, 'lineage_enrichment.pdf'))


# 1b, alternative. Raw pooled region VAF, the share of region cuts carrying alt reads, and the
# test, all SNVs ordered by pooled heart VAF (highest first); hits marked, enriched SNVs in bold
cols_af = np.argsort(-np.nan_to_num(heart_af))
region_frac = region_tab('frac_cuts_alt')
fig, axs = plt.subplots(3, 1, figsize=(13, 7), sharex=True)
for ax, values, cmap, vmax, label in [
        (axs[0], region_af.values, cmap_af, vmax_af, 'VAF'),
        (axs[1], 100 * region_frac.values, matplotlib.colormaps['Blues'].with_extremes(bad='#dddddd'), 100,
         '% samples'),
        (axs[2], logp.values, cmap_p, vmax_p, '-log10p')]:
    im = ax.imshow(np.ma.masked_invalid(values[:, cols_af]), cmap=cmap, aspect='auto', vmin=0, vmax=vmax)
    for yi in range(len(REGIONS)):
        for xi, j in enumerate(cols_af):
            if region_call.values[yi, j]:
                ax.text(xi, yi, '*', ha='center', va='center', fontsize=7,
                        color='white' if ax is axs[2] or values[yi, j] > vmax / 2 else 'k')
    cb = fig.colorbar(im, ax=ax, pad=.005, fraction=.015)
    cb.set_label(label, fontsize=7)
    cb.ax.tick_params(labelsize=6)
plu.format_ax(ax=axs[0], yticks=REGIONS, xticks=[])
axs[0].set_title('Region VAF (pooled AD / DP)', pad=16)
axs[0].text(1, 1.32, f'{int(enrich_perm.n_enriched_snvs)} enriched SNVs vs '
            f'{enrich_perm.expected_by_chance:.1f} expected by chance '
            f'({int(enrich_perm.n_permutations)} cut-label permutations; global p < '
            f'{max(enrich_perm.p_global, 1 / enrich_perm.n_permutations):.3f}, '
            f'empirical FDR {enrich_perm.empirical_FDR:.2f})',
            transform=axs[0].transAxes, ha='right', va='bottom', fontsize=7)
plu.format_ax(ax=axs[1], yticks=REGIONS, xticks=[],
              title='% samples >0 ALT reads')
plu.format_ax(ax=axs[2], yticks=REGIONS, xticks=[], title='-log10p binomial test')
axs[2].set_xticks(range(len(cols_af)))
axs[2].set_xticklabels([m.rsplit('_', 2)[0].replace('_', ':') for m in muts[cols_af]],
                       rotation=90, fontsize=4.8, ha='center', va='top')
for lab, j in zip(axs[2].get_xticklabels(), cols_af):
    lab.set_fontweight('bold' if is_enriched[j] else 'normal')
axs[2].tick_params(axis='x', length=2, width=.4, pad=1.5)
fig.subplots_adjust(left=.04, right=.96, top=.89, bottom=.14, hspace=.25)
fig.savefig(os.path.join(path_figures, 'lineage_enrichment_af.pdf'))


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
plu.format_ax(ax=ax, xticks=[o.replace(' ', '\n') for o in order], ylabel='n SNVs', rotx=0,
              reduced_spines=True, title=f'k={K_NEIGHBOURS} neighbours')
ax.text(.5, -.32, f'k=4: {alt_counts[4].get("Extended patches", 0)}/'
        f'{alt_counts[4].get("Scattered", 0)}   k=5: '
        f'{alt_counts[5].get("Extended patches", 0)}/{alt_counts[5].get("Scattered", 0)}',
        transform=ax.transAxes, ha='center', fontsize=6, color='#555555')
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'enrichment_type.pdf'))


##


# 1c, companion. One example of each pattern, in space
def pick_example(label):
    """
    Clearest case of a class: the tightest cluster of present cuts among mutations
    with at least 8 of them (extended), the most dispersed among those with at
    least 5 (scattered). Tightness is the mean distance of each present cut to its
    nearest other present cut.
    """
    min_cuts = 8 if label == 'Extended patches' else 5
    nn = {}
    for j in enriched_idx:
        idx = np.where(present[:, j])[0]
        if patch_class[muts[j]] == label and len(idx) >= min_cuts:
            nn[j] = D_space[np.ix_(idx, idx)].min(1).mean()
    if not nn:
        return None
    return (min if label == 'Extended patches' else max)(nn, key=nn.get)


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


# 1c, read-count view. Two enriched lineages spread across their region, and two strongly
# enriched but confined to a few cuts. Cuts with alternate reads filled by their VAF, cuts
# without open. No genotype calls involved.
EXAMPLES_SPREAD_LOCAL = [
    ('Enriched, spread', 'chr5_132654211_G_A', 'RS'),
    ('Enriched, spread', 'chr11_43178783_G_T', 'RS'),
    ('Enriched, local', 'chr11_112067303_A_G', 'LV'),
    ('Enriched, local', 'chrX_131566302_T_G', 'LV'),
]
VAF_MAX_3D = .2


def draw_counts_3d(ax, j):
    carries = AD_df.values[:, j] > 0
    ax.scatter(coords['x'][~carries], coords['y'][~carries], coords['z'][~carries], s=POINT_SIZE,
               facecolor='white', edgecolor='k', linewidth=.3, depthshade=False, zorder=4)
    ax.scatter(coords['x'][carries], coords['y'][carries], coords['z'][carries], s=POINT_SIZE,
               c=AF[carries, j], cmap='afmhot_r', vmin=0, vmax=VAF_MAX_3D,
               edgecolor='k', linewidth=.3, depthshade=False, zorder=5)
    style_3d(ax, coords)


def example_figure(examples, name):
    """
    Sampling layout next to the read-count view of each example, saved as <name>.pdf.
    """
    examples = [(lab, m, r) for lab, m, r in examples if m in set(muts)]
    fig = plt.figure(figsize=(3 * (len(examples) + 1), 3.3))
    ax = fig.add_subplot(1, len(examples) + 1, 1, projection='3d')
    ax.computed_zorder = False
    ax.scatter(coords['x'], coords['y'], coords['z'], s=POINT_SIZE,
               c=[COLORS[r] for r in reg], edgecolor='white', linewidth=.25, depthshade=False, zorder=5)
    style_3d(ax, coords)
    ax.set_title('Sampling', fontsize=8)
    plu.add_legend(colors=COLORS, label='Region', ax=ax, ticks_size=5, artists_size=4,
                   label_size=6, loc='upper center', bbox_to_anchor=(.5, .08), ncols=5)
    etab = enrich.set_index(['mutation_id', 'region'])
    for i, (lab, m, r) in enumerate(examples):
        j = list(muts).index(m)
        row = etab.loc[(m, r)]
        ax = fig.add_subplot(1, len(examples) + 1, i + 2, projection='3d')
        ax.computed_zorder = False
        draw_counts_3d(ax, j)
        ax.set_title(f'{lab} ({r})\n{m.rsplit("_", 2)[0].replace("_", ":")}\n'
                     f'VAF {row.AF:.1%} vs {row.AF_rest:.1%} rest | '
                     f'{row.n_cuts_alt}/{row.n_cuts} cuts, {row.n_blocks_alt}/{row.n_blocks} blocks | '
                     f'-log10p {-np.log10(max(row.p_binom, 1e-30)):.1f}', fontsize=6.5)
    plu.add_cbar(np.array([0, VAF_MAX_3D]), ax=ax, label='VAF', palette='afmhot_r', vmin=0,
                 vmax=VAF_MAX_3D)
    fig.text(.5, .02, 'filled = cut with >=1 alternate read (colour: VAF), open = no alternate read',
             ha='center', fontsize=6, color='#555555')
    fig.subplots_adjust(left=.01, right=.9, top=.8, bottom=.08, wspace=.02)
    fig.savefig(os.path.join(path_figures, f'{name}.pdf'))


example_figure([x for x in EXAMPLES_SPREAD_LOCAL if x[0] == 'Enriched, spread'], 'example_spread')
example_figure([x for x in EXAMPLES_SPREAD_LOCAL if x[0] == 'Enriched, local'], 'example_local')


##


# 1d. How far each enriched lineage swept, cut by cut
sweep = []
for j in enriched_idx:
    for i in np.where(present[:, j])[0]:
        sweep.append(dict(mutation=muts[j].rsplit('_', 2)[0].replace('_', ':'), cell_fraction=2 * AF[i, j],
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
plu.format_ax(ax=ax, xlabel='Enriched lineages', ylabel='Cell fraction', rotx=45,
              xticks_size=6, reduced_spines=True,
              title='Sweep of enriched lineages (black line: mean)')
plt.setp(ax.get_xticklabels(), ha='right', rotation_mode='anchor')
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'lineage_sweep.pdf'))


##


# 2a. Region dendrogram beside the region VAF matrix. Profiles are pooled AD/DP per region
# (the same values as the heatmap); clade support (SNV jackknife / block bootstrap) and the
# consensus call come from 8.lineage_analysis.py (HIERARCHY.tsv)
hierarchy = pd.read_csv(os.path.join(path_results, 'HIERARCHY.tsv'), sep='\t')
D_reg = pairwise_distances(region_af.values, metric='cosine')
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


clade_support = hierarchy.query('profile == "pooled" and snv_set == "All"').set_index('clade')

leaf_order = [REGIONS[i] for i in leaves_list(Z_reg)][::-1]   # top to bottom
# Mutations: region-specific ones (>= SPEC_MIN of their summed regional VAF in one
# region, chance is 1/5) grouped by that region in dendrogram order, sharpest first;
# the broadly shared ones last, by decreasing mean VAF
SPEC_MIN = .35
share = region_af / region_af.sum(0)
dominant = share.idxmax(0).reindex(muts)
specificity = share.max(0).reindex(muts).values
specific = specificity >= SPEC_MIN
groups = []
for r in leaf_order:
    sel = np.where(specific & (dominant.values == r))[0]
    groups.append((r, sel[np.argsort(-specificity[sel])]))
shared = np.where(~specific)[0]
groups.append(('Shared', shared[np.argsort(-np.nan_to_num(region_af.values[:, shared].mean(0)))]))
mut_order = np.concatenate([g for _, g in groups]).astype(int)
block_edges = np.cumsum([len(g) for _, g in groups])[:-1] - .5
group_colors = {**COLORS, 'Shared': '#999999'}


def group_strip(ax, y0, height, grp=None):
    """
    Colour strip over the mutation groups, with their names and sizes.
    """
    start, raised = 0, False
    for name, g in (groups if grp is None else grp):
        ax.add_patch(plt.Rectangle((start - .5, y0), len(g), height, color=group_colors[name],
                                   clip_on=False, linewidth=0))
        narrow = len(g) < 6
        ax.text(start - .5 + len(g) / 2, y0 - .3 * height - (height if narrow and raised else 0),
                f'{name} ({len(g)})', ha='center', va='bottom', fontsize=6)
        raised = narrow and not raised
        start += len(g)



fig = plt.figure(figsize=(10, 3.3))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 3], wspace=.08)

ax = fig.add_subplot(gs[0, 0])
dn = dendrogram(Z_reg, labels=REGIONS, orientation='left', color_threshold=0,
                above_threshold_color='k', ax=ax)
for grp, height in clade_set(Z_reg, REGIONS):
    ys = [dn['ivl'].index(m) for m in grp]
    row = clade_support.loc['+'.join(sorted(grp))]
    ax.text(height - .004, 10 * np.mean(ys) + 5,
            f'{100 * row.support_snv:.0f} / {100 * row.support_block:.0f}',
            fontsize=6, va='center', ha='left', color='#c0392b' if row.consensus else '#999999')
ax.text(0, -.36, 'support %: SNV jackknife / block bootstrap\n(grey: block support < 50%, unresolved)',
        transform=ax.transAxes, fontsize=5.5, color='#555555', va='top')
plu.format_ax(ax=ax, xlabel='Cosine distance', reduced_spines=True)
ax.set_yticks([])
for i, r in enumerate(dn['ivl']):       # leaf i sits at y = 10 * i + 5
    ax.text(-.004, 10 * i + 5, r, ha='left', va='center', fontsize=8, clip_on=False)
ax = fig.add_subplot(gs[0, 1])
heat_rows = [REGIONS.index(r) for r in leaf_order]
ax.imshow(np.ma.masked_invalid(region_af.values[np.ix_(heat_rows, mut_order)]),
          cmap=cmap_af, aspect='auto', vmin=0, vmax=vmax_af)
for e in block_edges:
    ax.axvline(e, color='k', linewidth=.6)
ax.set_yticks([])
group_strip(ax, -1.1, .5)
plu.format_ax(ax=ax, xticks=[], xlabel=f'SNVs (n={len(muts)})')

fig.subplots_adjust(left=.02, right=.95, top=.84, bottom=.33)
pos = ax.get_position()
cax = fig.add_axes((pos.x1 - .12, pos.y0 - .02 - .03, .12, .03))
cb = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0, vmax_af), cmap='afmhot_r'),
                  cax=cax, orientation='horizontal')
cb.set_label('VAF', fontsize=6, labelpad=1)
cb.ax.tick_params(labelsize=5, length=2, pad=1)
fig.savefig(os.path.join(path_figures, 'phylogenetic_relationships.pdf'))


# 2a, enriched SNVs only. Same construction restricted to the enriched SNVs; all SNVs clustered
# together, a colour strip giving the region where each is enriched (its strongest hit)
enr_cols = np.where(is_enriched)[0]
D_enr = pairwise_distances(region_af.values[:, enr_cols], metric='cosine')
Z_enr = linkage(squareform(D_enr, checks=False), method='average')
hit_p = np.where(region_call.values, region_p.values, np.inf)
enr_home = pd.Series([REGIONS[i] for i in hit_p.argmin(0)], index=range(len(muts)))

fig = plt.figure(figsize=(10, 3.7))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 3], wspace=.08)
ax = fig.add_subplot(gs[0, 0])
dn = dendrogram(Z_enr, labels=REGIONS, orientation='left', color_threshold=0,
                above_threshold_color='k', ax=ax)
for grp, height in clade_set(Z_enr, REGIONS):
    ys = [dn['ivl'].index(m) for m in grp]
    ax.text(height - .004, 10 * np.mean(ys) + 5, f'{height:.3f}', fontsize=6, va='center',
            ha='left', color='#555555')
plu.format_ax(ax=ax, xlabel='Cosine distance', reduced_spines=True)
ax.set_yticks([])
for i, r in enumerate(dn['ivl']):
    ax.text(-.004, 10 * i + 5, r, ha='left', va='center', fontsize=8, clip_on=False)
leaf_order_enr = dn['ivl'][::-1]                                   # top to bottom
# All enriched SNVs clustered together on their regional VAF profiles
enr_order = enr_cols[leaves_list(linkage(region_af.values[:, enr_cols].T, method='average',
                                         metric='cosine'))]

ax = fig.add_subplot(gs[0, 1])
im = ax.imshow(np.ma.masked_invalid(region_af.values[np.ix_([REGIONS.index(r) for r in leaf_order_enr],
                                                            enr_order)]),
               cmap=cmap_af, aspect='auto', vmin=0, vmax=vmax_af)
ax.set_yticks([])
for xi, j in enumerate(enr_order):                 # strip: region where each SNV is enriched
    ax.add_patch(plt.Rectangle((xi - .5, -1.1), 1, .5, color=COLORS[enr_home[j]], clip_on=False,
                               linewidth=0))
plu.add_legend(colors={r: COLORS[r] for r in REGIONS}, label='Enriched in', ax=ax, ticks_size=6,
               artists_size=5, label_size=6, loc='lower left', bbox_to_anchor=(0, 1.2), ncols=5)
plu.format_ax(ax=ax, xticks=[], xlabel=f'Enriched SNVs (n={len(enr_cols)})')
fig.subplots_adjust(left=.02, right=.95, top=.74, bottom=.3)
pos = ax.get_position()
cax = fig.add_axes((pos.x1 - .12, pos.y0 - .02 - .03, .12, .03))
cb = fig.colorbar(im, cax=cax, orientation='horizontal')
cb.set_label('VAF', fontsize=6, labelpad=1)
cb.ax.tick_params(labelsize=5, length=2, pad=1)
fig.savefig(os.path.join(path_figures, 'phylogenetic_relationships_enriched.pdf'))


# 2a, clean markers. SNVs whose carrier regions (pooled VAF >= MARKER_FRAC of the SNV's highest
# regional VAF) are exactly one clade of the all-SNV region tree. They are selected to match
# that tree, so the tree built from them illustrates the markers, not independent support.
MARKER_FRAC = .5
tree_clades = [c for c, _h in clade_set(Z_reg, REGIONS)]
carriers = [frozenset(np.array(REGIONS)[region_af.values[:, j] >= MARKER_FRAC * np.nanmax(region_af.values[:, j])])
            for j in range(len(muts))]
marker_groups = []
for c in sorted(tree_clades, key=len):
    sel = np.array([j for j in range(len(muts)) if carriers[j] == c], dtype=int)
    if len(sel):
        marker_groups.append(('+'.join(r for r in REGIONS if r in c), sel))
exclusive_groups = []
for r in REGIONS:
    sel = np.array([j for j in range(len(muts)) if carriers[j] == frozenset({r})], dtype=int)
    if len(sel):
        exclusive_groups.append((r, sel))


def marker_figure(groups, name, title):
    marker_cols = np.concatenate([g for _, g in groups]).astype(int)
    Z_mk = linkage(squareform(pairwise_distances(region_af.values[:, marker_cols], metric='cosine'),
                              checks=False), method='average')

    fig = plt.figure(figsize=(3.2 + .28 * len(marker_cols), 3.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[3.2, .28 * len(marker_cols)], wspace=.08)
    ax = fig.add_subplot(gs[0, 0])
    dn = dendrogram(Z_mk, labels=REGIONS, orientation='left', color_threshold=0,
                    above_threshold_color='k', ax=ax)
    for grp, height in clade_set(Z_mk, REGIONS):
        ys = [dn['ivl'].index(m) for m in grp]
        ax.text(height, 10 * np.mean(ys) + 5 + 1.5, f'{height:.2f}', fontsize=6, va='bottom',
                ha='center', color='#555555', bbox=dict(facecolor='white', edgecolor='none', pad=.3))
    plu.format_ax(ax=ax, xlabel='Cosine distance', reduced_spines=True)
    ax.set_yticks([])
    for i, r in enumerate(dn['ivl']):
        ax.text(-.01, 10 * i + 5, r, ha='left', va='center', fontsize=8, clip_on=False)
    rows_mk = [REGIONS.index(r) for r in dn['ivl'][::-1]]
    ax = fig.add_subplot(gs[0, 1])
    im = ax.imshow(np.ma.masked_invalid(region_af.values[np.ix_(rows_mk, marker_cols)]), cmap=cmap_af,
                   aspect='auto', vmin=0, vmax=vmax_af)
    for e in np.cumsum([len(g) for _, g in groups])[:-1] - .5:
        ax.axvline(e, color='k', linewidth=.8)
    start = 0
    for lab, g in groups:
        ax.text(start - .5 + len(g) / 2, -.75, lab.replace(' ', '\n') if len(g) < 2 else lab,
                ha='center', va='bottom', fontsize=6.5)
        start += len(g)
    ax.set_yticks([])
    ax.set_xticks(range(len(marker_cols)))
    ax.set_xticklabels([muts[j].rsplit('_', 2)[0].replace('_', ':') for j in marker_cols], rotation=45,
                       ha='right', rotation_mode='anchor', fontsize=7.5)
    ax.set_title(f'{title} (n={len(marker_cols)}; carrier regions >= {MARKER_FRAC:.0%} of max VAF)',
                 fontsize=7, pad=22)
    fig.subplots_adjust(left=.02, right=.95, top=.82, bottom=.34)
    pos = ax.get_position()
    cax = fig.add_axes((pos.x1 - .15, .09, .15, .03))
    cb = fig.colorbar(im, cax=cax, orientation='horizontal')
    cb.set_label('VAF', fontsize=6, labelpad=1)
    cb.ax.tick_params(labelsize=5, length=2, pad=1)
    fig.savefig(os.path.join(path_figures, f'{name}.pdf'))


marker_figure(marker_groups, 'phylogenetic_relationships_markers', 'Markers of each clade')


def markers_region_and_samples(groups, name):
    """
    Region dendrogram and region VAF heatmap (top), the same SNVs per LCM sample (bottom, by
    region then block); a strip gives each SNV's lineage class.
    """
    cols = np.concatenate([g for _, g in groups]).astype(int)
    Z_ = linkage(squareform(pairwise_distances(region_af.values[:, cols], metric='cosine'),
                            checks=False), method='average')
    n = len(cols)
    fig = plt.figure(figsize=(2.6 + .2 * n, 8))
    gs = fig.add_gridspec(2, 2, width_ratios=[2.6, .2 * n], height_ratios=[1, 3.2],
                          wspace=.1, hspace=.06)

    ax = fig.add_subplot(gs[0, 0])
    dn = dendrogram(Z_, labels=REGIONS, orientation='left', color_threshold=0,
                    above_threshold_color='k', ax=ax)
    for grp, height in clade_set(Z_, REGIONS):
        ys = [dn['ivl'].index(m) for m in grp]
        ax.text(height, 10 * np.mean(ys) + 8, f'{height:.2f}', fontsize=6, va='bottom', ha='center',
                color='#555555', bbox=dict(facecolor='white', edgecolor='none', pad=.3))
    plu.format_ax(ax=ax, xlabel='Cosine distance', reduced_spines=True)
    ax.xaxis.label.set_fontsize(8)
    ax.tick_params(axis='x', labelsize=6)
    ax.set_yticks([])
    for k, r in enumerate(dn['ivl']):
        ax.text(-.01, 10 * k + 5, r, ha='left', va='center', fontsize=8, clip_on=False)
    reg_order = dn['ivl'][::-1]                                      # top to bottom

    ax = fig.add_subplot(gs[0, 1])
    im = ax.imshow(np.ma.masked_invalid(region_af.values[np.ix_([REGIONS.index(r) for r in reg_order],
                                                                cols)]),
                   cmap=cmap_af, aspect='auto', vmin=0, vmax=vmax_af)
    for e in np.cumsum([len(g) for _, g in groups])[:-1] - .5:
        ax.axvline(e, color='k', linewidth=.8)
    for xi, j in enumerate(cols):                                    # lineage class strip
        ax.add_patch(plt.Rectangle((xi - .5, -1.05), 1, .45, color=CLASS_COLORS[klass[j]],
                                   clip_on=False, linewidth=0))
    start = 0
    for lab, g in groups:
        ax.text(start - .5 + len(g) / 2, -1.2, lab.replace(' ', '\n') if len(g) < 2 else lab,
                ha='center', va='bottom', fontsize=6.5)
        start += len(g)
    ax.set_xticks([])
    ax.set_yticks([])
    plu.add_legend(colors=CLASS_COLORS, label='SNV class', ax=ax, ticks_size=6, artists_size=5,
                   label_size=6, loc='lower left', bbox_to_anchor=(0, 1.32), ncols=3)
    cb = fig.colorbar(im, ax=ax, pad=.01, fraction=.02)
    cb.set_label('VAF (region)', fontsize=6)
    cb.ax.tick_params(labelsize=5)

    ax = fig.add_subplot(gs[1, 1])
    order = [i for r in reg_order for b in sorted(set(chunk)) for i in np.where((reg == r) & (chunk == b))[0]]
    M = AF[np.ix_(order, cols)]
    im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap_af, aspect='auto', vmin=0,
                   vmax=np.nanpercentile(M, 98))
    row_reg, row_blk = np.array(reg[order]), np.array(chunk.values[order] if hasattr(chunk, 'values')
                                                         else chunk[order])
    for e in np.cumsum([len(g) for _, g in groups])[:-1] - .5:
        ax.axvline(e, color='k', linewidth=.8)
    for y in np.where(row_reg[:-1] != row_reg[1:])[0] + .5:
        ax.axhline(y, color='k', linewidth=1)
    for y in np.where(row_blk[:-1] != row_blk[1:])[0] + .5:
        ax.axhline(y, color='#999999', linewidth=.3)
    for k, r in enumerate(row_reg):
        ax.add_patch(plt.Rectangle((-1.4, k - .5), .7, 1, color=COLORS[r], clip_on=False, linewidth=0))
    ax.set_yticks([])
    ax.set_ylabel(f'{len(order)} LCM samples, by region then block', fontsize=8, labelpad=24)
    ax.set_xticks(range(n))
    ax.set_xticklabels([muts[j].rsplit('_', 2)[0].replace('_', ':') for j in cols], rotation=45,
                       ha='right', rotation_mode='anchor', fontsize=6.5)
    cb = fig.colorbar(im, ax=ax, pad=.01, fraction=.02)
    cb.set_label('VAF (sample)', fontsize=6)
    cb.ax.tick_params(labelsize=5)
    ax_leg = fig.add_subplot(gs[1, 0])
    ax_leg.axis('off')
    plu.add_legend(colors=COLORS, label='Region (samples)', ax=ax_leg, ticks_size=7, artists_size=6,
                   label_size=7, loc='center', bbox_to_anchor=(.5, .5), ncols=1)
    fig.subplots_adjust(left=.03, right=.95, top=.86, bottom=.1)
    fig.savefig(os.path.join(path_figures, f'{name}.pdf'))



# 2a, block level. The same construction on the 14 tissue blocks (pooled AD/DP per block).
# Support: SNV jackknife (80%) / cuts resampled with replacement within each block. Region
# assortativity: are same-region blocks closer than different-region ones, tested by
# shuffling region labels across blocks.
blk_names = np.array(sorted(set(chunk)))
blk_reg = pd.Series(reg, index=chunk).groupby(level=0).first().reindex(blk_names).values
blk_cuts = [np.where(chunk == b)[0] for b in blk_names]
AD_v, DP_v = AD_df.values.astype(float), DP_df.values.astype(float)


def block_profiles(cut_sets, cols):
    return np.vstack([AD_v[np.ix_(c, cols)].sum(0) / DP_v[np.ix_(c, cols)].sum(0) for c in cut_sets])


all_snvs = np.arange(len(muts))
P_blk = block_profiles(blk_cuts, all_snvs)
D_blk = pairwise_distances(P_blk, metric='cosine')
Z_blk = linkage(squareform(D_blk, checks=False), method='average')
obs_clades = dict(clade_set(Z_blk, list(blk_names)))
sup_snv = {c: 0 for c in obs_clades}
sup_cut = {c: 0 for c in obs_clades}
for _ in range(N_BOOT):
    cols_b = rng.choice(all_snvs, int(CHAR_FRACTION * len(all_snvs)), replace=False)
    Zb = linkage(squareform(pairwise_distances(block_profiles(blk_cuts, cols_b), metric='cosine'),
                            checks=False), method='average')
    for c, _h in clade_set(Zb, list(blk_names)):
        if c in sup_snv:
            sup_snv[c] += 1
    cuts_b = [rng.choice(c, len(c), replace=True) for c in blk_cuts]
    Zb = linkage(squareform(pairwise_distances(block_profiles(cuts_b, all_snvs), metric='cosine'),
                            checks=False), method='average')
    for c, _h in clade_set(Zb, list(blk_names)):
        if c in sup_cut:
            sup_cut[c] += 1
iu = np.triu_indices(len(blk_names), 1)


def region_gap(labels):
    same = (labels[:, None] == labels[None, :])[iu]
    return D_blk[iu][~same].mean() - D_blk[iu][same].mean()


gap_obs = region_gap(blk_reg)
gap_null = np.array([region_gap(rng.permutation(blk_reg)) for _ in range(10000)])
p_assort = (np.sum(gap_null >= gap_obs) + 1) / (len(gap_null) + 1)

fig = plt.figure(figsize=(10, 4.2))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 3], wspace=.2)
ax = fig.add_subplot(gs[0, 0])
dn = dendrogram(Z_blk, labels=list(blk_names), orientation='left', color_threshold=0,
                above_threshold_color='k', ax=ax)
for grp, height in obs_clades.items():
    if len(grp) > 7:
        continue
    ys = [dn['ivl'].index(m) for m in grp]
    a, b = sup_snv[grp] / N_BOOT, sup_cut[grp] / N_BOOT
    ax.text(height - .004, 10 * np.mean(ys) + 5, f'{100 * a:.0f}/{100 * b:.0f}', fontsize=5,
            va='center', ha='left', color='#c0392b' if b >= .5 else '#999999')
plu.format_ax(ax=ax, xlabel='Cosine distance', reduced_spines=True)
ax.set_yticks([])
for i, b in enumerate(dn['ivl']):
    r = blk_reg[list(blk_names).index(b)]
    ax.text(-.004, 10 * i + 5, b.replace('_Ventricle', ' V').replace('_septum', ' S')
            .replace('Centre S', 'Centre septum'), ha='left', va='center', fontsize=6.5,
            color=COLORS[r], fontweight='bold', clip_on=False)
ax.text(0, -.16, 'support %: SNV jackknife / cuts resampled within blocks\n(grey: < 50%)\n'
        f'same-region blocks closer than different-region: block-label permutation p = {p_assort:.2f}',
        transform=ax.transAxes, fontsize=5.5, color='#555555', va='top')
ax = fig.add_subplot(gs[0, 1])
leaf_rows = [list(blk_names).index(b) for b in dn['ivl']][::-1]
im = ax.imshow(np.ma.masked_invalid(P_blk[np.ix_(leaf_rows, mut_order)]), cmap=cmap_af, aspect='auto',
               vmin=0, vmax=np.nanpercentile(P_blk, 98))
for e in block_edges:
    ax.axvline(e, color='k', linewidth=.6)
ax.set_yticks([])
group_strip(ax, -1.6, .8)
plu.format_ax(ax=ax, xticks=[], xlabel=f'SNVs (n={len(muts)})')
fig.subplots_adjust(left=.02, right=.95, top=.86, bottom=.26)
pos = ax.get_position()
cax = fig.add_axes((pos.x1 - .12, pos.y0 - .05, .12, .025))
cb = fig.colorbar(im, cax=cax, orientation='horizontal')
cb.set_label('VAF', fontsize=6, labelpad=1)
cb.ax.tick_params(labelsize=5, length=2, pad=1)
fig.savefig(os.path.join(path_figures, 'block_hierarchy.pdf'))


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
    ax.add_patch(plt.Rectangle((-2.9, k - .5), 2, 1, color=COLORS[reg[i]],
                               clip_on=False, linewidth=0))
group_strip(ax, -3, 2)
plu.format_ax(ax=ax, xticks=[], yticks=[], xlabel=f'SNVs (n={len(muts)})',
              ylabel='62 LCM cuts, by region then chunk')
ax.yaxis.labelpad = 14
plu.add_legend(colors=COLORS, label='Region', ax=ax, ticks_size=6, artists_size=5,
               label_size=7, loc='upper center', bbox_to_anchor=(.5, -.06), ncols=5)
fig.subplots_adjust(left=.06, right=.97, top=.86, bottom=.14)
pos = ax.get_position()
cax = fig.add_axes((pos.x1 - .12, pos.y0 - .02 - .02, .12, .02))
cb = fig.colorbar(im, cax=cax, orientation='horizontal')
cb.set_label('VAF', fontsize=6, labelpad=1)
cb.ax.tick_params(labelsize=5, length=2, pad=1)
fig.savefig(os.path.join(path_figures, 'lcm_heatmap.pdf'))


##


# 2b, companion. The same view for the enriched mutations, grouped by the region
# where each is enriched, strongest first
enr_region = region_p.idxmin(0)
enr_groups = []
for r in leaf_order:
    sel = np.array([j for j in np.where(is_enriched)[0] if enr_region.iloc[j] == r], dtype=int)
    enr_groups.append((r, sel[np.argsort(-max_logp[sel])]))
enr_order = np.concatenate([g for _, g in enr_groups]).astype(int)
enr_edges = np.cumsum([len(g) for _, g in enr_groups])[:-1] - .5

fig, ax = plt.subplots(figsize=(7, 5.5))
im = ax.imshow(AF[np.ix_(cut_order, enr_order)], cmap='afmhot_r', aspect='auto',
               vmin=0, vmax=np.nanpercentile(AF, 99))
for e in chunk_edges:
    ax.axhline(e, color='#999999', linewidth=.3)
for e in region_edges:
    ax.axhline(e, color='k', linewidth=1.2)
for e in enr_edges:
    ax.axvline(e, color='k', linewidth=.6)
for k, i in enumerate(cut_order):
    ax.add_patch(plt.Rectangle((-1.2, k - .5), .6, 1, color=COLORS[reg[i]],
                               clip_on=False, linewidth=0))
group_strip(ax, -3, 2, grp=[(r, g) for r, g in enr_groups if len(g)])
plu.format_ax(ax=ax, yticks=[], ylabel='62 LCM cuts, by region then chunk',
              xlabel=f'{len(enr_order)} enriched mutations, by enriched region')
ax.set_xticks(range(len(enr_order)))
ax.set_xticklabels([muts[j].rsplit('_', 2)[0].replace('_', ':') for j in enr_order],
                   rotation=45, ha='right', rotation_mode='anchor', fontsize=6.5)
ax.yaxis.labelpad = 14
cb = fig.colorbar(im, ax=ax, pad=.01, fraction=.03)
cb.set_label('VAF', fontsize=7)
cb.ax.tick_params(labelsize=6)
fig.subplots_adjust(left=.08, right=.92, top=.9, bottom=.2)
fig.savefig(os.path.join(path_figures, 'lcm_heatmap_enriched.pdf'))


##


# 2c. The region distances behind the dendrogram
fig, ax = plt.subplots(figsize=(3.6, 3))
D_show = pd.DataFrame(D_reg, index=REGIONS, columns=REGIONS).loc[leaf_order, leaf_order]
plu.plot_heatmap(D_show, ax=ax, annot=True, fmt='.3f', cb=True)
ax.set_aspect('equal')
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'distances.pdf'))


# 2d. Clade support across profile constructions and lineage classes: colour = block bootstrap
# support, text = SNV jackknife / block bootstrap; blank = clade not in that tree
PROFILE_LABELS = {'pooled': 'Pooled AD/DP', 'equal_block': 'Equal block weight',
                  'mean_cut': 'Mean of cuts'}
SET_ORDER = ['All', 'Pre-gastrulation', 'Heart-specific', 'Other shared']
in_tree = hierarchy.query('in_tree')
col_keys = [(k, s_) for k in PROFILE_LABELS for s_ in SET_ORDER]
clades = (in_tree.groupby('clade').support_block.max().sort_values(ascending=False).index.tolist())
blk_sup = np.full((len(clades), len(col_keys)), np.nan)
labels_sup = np.full(blk_sup.shape, '', dtype=object)
for _, row in in_tree.iterrows():
    yi, xi = clades.index(row.clade), col_keys.index((row.profile, row.snv_set))
    blk_sup[yi, xi] = 100 * row.support_block
    labels_sup[yi, xi] = f'{100 * row.support_snv:.0f}/{100 * row.support_block:.0f}'
fig, ax = plt.subplots(figsize=(9, .35 * len(clades) + 1.8))
im = ax.imshow(np.ma.masked_invalid(blk_sup), cmap=matplotlib.colormaps['Blues'], vmin=0, vmax=100,
               aspect='auto')
for yi in range(len(clades)):
    for xi in range(len(col_keys)):
        if labels_sup[yi, xi]:
            ax.text(xi, yi, labels_sup[yi, xi], ha='center', va='center', fontsize=5.5,
                    color='white' if blk_sup[yi, xi] > 60 else 'k')
for x in (3.5, 7.5):
    ax.axvline(x, color='k', linewidth=.8)
ax.set_yticks(range(len(clades)))
ax.set_yticklabels(clades, fontsize=7)
ax.set_xticks(range(len(col_keys)))
n_set = hierarchy.drop_duplicates('snv_set').set_index('snv_set').n_snvs
ax.set_xticklabels([f'{s_} ({n_set[s_]})' for _, s_ in col_keys], rotation=45, ha='right',
                   rotation_mode='anchor', fontsize=6)
for i, k in enumerate(PROFILE_LABELS):
    ax.text(4 * i + 1.5, -.9, PROFILE_LABELS[k], ha='center', va='bottom', fontsize=7)
ax.set_ylim(len(clades) - .5, -.5)
cb = fig.colorbar(im, ax=ax, pad=.01, fraction=.03)
cb.set_label('Block bootstrap support (%)', fontsize=7)
cb.ax.tick_params(labelsize=6)
ax.set_title('Clades per tree (text: SNV jackknife / block bootstrap %)', fontsize=8, pad=18)
fig.subplots_adjust(left=.16, right=.93, top=.85, bottom=.28)
fig.savefig(os.path.join(path_figures, 'hierarchy_support.pdf'))


##


# 3. Septum/ventricle asymmetry: the statistic, its test and its robustness
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

fig = plt.figure(figsize=(10, 2.9))
gs = fig.add_gridspec(1, 3, wspace=.3)

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

fig.subplots_adjust(left=.07, right=.97, top=.82, bottom=.18)
fig.savefig(os.path.join(path_figures, 'septum_ventricle_asymmetry.pdf'))

# Read support per SNV: mean alternate reads and mean coverage across the LCM samples
from scipy.stats import gaussian_kde
mean_ad = AD_df.mean(0).reindex(muts).values
mean_dp = DP_df.mean(0).reindex(muts).values
n_alt_samples = (AD_df > 0).sum(0).reindex(muts).values.astype(float)
fig, axs = plt.subplots(1, 3, figsize=(9.5, 2.6))
for ax, v, label in [(axs[0], mean_ad, 'Mean AD per LCM sample'),
                     (axs[1], mean_dp, 'Mean DP per LCM sample'),
                     (axs[2], n_alt_samples, 'LCM samples with AD >= 1')]:
    grid = np.linspace(0, v.max() * 1.1, 300)
    ax.fill_between(grid, gaussian_kde(v)(grid), color='#9e9e9e', alpha=.5, linewidth=0)
    ax.plot(grid, gaussian_kde(v)(grid), color='k', linewidth=1)
    ax.plot(v, np.full(len(v), -.03 * gaussian_kde(v)(grid).max()), '|', color='k', markersize=5,
            alpha=.6)
    ax.axvline(np.median(v), color='#c0392b', linestyle='--', linewidth=.8)
    ax.text(np.median(v), ax.get_ylim()[1], f' median {np.median(v):.1f}', color='#c0392b',
            fontsize=6.5, va='top', ha='left')
    plu.format_ax(ax=ax, xlabel=label, ylabel='Density', reduced_spines=True)
fig.suptitle(f'Read support across {len(muts)} SNVs ({len(cuts)} LCM samples)', fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'snv_read_support.pdf'))


# Clade markers and region-exclusive SNVs: region tree and heatmap over the per-sample heatmap
markers_region_and_samples(marker_groups + exclusive_groups, 'markers_exclusive')

print('written: embryonic_lineage_features, genotype_call_composition, '
      'embryonic_lineage_carriers, lineage_enrichment, enrichment_type, example_enriched, example_spread, example_local, lineage_sweep, '
      'phylogenetic_relationships, lcm_heatmap, distances, septum_ventricle_asymmetry')
