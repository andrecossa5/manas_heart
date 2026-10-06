"""
Final figure on an A4 page, and its supplementary panels. Reads results tables only:
results/GENOTYPES_TRUE.tsv.gz, LINEAGE_SUMMARY.tsv, REGION_ENRICHMENT.tsv (5.lineage_analysis.py),
SEPTUM_LV_CUTS.tsv, SEPTUM_LV_CUTS_TESTS.tsv (6.septum_ventricle.py), and the sample coordinates.
Enriched = BH q < 0.1 of the region-vs-rest test; its permutation test is run here.

Main (figures/main_figure.pdf, and each panel alone as main_<letter>.pdf, cut from the
same render so sizes and fonts are identical):
  a  placeholder for the experimental-design cartoon
  b  sampling in 3D: analysed samples by region, the other heart samples in grey
  c  SNV categories from the tree assignment (+ placenta reads)
  d  SNVs ranked by their strongest region enrichment (-log10 p), size = that region's VAF
  e  one example per region and one across the heart, read-count 3D view
  f  septal samples, cosine similarity to pooled LV vs RV (all SNVs)
  g  two SNVs, cell fraction per sample by region

Supplementary: supp_c, supp_d, supp_fi, supp_fii, supp_fiii, supp_fiv.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.transforms import Bbox
from mpl_toolkits.mplot3d.proj3d import proj_transform
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
from scipy.stats import binom
from sklearn.metrics import pairwise_distances
from statsmodels.stats.multitest import multipletests
import plotting_utils as plu

matplotlib.use('Agg')
plu.set_rcParams()
FS_SMALL, FS_TINY = 8, 7                    # plu defaults for ticks; annotations


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
OTHER_COLOR = '#d9d9d9'                      # heart samples not analysed
CATEGORIES = ['Heart', 'Gut+Blood', 'Blood', 'Placenta', 'Gut']      # most to least abundant
CAT_COLORS = {'Heart': '#FF0000', 'Gut+Blood': '#00A08A', 'Blood': '#F2AD00', 'Placenta': '#F98400',
              'Gut': '#5BBCD6'}                  # Darjeeling1
TREE_TO_CAT = {'Gut': 'Gut', 'Blood': 'Blood', 'Blood_shared': 'Blood', 'Blood_Gut': 'Gut+Blood'}
E_ORDER = ['All heart', 'LV', 'RV', 'LS', 'RS', 'CS']
G_ORDER = ['LV', 'LS', 'RS', 'CS', 'RV']
E_OVERRIDE = {}                              # slot -> mutation_id, replaces the automatic pick
G_OVERRIDE = {}                              # 'LV-LS/RS' or 'RV-CS' -> mutation_id
B_POINT = 22                                 # 3D point sizes: panel b, panel e
E_POINT = 10
VAF_MAX_3D = .2
A4 = (8.27, 11.69)


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


def style_3d(ax, xyz, elev=20, azim=50):
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
    ax.view_init(elev=elev, azim=azim)
    lo, hi = xyz.min(), xyz.max()
    ends = {
        'x': (hi['x'] + (hi['x'] - lo['x']) * .15, lo['y'], lo['z']),
        'y': (lo['x'], hi['y'] + (hi['y'] - lo['y']) * .15, lo['z']),
        'z': (lo['x'], lo['y'], hi['z'] + (hi['z'] - lo['z']) * .15),
    }
    for label, (xe, ye, ze) in ends.items():
        ax.add_artist(Arrow3D(lo['x'], lo['y'], lo['z'], xe, ye, ze, mutation_scale=6,
                              lw=.4, arrowstyle='-|>', color='k', shrinkA=0, shrinkB=0))
        ax.text(lo['x'] + (xe - lo['x']) * 1.07, lo['y'] + (ye - lo['y']) * 1.07,
                lo['z'] + (ze - lo['z']) * 1.07, label, fontsize=6, ha='center', va='center')


def site(m):
    return m.rsplit('_', 2)[0].replace('_', ':')


def p_text(p):
    return 'p < 0.001' if p < .001 else f'p = {p:.3f}'


##


path_main = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))     # repository root
path_results = os.path.join(path_main, 'results')
path_out = os.path.join(path_main, 'figures')
os.makedirs(path_out, exist_ok=True)

geno = pd.read_csv(os.path.join(path_results, 'GENOTYPES_TRUE.tsv.gz'), sep='\t')
summary = pd.read_csv(os.path.join(path_results, 'LINEAGE_SUMMARY.tsv'), sep='\t').set_index('mutation_id')
enrich = pd.read_csv(os.path.join(path_results, 'REGION_ENRICHMENT.tsv'), sep='\t')
lv_cuts = pd.read_csv(os.path.join(path_results, 'SEPTUM_LV_CUTS.tsv'), sep='\t')
lv_tests = pd.read_csv(os.path.join(path_results, 'SEPTUM_LV_CUTS_TESTS.tsv'), sep='\t')
xyz_all = pd.read_csv(os.path.join(path_main, 'data', 'Heart_final_coorindates_135.csv')).set_index('name')

geno['reg'] = geno['region'].map(REGION_ABBR)
AF_df = geno.pivot(index='Sample_ID', columns='mutation_id', values='AF')
AD_df = geno.pivot(index='Sample_ID', columns='mutation_id', values='AD_alt')
DP_df = geno.pivot(index='Sample_ID', columns='mutation_id', values='DP')
STATE = geno.pivot(index='Sample_ID', columns='mutation_id', values='state_af10')
cuts, muts = AF_df.index, AF_df.columns
AF = AF_df.values
meta = geno.drop_duplicates('Sample_ID').set_index('Sample_ID').loc[cuts]
reg = meta['reg'].values
present = (STATE == 'present').values
coords = xyz_all.loc[cuts, ['x', 'y', 'z']]
analysed = xyz_all.index.isin(cuts)

# SNV categories: tree assignment; unassigned SNVs with placenta reads above error (rule of
# 5.lineage_analysis.py: binomial at error 5e-4, p < 0.01, placenta VAF > 20% of heart VAF) are Placenta
heart_af = (AD_df.sum() / DP_df.sum()).reindex(muts).values
s = summary.reindex(muts)
placenta_af = s['placenta_AD'] / s['placenta_DP'].replace(0, np.nan)
placenta_signal = ((binom.sf(s['placenta_AD'] - 1, s['placenta_DP'], 5e-4) < .01)
                   & (placenta_af.fillna(0).values > .2 * heart_af))
category = np.where(s['tree_assignment'] == 'Unassigned',
                    np.where(placenta_signal, 'Placenta', 'Heart'),
                    s['tree_assignment'].map(TREE_TO_CAT).values)
cat_counts = pd.Series(category).value_counts().reindex(CATEGORIES).fillna(0).astype(int)
assert cat_counts.sum() == len(muts) == 123

# Region x mutation tables
region_tab = lambda col: enrich.pivot_table(index='region', columns='mutation_id', values=col).reindex(REGIONS)[muts]
region_af = region_tab('AF')
region_p = region_tab('p_binom')
region_frac = region_tab('frac_cuts_alt')
ENRICH_Q = .1
region_call = (region_tab('q_binom') < ENRICH_Q).fillna(False)              # enriched: BH q < 0.1
logp = -np.log10(region_p.clip(lower=1e-30))
is_enriched = region_call.any().values

# Global significance of the number of enriched SNVs: the same rule (read-level binomial of each
# region against the rest of the heart, rest VAF floored at the site background, BH over the
# 5 x 123 tests) rerun on sample-label permutations across regions
A_c, D_c = AD_df.values.astype(float), DP_df.values.astype(float)
site_bg = summary['site_background'].reindex(muts).values
N_PERM = 1000


def enrichment_calls(labels):
    p = []
    for r in REGIONS:
        k = labels == r
        a, d = A_c[k].sum(0), D_c[k].sum(0)
        p.append(binom.sf(a - 1, d, np.maximum((A_c.sum(0) - a) / (D_c.sum(0) - d), site_bg)))
    return multipletests(np.ravel(p), method='fdr_bh')[1].reshape(len(REGIONS), -1) < ENRICH_Q


assert (enrichment_calls(reg) == region_call.values).all(), 'does not reproduce REGION_ENRICHMENT.tsv'
perm_rng = np.random.default_rng(1)
n_null = np.array([enrichment_calls(perm_rng.permutation(reg)).any(0).sum() for _ in range(N_PERM)])
enrich_perm = pd.Series({'n_enriched_snvs': int(is_enriched.sum()), 'n_permutations': N_PERM,
                         'expected_by_chance': n_null.mean(),
                         'p_global': (np.sum(n_null >= is_enriched.sum()) + 1) / (N_PERM + 1)})

print(f'{len(muts)} SNVs x {len(cuts)} cuts | categories: ' +
      ', '.join(f'{k} {v}' for k, v in cat_counts.items()))
print(f'  not Heart (pre-date heart formation): {len(muts) - cat_counts["Heart"]}/{len(muts)}')
print(f'  enriched SNVs: {int(is_enriched.sum())}/{len(muts)} ({enrich_perm.expected_by_chance:.1f} expected, '
      f'permutation p = {enrich_perm.p_global:.4f}, {int(region_call.values.sum())} SNV x region hits)')


##


# Example picks. e: per region, the SNV enriched there with the lowest p among those with alt
# reads in > 30% of the region's samples; All heart, the non-enriched SNV with alt reads in most
# samples (ties: higher heart VAF). g: SNVs whose pooled VAF is high in one set of regions and
# low in the other, score = min over the high set - max over the low set.
n_carriers = (AD_df.values > 0).sum(0)
by_carriers = np.lexsort((-heart_af, -n_carriers))
e_candidates = {'All heart': list(muts[by_carriers][~is_enriched[by_carriers]])}
for r in REGIONS:
    i = REGIONS.index(r)
    sel = np.where(region_call.values[i] & (region_frac.values[i] > .3))[0]
    e_candidates[r] = list(muts[sel[np.argsort(-logp.values[i, sel])]])
G_SETS = {'LV-LS/RS': (['LV', 'LS', 'RS'], ['RV', 'CS']), 'RV-CS': (['RV', 'CS'], ['LV', 'LS', 'RS'])}
g_candidates = {}
for k, (high, low) in G_SETS.items():
    score = (region_af.loc[high].min(0) - region_af.loc[low].max(0)).values
    g_candidates[k] = list(muts[np.argsort(-score)])
E_PICK = {k: E_OVERRIDE.get(k, e_candidates[k][0]) for k in E_ORDER}
G_PICK = {k: G_OVERRIDE.get(k, g_candidates[k][0]) for k in G_SETS}
print('\nexample candidates (top 5; first is used unless overridden)')
for k in E_ORDER:
    print(f'  e {k}: ' + ', '.join(e_candidates[k][:5]))
for k in G_SETS:
    print(f'  g {k}: ' + ', '.join(g_candidates[k][:5]))


##


# Panels. Each draws into axes it is given.
def panel_placeholder(ax):
    ax.set_xticks([]); ax.set_yticks([])
    ax.text(.5, .5, 'Experimental design\n(cartoon)', ha='center', va='center', fontsize=FS_SMALL, color='#999999',
            transform=ax.transAxes)
    for sp in ax.spines.values():
        sp.set(color='#bbbbbb', linewidth=.6)


def panel_b(ax):
    """
    Analysed samples by region, the other heart samples in grey behind them.
    """
    ax.computed_zorder = False
    o = xyz_all.loc[~analysed]
    ax.scatter(o['x'], o['y'], o['z'], s=B_POINT, c=OTHER_COLOR, edgecolor='white', linewidth=.25,
               depthshade=False, zorder=4)
    ax.scatter(coords['x'], coords['y'], coords['z'], s=B_POINT, c=[COLORS[r] for r in reg],
               edgecolor='white', linewidth=.25, depthshade=False, zorder=5)
    style_3d(ax, xyz_all[['x', 'y', 'z']])
    plu.add_legend(colors={**COLORS, 'Other': OTHER_COLOR}, label='Region', ax=ax, ticks_size=FS_TINY,
                   artists_size=FS_TINY - 1, label_size=FS_TINY, loc='center left', bbox_to_anchor=(.97, .5), ncols=1)


def panel_c(ax):
    ax.bar(range(len(CATEGORIES)), cat_counts.values, color=[CAT_COLORS[k] for k in CATEGORIES], width=.85,
           edgecolor='k', linewidth=.5)
    for i, v in enumerate(cat_counts.values):
        ax.text(i, v + 1, str(v), ha='center', fontsize=FS_SMALL)
    plu.format_ax(ax=ax, xticks=CATEGORIES, ylabel='n SNVs', rotx=40, reduced_spines=True)
    plt.setp(ax.get_xticklabels(), ha='right', rotation_mode='anchor')


best_region = logp.values.argmax(0)                 # per SNV, its most enriched region
d_tab = pd.DataFrame({
    'region': [REGIONS[i] for i in best_region],
    'logp': logp.values[best_region, np.arange(len(muts))],
    'vaf': region_af.values[best_region, np.arange(len(muts))],
    'hit': region_call.values[best_region, np.arange(len(muts))],
}).sort_values('logp', ascending=False).reset_index(drop=True)
VAF_SIZE = lambda v: 4 + 96 * (np.clip(v, .01, .2) - .01) / .19         # 1% -> 20% VAF


SIG_LOGP = -np.log10(enrich.loc[enrich['q_binom'] < ENRICH_Q, 'p_binom'].max())      # weakest p with BH q < 0.1


def panel_d(ax):
    """
    SNVs ranked by their strongest region-vs-rest enrichment; red = enriched there
    (BH q < 0.1 on the p-value alone), black = not; size = that region's pooled VAF (1-20%).
    Red line: the p-value at which q reaches 0.1.
    """
    for hit, col in [(False, 'k'), (True, '#c0392b')]:
        t = d_tab[d_tab.hit == hit]
        ax.scatter(t.index, t.logp, s=VAF_SIZE(t.vaf.values), color=col, edgecolor='none', alpha=.5, zorder=3)
    ax.axhline(SIG_LOGP, color='#c0392b', linestyle='--', linewidth=.8, zorder=1)
    ax.text(.02, SIG_LOGP - .35, f'FDR 10%, p = {10 ** -SIG_LOGP:.3f}', transform=ax.get_yaxis_transform(),
            ha='left', va='top', fontsize=FS_TINY, color='#c0392b')
    hs = [ax.scatter([], [], s=VAF_SIZE(v), color='#888888', edgecolor='none') for v in (.01, .05, .1, .2)]
    leg = ax.legend(hs, ['1%', '5%', '10%', '20%'], title='Region VAF', frameon=False, loc='upper right',
                    bbox_to_anchor=(1, 1), fontsize=FS_TINY, title_fontsize=FS_TINY, labelspacing=.6)
    ax.add_artist(leg)
    hs = [ax.scatter([], [], s=VAF_SIZE(.05), color=c, edgecolor='none') for c in ('#c0392b', 'k')]
    ax.legend(hs, ['enriched', 'not enriched'], title='Status', frameon=False, loc='upper right',
              bbox_to_anchor=(.62, 1), fontsize=FS_TINY, title_fontsize=FS_TINY, labelspacing=.6)
    plu.format_ax(ax=ax, xlabel='SNVs', ylabel='Max enrichment score', reduced_spines=True,
                  title=f'{int(is_enriched.sum())}/{len(muts)} enriched SNVs')
    ins = ax.inset_axes([.5, .33, .48, .2])                     # enriched SNVs per region (strongest region)
    n_reg = d_tab[d_tab.hit].region.value_counts().reindex(REGIONS).fillna(0).astype(int).sort_values(ascending=False)
    ins.bar(range(len(REGIONS)), n_reg.values, color='#c0392b', alpha=.5, width=.7)
    for i, v in enumerate(n_reg.values):
        ins.text(i, v, str(v), ha='center', va='bottom', fontsize=FS_TINY - 1)
    plu.format_ax(ax=ins, xticks=list(n_reg.index), ylabel='n', reduced_spines=True, xticks_size=FS_TINY - 1,
                  yticks_size=FS_TINY - 1, ylabel_size=FS_TINY)


def draw_counts_3d(ax, j):
    """
    Samples with >=1 alternate read filled by VAF, others open.
    """
    carries = AD_df.values[:, j] > 0
    ax.scatter(coords['x'][~carries], coords['y'][~carries], coords['z'][~carries], s=E_POINT,
               facecolor='white', edgecolor='k', linewidth=.3, depthshade=False, zorder=4)
    ax.scatter(coords['x'][carries], coords['y'][carries], coords['z'][carries], s=E_POINT,
               c=AF[carries, j], cmap='afmhot_r', vmin=0, vmax=VAF_MAX_3D,
               edgecolor='k', linewidth=.3, depthshade=False, zorder=5)
    style_3d(ax, coords)


def panel_e(axs):
    for ax, slot in zip(axs, E_ORDER):
        ax.computed_zorder = False
        j = list(muts).index(E_PICK[slot])
        draw_counts_3d(ax, j)
        carries = AD_df.values[:, j] > 0
        ax.set_title(f'{site(muts[j])} ({slot})\nVAF {100 * AF[carries, j].mean():.1f}%, n={int(carries.sum())} samples',
                     fontsize=FS_TINY, pad=0)
    plu.add_cbar(np.array([0, VAF_MAX_3D]), ax=axs[5], label='VAF', palette='afmhot_r', vmin=0, vmax=VAF_MAX_3D,
                 label_size=FS_TINY, ticks_size=FS_TINY - 1)


SET_LABELS = {'All SNVs': 'All SNVs', 'Pre-gastrulation + Other shared': 'Shared',
              'Heart-specific': 'Heart-specific', 'Physical space': 'Physical space'}
SHOW_GROUPS = ['Septum', 'LS', 'CS', 'RS']
_mol = lv_cuts[lv_cuts.snv_set != 'Physical space']
_sim = 1 - _mol[['d_LV', 'd_RV']].values                       # cosine similarity = 1 - cosine distance
LIM_MOL = (np.floor(_sim.min() * 20) / 20, np.ceil(_sim.max() * 20) / 20)
_sp = lv_cuts[lv_cuts.snv_set == 'Physical space']
LIM_SP = (np.floor((_sp[['d_LV', 'd_RV']].values.min() - 30) / 50) * 50,
          np.ceil((_sp[['d_LV', 'd_RV']].values.max() + 30) / 50) * 50)      # padded: no point on the axes


ROW_PITCH = 8.6                                 # points between the rows of the stats block (as in main f)


def septum_scatter(ax, set_name, legend=True, x_text=.44):
    d = lv_cuts[lv_cuts.snv_set == set_name]
    lim = LIM_SP if set_name == 'Physical space' else LIM_MOL
    ax.plot(lim, lim, color='k', linewidth=.6, linestyle='--', zorder=1)
    physical = set_name == 'Physical space'
    x, y = (d.d_RV, d.d_LV) if physical else (1 - d.d_RV, 1 - d.d_LV)           # molecular: cosine similarity
    ax.scatter(x, y, s=20, c=[COLORS[r] for r in d.region], edgecolor='k', linewidth=.3, zorder=3)
    sub = lv_tests[lv_tests.snv_set == set_name].set_index('septal_group')
    tests = [sub.loc[g] for g in SHOW_GROUPS]
    lines = [('Closer to LV:', 'k')] + [
        (f'{g} {int(r.cuts_closer_LV)}/{int(r.n_cuts)}, {p_text(r.p_wilcoxon)}', COLORS.get(g, 'k'))
        for g, r in zip(SHOW_GROUPS, tests)]
    for i, (txt, col) in enumerate(lines):      # molecular: above the diagonal = closer to LV, text lower right
        xy, dy, va = ((.1, .9), -ROW_PITCH * i, 'top') if physical else \
                     ((x_text, .06), ROW_PITCH * (len(lines) - 1 - i), 'bottom')   # physical: below = closer, upper left
        ax.annotate(txt, xy=xy, xycoords='axes fraction', xytext=(0, dy), textcoords='offset points',
                    fontsize=FS_TINY - 1, ha='left', va=va, color=col)
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect('equal')
    unit = 'Mean Euclidean distance' if physical else 'Similarity'
    n = '' if physical else f' (n={int(d.n_snvs.iloc[0])})'
    plu.format_ax(ax=ax, xlabel=f'{unit} to RV', ylabel=f'{unit} to LV', title=f'{SET_LABELS[set_name]}{n}',
                  reduced_spines=True)
    if legend:
        plu.add_legend(colors={r: COLORS[r] for r in ['LS', 'CS', 'RS']}, label='Septal section', ax=ax,
                       ticks_size=FS_TINY - 1, artists_size=FS_TINY - 2, label_size=FS_TINY - 1, loc='lower right' if physical else 'upper left',
                       bbox_to_anchor=(1, .02) if physical else (0, .98))


def panel_g(axs):
    for k, (ax, key) in enumerate(zip(axs, G_SETS)):
        j = list(muts).index(G_PICK[key])
        df = pd.DataFrame({'region': reg, 'cf': 2 * AF[:, j]})
        plu.box(df, x='region', y='cf', x_order=G_ORDER, color='white', width=.6, ax=ax)
        plu.strip(df, x='region', y='cf', x_order=G_ORDER, categorical_cmap=COLORS, size=5, ax=ax)
        plu.format_ax(ax=ax, xticks=G_ORDER, xlabel='', ylabel='Cell fraction (2 x VAF)' if k == 0 else '',
                      title=site(muts[j]), reduced_spines=True)


##


# Main figure on an A4 page, laid out in inches from the top-left corner (proportions of
# mock_figure.pdf); the content fills the page from the top, not its full length.
W, H = A4


def add_ax(fig, x, top, w, h, **kw):
    return fig.add_axes([x / W, 1 - (top + h) / H, w / W, h / H], **kw)


T1, T2, T3 = .45, 3.0, 6.5                     # row tops
fig = plt.figure(figsize=A4)
panels = {
    'a': [add_ax(fig, .55, T1 + .1, 1.45, 1.9)],
    'b': [add_ax(fig, 1.95, T1 - .15, 3.5, 2.85, projection='3d')],
    'c': [add_ax(fig, 6.45, T1 + .15, 1.5, 1.55)],
    'd': [add_ax(fig, .8, T2 + .4, 2.1, 2.4)],
    'e': [add_ax(fig, 3.12 + 1.5 * (k % 3), T2 + .2 + 1.6 * (k // 3), 1.5, 1.4, projection='3d') for k in range(6)],
    'f': [add_ax(fig, .8, T3 + .3, 2.0, 2.0)],
    'g': [add_ax(fig, 3.9, T3 + .3, 1.85, 2.0), add_ax(fig, 6.1, T3 + .3, 1.85, 2.0)],
}
panel_placeholder(*panels['a'])
panel_b(*panels['b'])
panel_c(*panels['c'])
panel_d(*panels['d'])
panel_e(panels['e'])
septum_scatter(*panels['f'], 'All SNVs')
panel_g(panels['g'])
LETTERS = {'a': (.2, T1), 'b': (2.3, T1), 'c': (5.85, T1), 'd': (.2, T2), 'e': (3.12, T2), 'f': (.2, T3), 'g': (3.2, T3)}
for letter, (x, top) in LETTERS.items():
    fig.text(x / W, 1 - top / H, letter, fontsize=12, fontweight='bold', va='top')
fig.savefig(os.path.join(path_out, 'main_figure.pdf'))
fig.savefig(os.path.join(path_out, 'main_figure.png'), dpi=200)

# Each panel alone: the same render, cropped to the panel's artists
renderer = fig.canvas.get_renderer()
for letter, axs in panels.items():
    bb = Bbox.union([ax.get_tightbbox(renderer) for ax in axs]).transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(os.path.join(path_out, f'main_{letter}.pdf'), bbox_inches=bb.padded(.05))
plt.close(fig)


##


# Supp c. Carrier samples and pooled VAF per category: box (white) and strip as in main g, colours as c
n_present = present.sum(0).astype(float)


def category_box(ax, values, ylabel, first):
    df = pd.DataFrame({'category': category, 'v': values})
    plu.box(df, x='category', y='v', x_order=CATEGORIES, color='white', width=.6, ax=ax)
    plu.strip(df, x='category', y='v', x_order=CATEGORIES, categorical_cmap=CAT_COLORS, size=5, ax=ax)
    plu.format_ax(ax=ax, xticks=CATEGORIES, xlabel='', ylabel=ylabel, rotx=40, reduced_spines=True)
    plt.setp(ax.get_xticklabels(), ha='right', rotation_mode='anchor')


W_C, H_C = 6.2, 3.2
fig = plt.figure(figsize=(W_C, H_C))
axs = [fig.add_axes([x / W_C, .85 / H_C, 2.3 / W_C, 2.1 / H_C]) for x in (.75, 3.75)]
category_box(axs[0], n_present, 'n carrier samples', True)
category_box(axs[1], heart_af, 'Pooled VAF, all heart', False)
fig.savefig(os.path.join(path_out, 'supp_c.pdf'))
plt.close(fig)


# Supp d. Pooled region VAF, % samples with alt reads, and the enrichment test; SNVs by pooled heart VAF
cols_af = np.argsort(-np.nan_to_num(heart_af))
grey_bad = lambda name: matplotlib.colormaps[name].with_extremes(bad='#dddddd')
fig, axs = plt.subplots(3, 1, figsize=(8.27, 4.8), sharex=True)
for ax, values, cmap, vmax, label, title in [
        (axs[0], region_af.values, grey_bad('afmhot_r'), np.nanpercentile(region_af.values, 98), 'VAF',
         'Region VAF (pooled AD / DP)'),
        (axs[1], 100 * region_frac.values, grey_bad('Blues'), 100, '% samples', '% samples with >0 alt reads'),
        (axs[2], logp.values, grey_bad('viridis'), np.nanpercentile(logp.values, 99), '-log10 p',
         'Region vs rest of heart, one-sided binomial (* enriched)')]:
    im = ax.imshow(np.ma.masked_invalid(values[:, cols_af]), cmap=cmap, aspect='auto', vmin=0, vmax=vmax)
    for yi in range(len(REGIONS)):
        for xi, j in enumerate(cols_af):
            if region_call.values[yi, j]:
                rgb = np.array(im.cmap(im.norm(values[yi, j])))[:3]
                ax.text(xi, yi, '*', ha='center', va='center', fontsize=FS_SMALL,
                        color='white' if rgb @ [.299, .587, .114] < .5 else 'k')
    cb = fig.colorbar(im, ax=ax, pad=.005, fraction=.015)
    cb.set_label(label, fontsize=FS_SMALL)
    cb.ax.tick_params(labelsize=FS_TINY)
    plu.format_ax(ax=ax, yticks=REGIONS, xticks=[], title=title)
axs[2].set_xticks(range(len(cols_af)))
axs[2].set_xticklabels([site(m) for m in muts[cols_af]], rotation=90, fontsize=4, ha='center', va='top')
for lab, j in zip(axs[2].get_xticklabels(), cols_af):
    lab.set_fontweight('bold' if is_enriched[j] else 'normal')
axs[2].tick_params(axis='x', length=2, width=.4, pad=1.5)
fig.subplots_adjust(left=.04, right=.95, top=.95, bottom=.2, hspace=.3)
fig.savefig(os.path.join(path_out, 'supp_d.pdf'))
plt.close(fig)


# Supp fi. Shared and heart-specific SNVs, as panel f
fig, axs = plt.subplots(1, 2, figsize=(6, 3))
septum_scatter(axs[0], 'Pre-gastrulation + Other shared', legend=False, x_text=.5)
septum_scatter(axs[1], 'Heart-specific', x_text=.5)
fig.subplots_adjust(left=.09, right=.98, top=.92, bottom=.14, wspace=.3)
fig.savefig(os.path.join(path_out, 'supp_fi.pdf'))
plt.close(fig)


# Supp fii. Share of septal samples closer to LV / RV, per SNV set
STACK_SETS = ['All SNVs', 'Heart-specific', 'Pre-gastrulation + Other shared']
fig, axs = plt.subplots(1, 3, figsize=(6.5, 2.6), sharey=True)
for c, set_name in enumerate(STACK_SETS):
    ax = axs[c]
    for i, g in enumerate(SHOW_GROUPS):
        row = lv_tests.query('snv_set == @set_name and septal_group == @g').iloc[0]
        n, n_lv = int(row.n_cuts), int(row.cuts_closer_LV)
        for bottom, n_, col, lab in [(0, n_lv, COLORS['LV'], 'Closer to LV'), (n_lv / n, n - n_lv, COLORS['RV'], 'Closer to RV')]:
            ax.bar(i, n_ / n, bottom=bottom, color=col, width=.7, edgecolor='white', linewidth=.6,
                   label=lab if (c == 0 and i == 0) else None)
            if n_:
                ax.text(i, bottom + n_ / n / 2, str(n_), ha='center', va='center', fontsize=FS_SMALL, color='white')
    ax.axhline(.5, color='k', linestyle='--', linewidth=.6)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, .25, .5, .75, 1])
    ax.set_yticklabels(['0', '25', '50', '75', '100'])
    n_snvs = int(lv_cuts.query('snv_set == @set_name').n_snvs.iloc[0])
    plu.format_ax(ax=ax, xticks=[f'{g}\n(n={int(lv_tests.query("snv_set == @set_name and septal_group == @g").n_cuts.iloc[0])})'
                                 for g in SHOW_GROUPS],
                  ylabel='% septal samples' if c == 0 else None, title=f'{SET_LABELS[set_name]} (n={n_snvs})',
                  reduced_spines=True)
axs[0].legend(frameon=False, fontsize=FS_SMALL, loc='lower left', bbox_to_anchor=(0, 1.1), ncols=2)
fig.subplots_adjust(left=.09, right=.98, top=.8, bottom=.17, wspace=.1)
fig.savefig(os.path.join(path_out, 'supp_fii.pdf'))
plt.close(fig)


# Supp fiii. Physical distance
fig, ax = plt.subplots(figsize=(3.8, 3.6))
septum_scatter(ax, 'Physical space')
fig.subplots_adjust(left=.22, right=.96, top=.92, bottom=.17)
fig.savefig(os.path.join(path_out, 'supp_fiii.pdf'))
plt.close(fig)


# Supp fiv. Samples x SNVs VAF, both clustered (average linkage, cosine)
row_order = leaves_list(linkage(squareform(pairwise_distances(AF, metric='cosine'), checks=False), method='average'))
col_order = leaves_list(linkage(squareform(pairwise_distances(AF.T, metric='cosine'), checks=False), method='average'))
fig, ax = plt.subplots(figsize=(8.27, 4.6))
im = ax.imshow(AF[np.ix_(row_order, col_order)], cmap='afmhot_r', aspect='auto', vmin=0, vmax=np.nanpercentile(AF, 99))
for k, i in enumerate(row_order):
    ax.add_patch(plt.Rectangle((-4.2, k - .5), 3, 1, color=COLORS[reg[i]], clip_on=False, linewidth=0))
for k, j in enumerate(col_order):
    ax.add_patch(plt.Rectangle((k - .5, -4.2), 1, 3, color=CAT_COLORS[category[j]], clip_on=False, linewidth=0))
plu.format_ax(ax=ax, xticks=[], yticks=[], xlabel=f'SNVs (n={len(muts)})', ylabel=f'LCM samples (n={len(cuts)})')
ax.yaxis.labelpad = 22
plu.add_legend(colors=COLORS, label='Region', ax=ax, ticks_size=FS_SMALL, artists_size=FS_SMALL - 1, label_size=FS_SMALL,
               loc='upper left', bbox_to_anchor=(1.01, 1))
plu.add_legend(colors=CAT_COLORS, label='SNV category', ax=ax, ticks_size=FS_SMALL, artists_size=FS_SMALL - 1,
               label_size=FS_SMALL, loc='upper left', bbox_to_anchor=(1.01, .55))
fig.subplots_adjust(left=.07, right=.86, top=.92, bottom=.17)
pos = ax.get_position()                         # small horizontal VAF bar, right edge on the heatmap's
cax = fig.add_axes([pos.x1 - .14, pos.y0 - .085, .14, .022])
cb = fig.colorbar(im, cax=cax, orientation='horizontal')
cb.set_label('VAF', fontsize=FS_TINY, labelpad=1)
cb.ax.tick_params(labelsize=FS_TINY - 1, length=2, pad=1)
fig.savefig(os.path.join(path_out, 'supp_fiv.pdf'))
plt.close(fig)

print('\nwritten to figures/final: main_figure (+ main_a..g), '
      'supp_c, supp_d, supp_fi, supp_fii, supp_fiii, supp_fiv')
