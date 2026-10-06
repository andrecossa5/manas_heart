"""
Lineage characterisation in the developing heart, and phylogenetic relationships
between heart regions.

Supersedes 7.genotyping.py for the TRUE mutation definition fixed in PLAN.md:
>=3 alt reads in >=2 LCM cuts, with >=5 alt reads summed over one chunk holding
such a cut; pericentromeric sites dropped.

Outputs (results/): TRUE_MUTATIONS.tsv, GENOTYPES_TRUE.tsv.gz, AF_MATRIX.tsv.gz,
ENRICHMENT_PERMUTATION.tsv, HIERARCHY.tsv, REGION_DISTANCES.tsv,
LINEAGE_SUMMARY.tsv, REGION_ENRICHMENT.tsv, REGION_CONTRASTS.tsv,
VENTRICULAR_LEAN.tsv, VENTRICULAR_LEAN_LOO.tsv.
"""

import os
import itertools
import numpy as np
import pandas as pd
from scipy.stats import binom, betabinom, norm, wilcoxon, mannwhitneyu
from scipy.optimize import minimize_scalar
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import pairwise_distances
from statsmodels.stats.multitest import multipletests


##


# Approximate hg38 centromere intervals (Mb), padded by PAD. Acrocentric p-arms
# (13, 14, 15, 21, 22) counted from 0. Replace with the UCSC gap track if exact
# boundaries matter.
CENTROMERES = {
    'chr1': (121.7, 125.1), 'chr2': (91.8, 96.0), 'chr3': (90.5, 93.7),
    'chr4': (49.7, 51.8), 'chr5': (46.5, 50.1), 'chr6': (58.5, 59.8),
    'chr7': (58.1, 61.0), 'chr8': (44.0, 45.9), 'chr9': (43.0, 45.5),
    'chr10': (39.6, 41.6), 'chr11': (51.1, 54.4), 'chr12': (34.7, 37.2),
    'chr13': (0, 18.1), 'chr14': (0, 18.2), 'chr15': (0, 19.7),
    'chr16': (36.3, 38.3), 'chr17': (22.7, 26.9), 'chr18': (15.4, 20.9),
    'chr19': (24.4, 27.2), 'chr20': (25.6, 30.4), 'chr21': (0, 13.0),
    'chr22': (0, 15.1), 'chrX': (58.1, 61.0), 'chrY': (0, 60),
}
PAD = 2.0

MIN_READS_CUT = 3       # alt reads making a cut count towards the TRUE rule
MIN_CUTS = 2            # cuts with MIN_READS_CUT needed
MIN_CHUNK_READS = 5     # alt reads summed over a chunk holding such a cut
ERR_INIT = 5e-4         # error rate used only to flag placenta-carrying sites
PRIOR_BASES = 1000      # shrinkage weight of the class rate, in bases
CLASS_PRIOR_BASES = 1000   # shrinkage of a class rate toward the pooled rate
MIN_DP_PRESENT = 10        # a present call needs this depth
PLACENTA_AF_FRACTION = 0.2
P_PRESENT = 0.01
MIN_AD_PRESENT = 2
AF_EXCLUDE = [0.05, 0.10, 0.20]     # absence stringencies carried through
AF_MAIN = 0.10                      # headline threshold
N_PERM = 2000
N_PERM_JOINT = 20000    # permutations of the three septal labels
N_BOOT = 1000
CHAR_FRACTION = 0.8     # characters resampled per bootstrap replicate
PLOIDY = 2              # donor is female (chrX depth ratio 1.02), so 2 everywhere

# Excluded from every analysis: total depth 20,085 reads over the 62 heart cuts, about ten
# times any other site (median 1,643), so most likely a collapsed repeat / mapping artefact
EXCLUDED = ['chr1_143233282_C_T']

REGION_ABBR = {
    'Left_septum': 'LS', 'Centre_septum': 'CS', 'Right_septum': 'RS',
    'Left_Ventricle': 'LV', 'Right_Ventricle': 'RV',
}
REGIONS = ['LS', 'CS', 'RS', 'LV', 'RV']


##


def is_pericentromeric(mutation_id):
    """
    True if the mutation falls within PAD Mb of a centromere interval.
    """
    chrom, pos = mutation_id.split('_')[0], int(mutation_id.split('_')[1]) / 1e6
    start, end = CENTROMERES.get(chrom, (-9, -9))
    return start - PAD <= pos <= end + PAD


def detection_limit(dp, power=0.95):
    """
    AF detectable with `power` probability by >=1 read at depth dp.
    """
    return 1 - (1 - power) ** (1 / np.maximum(dp, 1))


def perm_pvalue(null, obs):
    """
    Two-sided permutation p-value.
    """
    null = np.asarray(null)
    return min(1, 2 * min((null >= obs).mean(), (null <= obs).mean()))


def fit_concentration(alt, depth, mu):
    """
    MLE of the beta-binomial concentration s (alpha = mu*s, beta = (1-mu)*s)
    for one substitution class. Large s means no overdispersion beyond binomial.
    """
    alt, depth = np.asarray(alt, float), np.asarray(depth, float)

    def nll(log_s):
        s = np.exp(log_s)
        return -betabinom.logpmf(alt, depth, mu * s, (1 - mu) * s).sum()

    res = minimize_scalar(nll, bounds=(np.log(10), np.log(1e7)), method='bounded')
    return float(np.exp(res.x))


def tail_p(alt, depth, mu, s):
    """
    P(X >= alt) under beta-binomial with mean mu and concentration s, per cell.
    """
    return betabinom.sf(alt - 1, depth, mu * s, (1 - mu) * s)


def clade_set(Z, labels):
    """
    Groups of an average-linkage tree, as (leaf set, merge height) pairs, read
    off the linkage matrix. fcluster(maxclust) is not usable here: it returns at
    most k clusters, so on five leaves it silently drops the first merge.
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


##


# Paths
path_main = '/Users/cossa/Desktop/projects/manas_heart'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results')

df = pd.read_csv(os.path.join(path_results, 'ALLELIC_TABLE_FINAL.tsv.gz'), sep='\t')
muts = df.drop_duplicates('mutation_id').set_index('mutation_id')
samples = df.drop_duplicates('Sample_ID').set_index('Sample_ID')[['tissue', 'region', 'chunk']]

full = pd.read_csv(os.path.join(path_results, 'ALLELIC_TABLE.tsv.gz'), sep='\t')
full['mutation_id'] = (
    full['CHROM'] + '_' + full['POS'].astype(str) + '_' + full['REF'] + '_' + full['ALT']
)
full = (
    full
    .loc[full['mutation_id'].isin(muts.index) & full['Sample_ID'].isin(samples.index)]
    .drop_duplicates(['Sample_ID', 'mutation_id'])
)

heart = full.query('tissue == "heart"')
AD_all = heart.pivot(index='Sample_ID', columns='mutation_id', values='AD_alt')
DP_all = heart.pivot(index='Sample_ID', columns='mutation_id', values='DP')
region = samples['region'].loc[AD_all.index].map(REGION_ABBR)
chunk = samples['chunk'].loc[AD_all.index]
print(f'Heart cuts {len(AD_all)} in {chunk.nunique()} chunks | candidate mutations {AD_all.shape[1]}')


##


# 1. TRUE mutation set
strong = AD_all >= MIN_READS_CUT
chunk_totals = AD_all.groupby(chunk).sum()

passing = []
for m in AD_all.columns:
    cuts = strong.index[strong[m]]
    if len(cuts) < MIN_CUTS:
        continue
    if any(chunk_totals.loc[c, m] >= MIN_CHUNK_READS for c in chunk.loc[cuts].unique()):
        passing.append(m)

peri = pd.Series({m: is_pericentromeric(m) for m in passing})
true_muts = pd.Index([m for m in passing if not peri[m] and m not in EXCLUDED])
print(f'TRUE rule: {len(passing)} of {AD_all.shape[1]} | after dropping '
      f'{int(peri.sum())} pericentromeric and {len(EXCLUDED)} excluded: {len(true_muts)}')

AD = AD_all[true_muts]
DP = DP_all[true_muts]
A = AD.values.astype(float)
D = DP.values.astype(float)
AF = np.divide(A, D, out=np.zeros_like(A), where=D > 0)


##


# 2a. Per-site background: placenta reads shrunk to the substitution-class rate
placenta = (
    full.query('tissue == "placenta" and mutation_id in @true_muts')
    .groupby('mutation_id').agg(AD_alt=('AD_alt', 'sum'), DP=('DP', 'sum'))
    .reindex(true_muts).fillna(0)
)
heart_af = A.sum(0) / D.sum(0)
placenta_af = (placenta['AD_alt'] / placenta['DP'].replace(0, np.nan)).fillna(0)
carries = (
    (binom.sf(placenta['AD_alt'] - 1, placenta['DP'], ERR_INIT) < P_PRESENT)
    & (placenta_af.values > PLACENTA_AF_FRACTION * heart_af)
)
substitution = muts['SBS6'].reindex(true_muts).fillna('NA')

# Two-level shrinkage. C>A and C>G have zero placenta alt reads over ~1,090
# bases each, so their rate is unmeasured rather than low; a pseudo-count alone
# would leave them near zero and make a single read significant. Each class rate
# is therefore shrunk toward the rate pooled over all classes, and each site
# toward its class.
grouped = placenta.loc[~carries].groupby(substitution[~carries])
class_alt, class_dp = grouped['AD_alt'].sum(), grouped['DP'].sum()
pooled_rate = class_alt.sum() / class_dp.sum()
class_rate = (class_alt + CLASS_PRIOR_BASES * pooled_rate) / (class_dp + CLASS_PRIOR_BASES)
class_n = substitution[~carries].value_counts()
print(f'\nPlacenta background pooled over classes: {pooled_rate:.2e} '
      f'(classes with no observed alt read: {", ".join(class_alt.index[class_alt == 0]) or "none"})')
prior = substitution.map(class_rate).fillna(ERR_INIT).values
site_bg = np.where(
    carries, prior,
    (placenta['AD_alt'].values + PRIOR_BASES * prior) / (placenta['DP'].values + PRIOR_BASES)
).clip(min=1e-5)

print(f'\nBackground: {int(carries.sum())} sites placenta-carrying (class rate used); '
      f'median {np.median(site_bg):.1e}, max {site_bg.max():.1e}')
print('  class rates (n sites): ' + ', '.join(
    f'{k} {class_rate[k]:.1e} (n={class_n.get(k, 0)})' for k in class_rate.index))

# Overdispersion per class, fitted on the placenta counts of non-carrying sites
concentration = {}
for cls, idx in substitution[~carries].groupby(substitution[~carries]).groups.items():
    sub = placenta.loc[idx]
    if len(sub) >= 5 and sub['DP'].sum() > 0:
        concentration[cls] = fit_concentration(sub['AD_alt'], sub['DP'], class_rate[cls])
    else:
        concentration[cls] = np.inf     # too few sites: stay binomial
fitted = {k: v for k, v in concentration.items() if np.isfinite(v)}
print('  beta-binomial concentration: ' + ', '.join(
    f'{k} {"binomial (too few sites)" if np.isinf(v) else f"{v:.2e}"}'
    for k, v in concentration.items()))
if fitted and min(fitted.values()) > 1e6:
    print('  no overdispersion detectable in the placenta counts (they are mostly 0-2 reads '
          'per site), so the beta-binomial is numerically identical to the binomial here')

site_s = substitution.map(concentration).fillna(np.inf).values


##


# 2b. Genotypes. Present is one rule; absence is reported at three stringencies.
p_site = np.empty_like(A)
finite = np.isfinite(site_s)
if finite.any():
    cols = np.where(finite)[0]
    p_site[:, cols] = tail_p(A[:, cols], D[:, cols],
                             site_bg[cols][None, :], site_s[cols][None, :])
if (~finite).any():
    cols = np.where(~finite)[0]
    p_site[:, cols] = binom.sf(A[:, cols] - 1, D[:, cols], site_bg[cols][None, :])

present = (A >= MIN_AD_PRESENT) & (p_site < P_PRESENT) & (D >= MIN_DP_PRESENT)
weak = A == 1
lod = detection_limit(D)

states = {}
for af0 in AF_EXCLUDE:
    st = np.full(A.shape, 'undetermined', dtype=object)
    st[present] = 'present'
    st[weak & ~present] = 'weak'
    st[(A == 0) & (lod <= af0)] = 'absent'
    states[af0] = st

single_read_pass = int(((A == 1) & (p_site < P_PRESENT) & (D >= MIN_DP_PRESENT)).sum())
print(f'\nGenotypes ({A.shape[0]} cuts x {A.shape[1]} mutations = {A.size} cells)')
print(f'  present {int(present.sum())} (p<{P_PRESENT}, AD>={MIN_AD_PRESENT}, DP>={MIN_DP_PRESENT}) | '
      f'weak {int(weak.sum())}')
print(f'  the AD>={MIN_AD_PRESENT} rule is redundant: {single_read_pass} single-read cells clear '
      f'the p-value at DP>={MIN_DP_PRESENT}')
for af0 in AF_EXCLUDE:
    st = states[af0]
    print(f'  exclude clones >{af0:.0%}: absent {int((st == "absent").sum()):5d} | '
          f'undetermined {int((st == "undetermined").sum()):5d}')
print(f'  mutations with >=1 present call: {int((present.sum(0) > 0).sum())} of {len(true_muts)}')

main_state = states[AF_MAIN]
determinate = np.isin(main_state, ['present', 'absent'])

# Sensitivity: the iterative all-sample background, rejected during planning
bg_iter = np.full(len(true_muts), np.median(site_bg))
for _ in range(10):
    pres_i = (A >= MIN_AD_PRESENT) & (binom.sf(A - 1, D, bg_iter[None, :]) < P_PRESENT)
    num = np.where(pres_i, 0, A).sum(0) + placenta['AD_alt'].values
    den = np.where(pres_i, 0, D).sum(0) + placenta['DP'].values
    new = np.clip(num / np.maximum(den, 1), site_bg, None)
    if np.allclose(new, bg_iter, rtol=1e-3):
        break
    bg_iter = new
lost = int((pres_i.sum(0) == 0).sum())
print(f'  sensitivity: iterative background would leave {lost} mutations with no present call')


##


# 3. 1a Pre-gastrulation vs heart-specific
# A mutation shared with gut (endoderm) as well as heart (mesoderm) arose before the
# germ layers split. Blood-only sharing is not enough: blood and heart are both
# mesoderm. Placenta reads are not used: at ~210x the number of sites with >=2
# placenta reads is what background noise alone gives (11 vs 9.3 expected).
CROSS_LAYER = ['Blood_Gut', 'Gut']
tree = muts['desc_samples_orgin'].reindex(true_muts)
lineage_class = np.select(
    [tree.isin(CROSS_LAYER).values, tree.eq('Unassigned').values],
    ['Pre-gastrulation', 'Heart-specific'], default='Other shared')
pre_gastrulation = lineage_class == 'Pre-gastrulation'
heart_specific = lineage_class == 'Heart-specific'
print('\n1a Lineage class (tree assignment): ' + ' | '.join(
    f'{k} {int((lineage_class == k).sum())}'
    for k in ['Pre-gastrulation', 'Heart-specific', 'Other shared']))
n_present = present.sum(0)
print('  cuts with a present call: pre-gastrulation median %.0f | heart-specific median %.0f' % (
    np.median(n_present[pre_gastrulation]), np.median(n_present[heart_specific])))


##


# 4. 1b Region coverage and enrichment
reg = region.values
in_region = {r: reg == r for r in REGIONS}

covered = np.array([
    all((present[in_region[r], j]).any() for r in REGIONS) for j in range(len(true_muts))
])
covered_reads = np.array([
    all((A[in_region[r], j] > 0).any() for r in REGIONS) for j in range(len(true_muts))
])
print(f'\n1b Detected in all five regions: {int(covered.sum())} by genotype | '
      f'{int(covered_reads.sum())} by >=1 raw read')

global_af = ((A.sum(0) + placenta['AD_alt'].values)
             / (D.sum(0) + placenta['DP'].values))

# Region vs rest of the heart, on read counts (hard genotype calls are not used).
# For every SNV x region:
#   - read-level binomial: region AD/DP against the rest-of-heart VAF (floored at the
#     site background). Descriptive: every read counts as a replicate.
#   - beta-binomial likelihood-ratio test with tissue blocks as the units (primary), and
#     with cuts as the units (sensitivity). Null: one VAF for the whole heart; alternative:
#     the region's VAF differs from the rest; one-sided via the signed root of the LR.
#     VAFs are bounded below by the site background. Dispersion is not estimated per test:
#     per-SNV estimates are shrunk to a trend on abundance, fitted across SNVs.
#   - enriched (a hit): read-level q < ENRICH_Q, and alternate reads in >= MIN_FRAC_CUTS_ALT of
#     the region's cuts, so a single-cut spike does not count. A SNV is enriched with >=1 hit.
#     Global significance: the whole rule is rerun on N_PERM_ENRICH shuffles of region labels
#     across cuts; statistic = number of enriched SNVs (global p and empirical FDR).
#   - replicated (reported, not used for the call): alternate reads in >=2 region blocks and
#     block-level p < 0.05 after dropping each region block in turn.
LOG_S_BOUNDS = (np.log(1.0), np.log(1e5))
MU_MAX = 0.5
ENRICH_Q = 0.1
MIN_FRAC_CUTS_ALT = 0.3
N_PERM_ENRICH = 1000
N_SIM_NULL = 20
reg = region.values
chunk_of_cut = chunk.values
blocks = np.array(sorted(set(chunk_of_cut)))
block_region = pd.Series(reg, index=chunk_of_cut).groupby(level=0).first().reindex(blocks).values
A_blk = np.vstack([A[chunk_of_cut == b].sum(0) for b in blocks])
D_blk = np.vstack([D[chunk_of_cut == b].sum(0) for b in blocks])


def fit_mu(a, d, s, floor):
    """
    MLE of the beta-binomial mean for counts a/d at fixed concentration s.
    Returns (mu, log-likelihood).
    """
    if len(a) == 0:
        return np.nan, 0.0
    f = lambda mu: -betabinom.logpmf(a, d, mu * s, (1 - mu) * s).sum()
    res = minimize_scalar(f, bounds=(floor, MU_MAX), method='bounded')
    return float(res.x), -float(res.fun)


def fit_log_s(a, d, floor):
    """
    Per-SNV concentration under the single-mean (null) model, profiled over the mean.
    """
    f = lambda ls: -fit_mu(a, d, np.exp(ls), floor)[1]
    return float(minimize_scalar(f, bounds=LOG_S_BOUNDS, method='bounded').x)


def shrunk_concentration(Au, Du):
    """
    Per-SNV log concentration regressed on logit(pooled VAF) across SNVs (estimates at the
    bounds excluded from the fit); each SNV takes its trend value.
    """
    ls = np.array([fit_log_s(Au[:, j], Du[:, j], site_bg[j]) for j in range(Au.shape[1])])
    mu = np.clip(Au.sum(0) / Du.sum(0), 1e-4, None)
    x = np.log(mu / (1 - mu))
    ok = (ls > LOG_S_BOUNDS[0] + .05) & (ls < LOG_S_BOUNDS[1] - .05)
    slope, icpt = np.polyfit(x[ok], ls[ok], 1)
    return np.exp(np.clip(icpt + slope * x, *LOG_S_BOUNDS)), ls, (icpt, slope, int(ok.sum()))


def bb_region_test(a, d, in_r, s, floor):
    """
    One-sided beta-binomial LRT, region (in_r) vs rest, at fixed concentration s.
    """
    _, ll0 = fit_mu(a, d, s, floor)
    mu_r, ll_r = fit_mu(a[in_r], d[in_r], s, floor)
    mu_o, ll_o = fit_mu(a[~in_r], d[~in_r], s, floor)
    stat = max(2 * (ll_r + ll_o - ll0), 0.0)
    return norm.sf(np.sign(mu_r - mu_o) * np.sqrt(stat))


s_blk, ls_blk_raw, trend_blk = shrunk_concentration(A_blk, D_blk)
s_cut, ls_cut_raw, trend_cut = shrunk_concentration(A, D)
print(f'  block-level concentration trend: log s = {trend_blk[0]:.2f} + {trend_blk[1]:.2f} logit(VAF) '
      f'({trend_blk[2]} SNVs off the bounds); median s {np.median(s_blk):.0f}')
print(f'  cut-level concentration trend:   log s = {trend_cut[0]:.2f} + {trend_cut[1]:.2f} logit(VAF) '
      f'({trend_cut[2]} SNVs off the bounds); median s {np.median(s_cut):.0f}')

rng = np.random.default_rng(0)
rows = []
for j, m in enumerate(true_muts):
    for r in REGIONS:
        rb = block_region == r
        rc = reg == r
        a_r, d_r = A[rc, j].sum(), D[rc, j].sum()
        a_o, d_o = A[~rc, j].sum(), D[~rc, j].sum()
        af_o = a_o / d_o
        blk_af = A_blk[rb, j] / D_blk[rb, j]
        # leave-one-block-out, block level
        lobo = []
        if rb.sum() >= 2:
            for b in np.where(rb)[0]:
                keep = np.arange(len(blocks)) != b
                lobo.append(bb_region_test(A_blk[keep, j], D_blk[keep, j], rb[keep], s_blk[j], site_bg[j]))
        rows.append(dict(
            mutation_id=m, region=r,
            AD=int(a_r), DP=int(d_r), AF=a_r / d_r,
            AD_rest=int(a_o), DP_rest=int(d_o), AF_rest=af_o,
            AF_diff=a_r / d_r - af_o,
            block_mean_AF=blk_af.mean(), n_blocks=int(rb.sum()),
            n_blocks_alt=int((A_blk[rb, j] > 0).sum()),
            n_cuts=int(rc.sum()), n_cuts_alt=int((A[rc, j] > 0).sum()),
            frac_cuts_alt=(A[rc, j] > 0).mean(),
            top_block_share=A_blk[rb, j].max() / a_r if a_r else np.nan,
            p_binom=binom.sf(a_r - 1, d_r, max(af_o, site_bg[j])),
            p_bb_block=bb_region_test(A_blk[:, j], D_blk[:, j], rb, s_blk[j], site_bg[j]),
            p_bb_cut=bb_region_test(A[:, j], D[:, j], rc, s_cut[j], site_bg[j]),
            p_lobo_max=max(lobo) if lobo else np.nan,
        ))
enrichment = pd.DataFrame(rows)
for col in ['binom', 'bb_block', 'bb_cut']:
    enrichment[f'q_{col}'] = multipletests(enrichment[f'p_{col}'], method='fdr_bh')[1]
enrichment['replicated'] = (enrichment['n_blocks_alt'] >= 2) & (enrichment['p_lobo_max'] < 0.05)
enrichment['enriched'] = (enrichment['q_binom'] < ENRICH_Q) & \
    (enrichment['frac_cuts_alt'] >= MIN_FRAC_CUTS_ALT)


def enrichment_calls(labels):
    """
    The hit rule on a labelling of the cuts: regions x SNVs boolean matrix.
    """
    p, frac = [], []
    for r in REGIONS:
        k = labels == r
        a, d = A[k].sum(0), D[k].sum(0)
        rest = np.maximum((A.sum(0) - a) / (D.sum(0) - d), site_bg)
        p.append(binom.sf(a - 1, d, rest))
        frac.append((A[k] > 0).mean(0))
    q = multipletests(np.ravel(p), method='fdr_bh')[1].reshape(len(REGIONS), -1)
    return (q < ENRICH_Q) & (np.array(frac) >= MIN_FRAC_CUTS_ALT)


calls = enrichment_calls(reg)
assert calls.sum() == enrichment['enriched'].sum()
n_obs = int(calls.any(0).sum())
perm_rng = np.random.default_rng(1)
n_null = np.array([enrichment_calls(perm_rng.permutation(reg)).any(0).sum()
                   for _ in range(N_PERM_ENRICH)])
enrich_perm = pd.Series({
    'n_enriched_snvs': n_obs, 'n_hits': int(calls.sum()), 'n_tests': calls.size,
    'n_permutations': N_PERM_ENRICH, 'expected_by_chance': n_null.mean(),
    'null_95th': np.percentile(n_null, 95),
    'p_global': (np.sum(n_null >= n_obs) + 1) / (N_PERM_ENRICH + 1),
    'empirical_FDR': n_null.mean() / max(n_obs, 1),
})
enrich_perm.to_csv(os.path.join(path_results, 'ENRICHMENT_PERMUTATION.tsv'), sep='\t', header=False)

hits = enrichment.query('enriched')
print(f'  {len(enrichment)} SNV x region tests | q<{ENRICH_Q}: read-level binomial '
      f'{int((enrichment.q_binom < ENRICH_Q).sum())}, beta-binomial cut level '
      f'{int((enrichment.q_bb_cut < ENRICH_Q).sum())}, block level '
      f'{int((enrichment.q_bb_block < ENRICH_Q).sum())}')
print(f'  hits (read-level q<{ENRICH_Q}, alt reads in >={MIN_FRAC_CUTS_ALT:.0%} of region cuts): '
      f'{len(hits)} | enriched SNVs {n_obs} | by region '
      f'{hits.region.value_counts().reindex(REGIONS).fillna(0).astype(int).to_dict()}')
print(f'  cut-label permutation ({N_PERM_ENRICH}): {n_obs} enriched SNVs vs {n_null.mean():.1f} expected '
      f'(95th pct {np.percentile(n_null, 95):.0f}) | global p = {enrich_perm.p_global:.4f} | '
      f'empirical FDR {enrich_perm.empirical_FDR:.2f}')

# Calibration of the block-level test: data simulated under each SNV's fitted null
sim_p = []
for j in range(len(true_muts)):
    mu0, _ = fit_mu(A_blk[:, j], D_blk[:, j], s_blk[j], site_bg[j])
    for _ in range(N_SIM_NULL):
        a_sim = betabinom.rvs(D_blk[:, j].astype(int), mu0 * s_blk[j], (1 - mu0) * s_blk[j], random_state=rng)
        for r in REGIONS:
            sim_p.append(bb_region_test(a_sim, D_blk[:, j], block_region == r, s_blk[j], site_bg[j]))
sim_p = np.array(sim_p)
print(f'  null calibration, block level ({N_SIM_NULL} simulations per SNV, discovery filter not '
      f'applied): p<0.05 {np.mean(sim_p < .05):.3f}, p<0.01 {np.mean(sim_p < .01):.3f}')

enrichment.to_csv(os.path.join(path_results, 'REGION_ENRICHMENT.tsv'), sep='\t', index=False)


##


# 5. 1c Diffuse vs scattered
chunks = chunk.values
uniq_chunks = sorted(set(chunks))
per_chunk = np.array([
    [present[(chunks == c), j].sum() for c in uniq_chunks] for j in range(len(true_muts))
])
k_calls = present.sum(0)
max_in_chunk = per_chunk.max(1)
n_chunks_hit = (per_chunk > 0).sum(1)

klass = np.where(k_calls == 0, 'no call',
                 np.where(max_in_chunk >= 2, 'diffuse',
                          np.where(n_chunks_hit >= 2, 'scattered', 'single cut')))

# Chance baseline: shuffle each mutation's k present calls across cuts
chunk_index = pd.Series(np.arange(len(uniq_chunks)), index=uniq_chunks).loc[chunks].values
exp_diffuse = np.zeros(len(true_muts))
for j in range(len(true_muts)):
    k = k_calls[j]
    if k < 2:
        continue
    hits = 0
    for _ in range(1000):
        pick = rng.choice(len(chunks), k, replace=False)
        if np.bincount(chunk_index[pick], minlength=len(uniq_chunks)).max() >= 2:
            hits += 1
    exp_diffuse[j] = hits / 1000

print('\n1c Spatial pattern of present calls')
print(pd.Series(klass).value_counts().to_string())
obs_d = (klass == 'diffuse')
print(f'  diffuse observed {int(obs_d.sum())} of {int((k_calls >= 2).sum())} testable; '
      f'expected by chance {exp_diffuse[k_calls >= 2].sum():.1f}')
strict = obs_d & (exp_diffuse < 0.5)
print(f'  diffuse beyond chance (expected probability < 0.5): {int(strict.sum())}')


##


# 6. 1d Local sweep
mean_af = np.array([AF[present[:, j], j].mean() if present[:, j].any() else np.nan
                    for j in range(len(true_muts))])
max_af = np.array([AF[present[:, j], j].max() if present[:, j].any() else np.nan
                   for j in range(len(true_muts))])

# Mean AF over present cuts is truncated from below by the calling rule: a cut
# can only be present at AF >= 2/DP. Pooling reads over the whole chunk that
# holds the present calls is not conditioned on per-cut detection.
chunk_af = np.full((len(uniq_chunks), len(true_muts)), np.nan)
for i, c in enumerate(uniq_chunks):
    rows = chunks == c
    dp = D[rows].sum(0)
    chunk_af[i] = np.divide(A[rows].sum(0), dp, out=np.full(len(true_muts), np.nan), where=dp > 0)
local_af = np.array([chunk_af[per_chunk[j] > 0, j].max() if (per_chunk[j] > 0).any() else np.nan
                     for j in range(len(true_muts))])

floor = np.where(present, 2 / np.maximum(D, 1), np.nan)
print('\n1d Local sweep (carrier fraction = %dx AF)' % PLOIDY)
print('  mean AF over present cuts: median %.3f (IQR %.3f-%.3f) -> carrier cells %.0f%%' % (
    np.nanmedian(mean_af), np.nanpercentile(mean_af, 25), np.nanpercentile(mean_af, 75),
    100 * PLOIDY * np.nanmedian(mean_af)))
print('  but %.0f%% of present cells sit exactly at the 2-read floor (median floor AF %.3f), '
      'so this mostly measures depth' % (100 * (A[present] == 2).mean(), np.nanmedian(floor)))
print('  chunk-pooled AF where present (not conditioned on per-cut detection):')
print('    median %.3f (IQR %.3f-%.3f) -> carrier cells %.0f%%' % (
    np.nanmedian(local_af), np.nanpercentile(local_af, 25), np.nanpercentile(local_af, 75),
    100 * PLOIDY * np.nanmedian(local_af)))
for name, sel in [('pre-gastrulation', pre_gastrulation), ('heart-specific', heart_specific),
                  ('diffuse', klass == 'diffuse'), ('scattered', klass == 'scattered')]:
    if sel.sum():
        print('    %-14s n=%3d | chunk-pooled AF median %.3f | max cut AF median %.3f' % (
            name, int(sel.sum()), np.nanmedian(local_af[sel]), np.nanmedian(max_af[sel])))


##


# Per-mutation summary table
summary = pd.DataFrame({
    'mutation_id': true_muts,
    'substitution': substitution.values,
    'n_cuts_3reads': strong[true_muts].sum().values,
    'n_present_calls': k_calls,
    'n_chunks_with_present': n_chunks_hit,
    'max_present_in_a_chunk': max_in_chunk,
    'spatial_class': klass,
    'p_diffuse_by_chance': exp_diffuse.round(3),
    'placenta_AD': placenta['AD_alt'].values.astype(int),
    'placenta_DP': placenta['DP'].values.astype(int),
    'tree_assignment': muts['desc_samples_orgin'].reindex(true_muts).values,
    'lineage_class': lineage_class,
    'site_background': site_bg,
    'global_AF': global_af,
    'mean_AF_present': mean_af,
    'max_AF_present': max_af,
    'detected_all_regions': covered,
}).set_index('mutation_id')
summary.to_csv(os.path.join(path_results, 'LINEAGE_SUMMARY.tsv'), sep='\t')

summary.reset_index()[['mutation_id', 'n_cuts_3reads', 'placenta_AD', 'placenta_DP',
                       'tree_assignment', 'lineage_class', 'site_background']].to_csv(
    os.path.join(path_results, 'TRUE_MUTATIONS.tsv'), sep='\t', index=False)

genotypes = pd.DataFrame({
    'Sample_ID': np.repeat(AD.index.values, len(true_muts)),
    'mutation_id': np.tile(true_muts.values, len(AD.index)),
    'AD_alt': A.ravel().astype(int), 'DP': D.ravel().astype(int),
    'AF': AF.ravel().round(4),
    'detection_limit_AF': lod.ravel().round(3),
    'site_background': np.tile(site_bg, len(AD.index)),
    'p_presence': p_site.ravel(),
}).assign(
    region=lambda x: x['Sample_ID'].map(samples['region']),
    chunk=lambda x: x['Sample_ID'].map(samples['chunk']),
)
for af0 in AF_EXCLUDE:
    genotypes[f'state_af{int(af0 * 100)}'] = states[af0].ravel()
genotypes.to_csv(os.path.join(path_results, 'GENOTYPES_TRUE.tsv.gz'), sep='\t', index=False)
pd.DataFrame(AF, index=AD.index, columns=true_muts).to_csv(
    os.path.join(path_results, 'AF_MATRIX.tsv.gz'), sep='\t')


##


# 7. 2a Region-level hierarchy of lineage composition, from read counts.
# Regional VAF profiles are built three ways: pooled AD/DP over the region (primary),
# equal-weight average of block VAFs, and mean of per-cut VAFs. Each is clustered by cosine
# distance and average linkage, for all SNVs and for each lineage class. Clade support comes
# from two resamplings: SNVs (80% without replacement) and whole blocks (with replacement
# within each region; CS has one block, so it never varies). A clade enters the consensus
# tree only with block support >= CONSENSUS_MIN; weaker branches are left unresolved.
CONSENSUS_MIN = 0.5
PROFILE_KINDS = ['pooled', 'equal_block', 'mean_cut']
blocks_of = {r: np.where(block_region == r)[0] for r in REGIONS}
cuts_of_block = [np.where(chunk_of_cut == b)[0] for b in blocks]


def build_profiles(kind, cols, blk=None):
    """
    Region x SNV VAF profile. blk maps each region to the block indices to use (with
    repeats under the block bootstrap); default all its blocks.
    """
    blk = blk or blocks_of
    out = []
    for r in REGIONS:
        b = blk[r]
        if kind == 'pooled':
            out.append(A_blk[np.ix_(b, cols)].sum(0) / D_blk[np.ix_(b, cols)].sum(0))
        elif kind == 'equal_block':
            out.append((A_blk[np.ix_(b, cols)] / D_blk[np.ix_(b, cols)]).mean(0))
        else:
            c = np.concatenate([cuts_of_block[k] for k in b])
            out.append(AF[np.ix_(c, cols)].mean(0))
    return np.vstack(out)


def tree_of(P):
    D_ = pairwise_distances(P, metric='cosine')
    return D_, linkage(squareform(D_, checks=False), method='average')


row_idx = {r: np.where(in_region[r])[0] for r in REGIONS}
snv_sets = {'All': np.arange(len(true_muts))}
snv_sets.update({k: np.where(lineage_class == k)[0]
                 for k in ['Pre-gastrulation', 'Heart-specific', 'Other shared']})

print('\n2a Region-level hierarchy (cosine, average linkage)')
tree_rows, dist_rows = [], []
for kind in PROFILE_KINDS:
    for set_name, cols in snv_sets.items():
        D_, Z = tree_of(build_profiles(kind, cols))
        observed = clade_set(Z, REGIONS)
        n_char = int(CHAR_FRACTION * len(cols))
        hits_snv, hits_blk = {}, {}
        for _ in range(N_BOOT):
            for grp, _h in clade_set(tree_of(build_profiles(
                    kind, rng.choice(cols, n_char, replace=False)))[1], REGIONS):
                hits_snv[grp] = hits_snv.get(grp, 0) + 1
            blk = {r: rng.choice(b, len(b), replace=True) for r, b in blocks_of.items()}
            for grp, _h in clade_set(tree_of(build_profiles(kind, cols, blk))[1], REGIONS):
                hits_blk[grp] = hits_blk.get(grp, 0) + 1
        for grp in set(dict(observed)) | {g for g, v in hits_blk.items() if v / N_BOOT >= .2}:
            tree_rows.append(dict(
                profile=kind, snv_set=set_name, n_snvs=len(cols), clade='+'.join(sorted(grp)),
                in_tree=grp in dict(observed), height=dict(observed).get(grp, np.nan),
                support_snv=hits_snv.get(grp, 0) / N_BOOT, support_block=hits_blk.get(grp, 0) / N_BOOT))
        for a, b in itertools.combinations(range(len(REGIONS)), 2):
            dist_rows.append(dict(profile=kind, snv_set=set_name, region_a=REGIONS[a],
                                  region_b=REGIONS[b], cosine_distance=D_[a, b]))
hierarchy = pd.DataFrame(tree_rows)
hierarchy['consensus'] = hierarchy['in_tree'] & (hierarchy['support_block'] >= CONSENSUS_MIN)
hierarchy.to_csv(os.path.join(path_results, 'HIERARCHY.tsv'), sep='\t', index=False)
pd.DataFrame(dist_rows).to_csv(os.path.join(path_results, 'REGION_DISTANCES.tsv'), sep='\t', index=False)

show = hierarchy.query('in_tree').copy()
show['support'] = (100 * show.support_snv).round().astype(int).astype(str) + ' / ' + \
    (100 * show.support_block).round().astype(int).astype(str)
print('  clades in each tree, support % (SNV jackknife / block bootstrap); '
      f'consensus needs block support >= {CONSENSUS_MIN:.0%}:')
print(show.pivot_table(index='clade', columns=['profile', 'snv_set'], values='support', aggfunc='first')
      .reindex(columns=[(k, s_) for k in PROFILE_KINDS for s_ in snv_sets]).fillna('')
      .to_string().replace('\n', '\n    '))
prim = hierarchy.query('profile == "pooled" and snv_set == "All" and in_tree')
print('  consensus (pooled, all SNVs): ' + (', '.join(prim.query('consensus').clade) or 'none resolved'))


##


# 8. 2b Region-pair tests at cut level, paired design
sim_af = 1 - pairwise_distances(AF, metric='cosine')
np.fill_diagonal(sim_af, np.nan)
pres_det = np.where(determinate, present, np.nan)
sim_gt = np.full((len(AD), len(AD)), np.nan)
for i in range(len(AD)):
    for j in range(len(AD)):
        if i == j:
            continue
        both = determinate[i] & determinate[j]
        if both.sum() == 0:
            continue
        pi, pj = present[i, both], present[j, both]
        union = (pi | pj).sum()
        sim_gt[i, j] = (pi & pj).sum() / union if union else np.nan

measures = {'AF cosine': sim_af, 'genotype Jaccard': sim_gt}


def normalised(S):
    """
    Divide each column by that cut's mean similarity to all others.
    """
    level = np.nanmean(S, axis=1)
    return S / level[None, :]


def mean_to(S, source, targets):
    sel = targets.copy()
    sel[source] = False
    return np.nanmean(S[source, sel])


print('\n2b Region pairs, cut level')
for name, S in measures.items():
    X = normalised(S)
    print(f'\n  {name}')
    obs_pairs = {}
    for a, b in itertools.combinations_with_replacement(REGIONS, 2):
        ia, ib = row_idx[a], row_idx[b]
        sub = S[np.ix_(ia, ib)]
        obs_pairs[f'{a}-{b}'] = np.nanmean(sub[~np.eye(len(ia), dtype=bool)] if a == b else sub)
    null = {k: [] for k in obs_pairs}
    for _ in range(N_PERM):
        lab = rng.permutation(reg)
        idx = {r: np.where(lab == r)[0] for r in REGIONS}
        for a, b in itertools.combinations_with_replacement(REGIONS, 2):
            sub = S[np.ix_(idx[a], idx[b])]
            null[f'{a}-{b}'].append(
                np.nanmean(sub[~np.eye(len(idx[a]), dtype=bool)] if a == b else sub))
    tab = pd.DataFrame({
        'similarity': pd.Series(obs_pairs).round(4),
        'null_mean': pd.Series({k: np.mean(v) for k, v in null.items()}).round(4),
        'p': pd.Series({k: perm_pvalue(null[k], obs_pairs[k]) for k in obs_pairs}).round(3),
    }).sort_values('similarity', ascending=False)
    print(tab.to_string())

    # Paired contrasts, partner-normalised
    contrasts = []
    for label, fixed, left, right in [
        ('LS-LV vs LS-RV', 'LS', 'LV', 'RV'),
        ('RS-RV vs RS-LV', 'RS', 'RV', 'LV'),
        ('septum-LV vs septum-RV', 'septum', 'LV', 'RV'),
    ]:
        rows_fixed = (np.where(np.isin(reg, ['LS', 'CS', 'RS']))[0] if fixed == 'septum'
                      else row_idx[fixed])
        l = np.array([mean_to(X, i, in_region[left]) for i in rows_fixed])
        r_ = np.array([mean_to(X, i, in_region[right]) for i in rows_fixed])
        ok = ~(np.isnan(l) | np.isnan(r_))
        p = wilcoxon(l[ok], r_[ok])[1] if ok.sum() >= 6 else np.nan
        contrasts.append(dict(contrast=label, n=int(ok.sum()),
                              left=round(np.nanmean(l), 4), right=round(np.nanmean(r_), 4), p=p))

    # Side matching: the two lateral contrasts above ask the same question of
    # opposite sides, so pool their per-cut differences into one test.
    d_ls = np.array([mean_to(X, i, in_region['LV']) - mean_to(X, i, in_region['RV'])
                     for i in row_idx['LS']])
    d_rs = np.array([mean_to(X, i, in_region['RV']) - mean_to(X, i, in_region['LV'])
                     for i in row_idx['RS']])
    pooled = np.r_[d_ls, d_rs]
    pooled = pooled[~np.isnan(pooled)]
    contrasts.append(dict(
        contrast='side matching (LS+RS pooled)', n=len(pooled),
        left=round(float((pooled > 0).sum()), 4), right=round(float((pooled < 0).sum()), 4),
        p=wilcoxon(pooled)[1]))

    # CS-LS vs CS-RS: the compared cuts (LS vs RS) sit on the free side, so each
    # cut's similarity to CS is divided by its similarity to everything else.
    pref = np.array([mean_to(X, i, in_region['CS']) / mean_to(X, i, ~in_region['CS'])
                     if not np.isnan(mean_to(X, i, ~in_region['CS'])) else np.nan
                     for i in range(len(AD))])
    v1 = pref[row_idx['LS']][~np.isnan(pref[row_idx['LS']])]
    v2 = pref[row_idx['RS']][~np.isnan(pref[row_idx['RS']])]
    obs = v1.mean() - v2.mean()
    lat = np.where(np.isin(reg, ['LS', 'RS']))[0]
    nulls = []
    for _ in range(N_PERM):
        perm = rng.permutation(reg[lat])
        nulls.append(np.nanmean(pref[lat][perm == 'LS']) - np.nanmean(pref[lat][perm == 'RS']))
    contrasts.append(dict(contrast='CS-LS vs CS-RS', n=len(v1) + len(v2),
                          left=round(v1.mean(), 4), right=round(v2.mean(), 4),
                          p=perm_pvalue(nulls, obs)))

    # CS-LS vs RS-LS: LS cuts are the fixed side, CS and RS the partners.
    to_cs = np.array([mean_to(X, i, in_region['CS']) for i in row_idx['LS']])
    to_rs = np.array([mean_to(X, i, in_region['RS']) for i in row_idx['LS']])
    ok = ~(np.isnan(to_cs) | np.isnan(to_rs))
    contrasts.append(dict(contrast='CS-LS vs RS-LS', n=int(ok.sum()),
                          left=round(np.nanmean(to_cs), 4), right=round(np.nanmean(to_rs), 4),
                          p=wilcoxon(to_cs[ok], to_rs[ok])[1] if ok.sum() >= 6 else np.nan))

    ct = pd.DataFrame(contrasts)
    ct['q'] = multipletests(ct['p'].fillna(1), method='fdr_bh')[1].round(3)
    print(ct.round(4).to_string(index=False))
    ct.insert(0, 'measure', name)
    ct.to_csv(os.path.join(path_results, 'REGION_CONTRASTS.tsv'), sep='\t',
              mode='a' if name != list(measures)[0] else 'w',
              header=name == list(measures)[0], index=False)

print('\nAll similarities are partner-normalised for per-cut level; cut-level '
      'permutation treats cuts as independent.')


##


# 9. Ventricular lean of the septal regions. Each septal cut is scored by how
# much more it resembles left- than right-ventricular cuts, after the partner
# cuts are normalised by their own overall similarity (section 8). The pattern
# is then tested jointly, and checked against dropping each tissue block.
S_af = 1 - pairwise_distances(AF, metric='cosine')
np.fill_diagonal(S_af, np.nan)
X_af = normalised(S_af)

lean = np.array([
    mean_to(X_af, i, in_region['LV']) - mean_to(X_af, i, in_region['RV'])
    for i in range(len(AD.index))
])
septal = np.isin(region.values, ['LS', 'CS', 'RS'])

print('\n9 Ventricular lean of septal cuts (normalised similarity to LV minus RV)')
for r in ['LS', 'RS', 'CS']:
    v = lean[row_idx[r]]
    print(f'  {r} n={len(v):2d} | mean {v.mean():+.4f} | median {np.median(v):+.4f} | '
          f'leaning LV {int((v > 0).sum())}/{len(v)}')

# Joint test: the three septal labels are shuffled across the septal cuts, and a
# replicate counts only if it reproduces the whole pattern at least as strongly.
observed = [lean[row_idx[r]].mean() for r in ['LS', 'RS', 'CS']]
separation = min(observed[0], observed[1]) - observed[2]
septal_idx = np.where(septal)[0]
hits = 0
for _ in range(N_PERM_JOINT):
    labels = region.values.copy()
    labels[septal_idx] = rng.permutation(region.values[septal_idx])
    means = [lean[labels == r].mean() for r in ['LS', 'RS', 'CS']]
    if (means[0] > 0 and means[1] > 0 and means[2] < 0
            and min(means[0], means[1]) - means[2] >= separation):
        hits += 1
print(f'  joint pattern (LS>0, RS>0, CS<0, separation >= {separation:.4f}): '
      f'{hits} of {N_PERM_JOINT} permutations, p = {hits / N_PERM_JOINT:.4f}')

mw = mannwhitneyu(lean[row_idx['CS']], lean[np.r_[row_idx['LS'], row_idx['RS']]])
print(f'  CS cuts vs lateral septal cuts: Mann-Whitney p = {mw[1]:.4f}')


def lean_fraction(target, drop=None, n=4, draws=300):
    """
    Fraction of draws in which a region's profile is closer to LV than to RV,
    with both ventricle profiles built from n cuts, optionally dropping a chunk.
    """
    keep = {r: [i for i in row_idx[r] if drop is None or chunk.values[i] != drop]
            for r in REGIONS}
    if min(len(keep[target]), len(keep['LV']), len(keep['RV'])) < n:
        return np.nan
    out = []
    for _ in range(draws):
        a = AF[rng.choice(keep[target], n, replace=False)].mean(0)
        lv = AF[rng.choice(keep['LV'], n, replace=False)].mean(0)
        rv = AF[rng.choice(keep['RV'], n, replace=False)].mean(0)
        out.append(pairwise_distances(np.vstack([a, rv]), metric='cosine')[0, 1]
                   - pairwise_distances(np.vstack([a, lv]), metric='cosine')[0, 1])
    return float(np.mean(np.array(out) > 0))


loo = pd.DataFrame(
    {r: {('none' if c is None else c): lean_fraction(r, drop=c)
         for c in [None] + sorted(set(chunk.values))}
     for r in ['LS', 'RS', 'CS']}
)
print('\n  leave-one-chunk-out, fraction of draws with LV closer than RV '
      '(>0.5 means leaning LV):')
print(loo.round(2).to_string().replace('\n', '\n  '))

pd.DataFrame({
    'Sample_ID': AD.index, 'region': region.values, 'chunk': chunk.values,
    'lean_LV_minus_RV': lean.round(5),
}).to_csv(os.path.join(path_results, 'VENTRICULAR_LEAN.tsv'), sep='\t', index=False)
loo.to_csv(os.path.join(path_results, 'VENTRICULAR_LEAN_LOO.tsv'), sep='\t')
np.save(os.path.join(path_results, 'ventricular_lean_null.npy'),
        np.array([hits, N_PERM_JOINT, separation]))
