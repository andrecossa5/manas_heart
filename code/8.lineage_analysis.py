"""
Lineage characterisation in the developing heart, and phylogenetic relationships
between heart regions.

Supersedes 7.genotyping.py for the TRUE mutation definition fixed in PLAN.md:
>=3 alt reads in >=2 LCM cuts, with >=5 alt reads summed over one chunk holding
such a cut; pericentromeric sites dropped.

Outputs (results/): TRUE_MUTATIONS.tsv, GENOTYPES_TRUE.tsv.gz, AF_MATRIX.tsv.gz,
LINEAGE_SUMMARY.tsv, REGION_ENRICHMENT.tsv, REGION_CONTRASTS.tsv,
VENTRICULAR_LEAN.tsv, VENTRICULAR_LEAN_LOO.tsv.
"""

import os
import itertools
import numpy as np
import pandas as pd
from scipy.stats import binom, betabinom, fisher_exact, wilcoxon, mannwhitneyu
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
N_PERM_FISHER = 200     # fewer: each permutation costs 620 Fisher tests
N_BOOT = 1000
CHAR_FRACTION = 0.8     # characters resampled per bootstrap replicate
PLOIDY = 2              # donor is female (chrX depth ratio 1.02), so 2 everywhere

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
true_muts = pd.Index([m for m in passing if not peri[m]])
print(f'TRUE rule: {len(passing)} of {AD_all.shape[1]} | after dropping '
      f'{int(peri.sum())} pericentromeric: {len(true_muts)}')

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


# 3. 1a Pre-existing vs heart-specific
tree_assigned = muts['desc_samples_orgin'].reindex(true_muts).ne('Unassigned')
pre_existing = (placenta['AD_alt'].values >= 2) | tree_assigned.values
print(f'\n1a Pre-existing (>=2 placenta reads or tree-assigned): {int(pre_existing.sum())} '
      f'({100 * pre_existing.mean():.0f}%) | heart-specific: {int((~pre_existing).sum())} '
      f'({100 * (~pre_existing).mean():.0f}%)')
print('  by evidence: placenta only %d | tree only %d | both %d' % (
    int(((placenta['AD_alt'].values >= 2) & ~tree_assigned.values).sum()),
    int(((placenta['AD_alt'].values < 2) & tree_assigned.values).sum()),
    int(((placenta['AD_alt'].values >= 2) & tree_assigned.values).sum())))
n_present = present.sum(0)
print('  cuts with a present call: pre-existing median %.0f | heart-specific median %.0f' % (
    np.median(n_present[pre_existing]), np.median(n_present[~pre_existing])))


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

# Test A: Fisher on determinate genotypes, region vs rest
rows = []
for j, m in enumerate(true_muts):
    for r in REGIONS:
        sel = in_region[r] & determinate[:, j]
        oth = (~in_region[r]) & determinate[:, j]
        a = int((main_state[sel, j] == 'present').sum())
        b = int((main_state[sel, j] == 'absent').sum())
        c = int((main_state[oth, j] == 'present').sum())
        d_ = int((main_state[oth, j] == 'absent').sum())
        if a + b == 0 or c + d_ == 0:
            continue
        rows.append(dict(mutation_id=m, region=r, present=a, absent=b,
                         present_rest=c, absent_rest=d_,
                         p=fisher_exact([[a, b], [c, d_]])[1]))
fisher_tab = pd.DataFrame(rows)
fisher_tab['q'] = np.nan
for m, g in fisher_tab.groupby('mutation_id'):
    fisher_tab.loc[g.index, 'q'] = multipletests(g['p'], method='fdr_bh')[1]
obs_hits = int((fisher_tab['q'] < 0.1).sum())


def count_hits(labels):
    """
    Number of mutation x region Fisher tests at p<0.05 under a labelling.
    """
    n = 0
    for j in range(len(true_muts)):
        det = determinate[:, j]
        for r in REGIONS:
            sel = (labels == r) & det
            oth = (labels != r) & det
            a = int((main_state[sel, j] == 'present').sum())
            b = int((main_state[sel, j] == 'absent').sum())
            c = int((main_state[oth, j] == 'present').sum())
            d_ = int((main_state[oth, j] == 'absent').sum())
            if a + b == 0 or c + d_ == 0:
                continue
            if fisher_exact([[a, b], [c, d_]])[1] < 0.05:
                n += 1
    return n


rng = np.random.default_rng(0)
obs_nominal = int((fisher_tab['p'] < 0.05).sum())
null_hits = np.array([count_hits(rng.permutation(reg)) for _ in range(N_PERM_FISHER)])
print(f'  Fisher (determinate calls only): {obs_hits} region hits at FDR<0.1; '
      f'{obs_nominal} at p<0.05 vs permuted mean {null_hits.mean():.1f} '
      f'(95th pct {np.percentile(null_hits, 95):.0f}), p = {(null_hits >= obs_nominal).mean():.3f}')

# Test B: density against the site's global AF
global_af = ((A.sum(0) + placenta['AD_alt'].values)
             / (D.sum(0) + placenta['DP'].values))
dens_rows = []
for j, m in enumerate(true_muts):
    for r in REGIONS:
        a = A[in_region[r], j].sum()
        d_ = D[in_region[r], j].sum()
        dens_rows.append(dict(mutation_id=m, region=r, AD=int(a), DP=int(d_),
                              AF=a / d_ if d_ else np.nan, global_AF=global_af[j],
                              p=binom.sf(a - 1, d_, global_af[j])))
dens_tab = pd.DataFrame(dens_rows)
dens_tab['q'] = multipletests(dens_tab['p'], method='fdr_bh')[1]
print(f'  Density vs global AF: {int((dens_tab["q"] < 0.1).sum())} mutation x region '
      f'hits at FDR<0.1, over {len(dens_tab)} tests')

enrichment = fisher_tab.merge(
    dens_tab[['mutation_id', 'region', 'AD', 'DP', 'AF', 'global_AF', 'p', 'q']],
    on=['mutation_id', 'region'], suffixes=('_fisher', '_density'), how='outer'
)
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
for name, sel in [('pre-existing', pre_existing), ('heart-specific', ~pre_existing),
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
    'pre_existing': pre_existing,
    'site_background': site_bg,
    'global_AF': global_af,
    'mean_AF_present': mean_af,
    'max_AF_present': max_af,
    'detected_all_regions': covered,
}).set_index('mutation_id')
summary.to_csv(os.path.join(path_results, 'LINEAGE_SUMMARY.tsv'), sep='\t')

summary.reset_index()[['mutation_id', 'n_cuts_3reads', 'placenta_AD', 'placenta_DP',
                       'tree_assignment', 'pre_existing', 'site_background']].to_csv(
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


# 7. 2a Region-level hierarchy, with character and cut bootstraps
def region_profiles(rows, cols):
    """
    Per-region mean of per-cut AF (averaging cuts, not pooling reads).
    """
    out = np.vstack([AF[np.ix_(rows[r], cols)].mean(0) for r in REGIONS])
    return out


def region_presence(rows, cols):
    """
    Region x mutation binary matrix: 1 if at least one present call.
    """
    return np.vstack([present[np.ix_(rows[r], cols)].any(0) for r in REGIONS]).astype(bool)


row_idx = {r: np.where(in_region[r])[0] for r in REGIONS}
all_cols = np.arange(len(true_muts))

D_cos = pairwise_distances(region_profiles(row_idx, all_cols), metric='cosine')
D_jac = pairwise_distances(region_presence(row_idx, all_cols), metric='jaccard')

print('\n2a Region-level structure')
print('  region x region cosine distance on mean AF profiles:')
print(pd.DataFrame(D_cos, index=REGIONS, columns=REGIONS).round(3).to_string()
      .replace('\n', '\n    '))
print('  mean distance to the other four regions: ' + ', '.join(
    f'{r} {np.delete(D_cos[i], i).mean():.3f} ({len(row_idx[r])} cuts)'
    for i, r in enumerate(REGIONS)))

# Regions are fixed; support comes from resampling characters (mutations) at
# CHAR_FRACTION without replacement, the usual jackknife of phylogenetics.
n_char = int(CHAR_FRACTION * len(true_muts))
for name, metric, builder in [('AF cosine', 'cosine', region_profiles),
                              ('Jaccard on region presence', 'jaccard', region_presence)]:
    Z = linkage(squareform(pairwise_distances(builder(row_idx, all_cols), metric=metric),
                           checks=False), method='average')
    observed = clade_set(Z, REGIONS)
    hits = {}
    for _ in range(N_BOOT):
        cols = rng.choice(len(true_muts), n_char, replace=False)
        Zb = linkage(squareform(pairwise_distances(builder(row_idx, cols), metric=metric),
                                checks=False), method='average')
        for grp, _h in clade_set(Zb, REGIONS):
            hits[grp] = hits.get(grp, 0) + 1
    print(f'  {name} ({n_char} of {len(true_muts)} characters per replicate, '
          f'{N_BOOT} replicates):')
    for grp, height in observed:
        print(f'    {"+".join(sorted(grp)):14s} merge height {height:.3f} | '
              f'support {100 * hits.get(grp, 0) / N_BOOT:3.0f}%')
    competing = {'+'.join(sorted(k)): round(100 * v / N_BOOT)
                 for k, v in sorted(hits.items(), key=lambda x: -x[1])
                 if k not in dict(observed)}
    print('    competing groups: ' + str(dict(list(competing.items())[:4])))
print('  regions are held fixed here, so this measures character support only, '
      'not how much the topology depends on which cuts were sampled')


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
