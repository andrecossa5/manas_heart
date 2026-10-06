"""
Lineage inputs of the final figure, from the final callset (results/ALLELIC_TABLE_FINAL.tsv.gz and
results/ALLELIC_TABLE.tsv.gz): TRUE mutation set, per-site background from the placenta, genotypes,
lineage classes, and the read-level region-vs-rest enrichment test.

TRUE mutation: >=3 alt reads in >=2 LCM cuts, with >=5 alt reads summed over one chunk holding
such a cut; pericentromeric sites dropped.

Outputs (results/): LINEAGE_SUMMARY.tsv, GENOTYPES_TRUE.tsv.gz, REGION_ENRICHMENT.tsv.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import binom, betabinom
from scipy.optimize import minimize_scalar
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
AF_MAIN = 0.10          # absent = no alt read at a depth that detects a clone at this AF

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


##


# Paths
path_main = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))     # repository root
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


# 2b. Genotypes. Present is one rule; absent needs the depth to detect a clone at AF_MAIN.
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

state = np.full(A.shape, 'undetermined', dtype=object)
state[present] = 'present'
state[weak & ~present] = 'weak'
state[(A == 0) & (lod <= AF_MAIN)] = 'absent'

print(f'\nGenotypes ({A.shape[0]} cuts x {A.shape[1]} mutations = {A.size} cells)')
print(f'  present {int(present.sum())} (p<{P_PRESENT}, AD>={MIN_AD_PRESENT}, DP>={MIN_DP_PRESENT}) | '
      f'weak {int(weak.sum())} | absent {int((state == "absent").sum())} | '
      f'undetermined {int((state == "undetermined").sum())}')
print(f'  mutations with >=1 present call: {int((present.sum(0) > 0).sum())} of {len(true_muts)}')


##


# 3. Lineage class
# A mutation shared with gut (endoderm) as well as heart (mesoderm) arose before the
# germ layers split. Blood-only sharing is not enough: blood and heart are both
# mesoderm. Placenta reads are not used: at ~210x the number of sites with >=2
# placenta reads is what background noise alone gives (11 vs 9.3 expected).
CROSS_LAYER = ['Blood_Gut', 'Gut']
tree = muts['desc_samples_orgin'].reindex(true_muts)
lineage_class = np.select(
    [tree.isin(CROSS_LAYER).values, tree.eq('Unassigned').values],
    ['Pre-gastrulation', 'Heart-specific'], default='Other shared')
print('\nLineage class (tree assignment): ' + ' | '.join(
    f'{k} {int((lineage_class == k).sum())}'
    for k in ['Pre-gastrulation', 'Heart-specific', 'Other shared']))


##


# 4. Region vs rest of the heart, on read counts (hard genotype calls are not used). For every
# SNV x region, a read-level one-sided binomial: region AD/DP against the rest-of-heart VAF
# (floored at the site background); BH over all tests. Every read counts as a replicate, so the
# test is descriptive: the figure scripts take q < 0.1 as enriched and test the number of
# enriched SNVs against sample-label permutations.
reg = region.values
rows = []
for j, m in enumerate(true_muts):
    for r in REGIONS:
        rc = reg == r
        a_r, d_r = A[rc, j].sum(), D[rc, j].sum()
        a_o, d_o = A[~rc, j].sum(), D[~rc, j].sum()
        af_o = a_o / d_o
        rows.append(dict(
            mutation_id=m, region=r,
            AD=int(a_r), DP=int(d_r), AF=a_r / d_r,
            AD_rest=int(a_o), DP_rest=int(d_o), AF_rest=af_o,
            AF_diff=a_r / d_r - af_o,
            n_cuts=int(rc.sum()), n_cuts_alt=int((A[rc, j] > 0).sum()),
            frac_cuts_alt=(A[rc, j] > 0).mean(),
            p_binom=binom.sf(a_r - 1, d_r, max(af_o, site_bg[j])),
        ))
enrichment = pd.DataFrame(rows)
enrichment['q_binom'] = multipletests(enrichment['p_binom'], method='fdr_bh')[1]
enrichment.to_csv(os.path.join(path_results, 'REGION_ENRICHMENT.tsv'), sep='\t', index=False)
print(f'\n{len(enrichment)} SNV x region tests | read-level binomial q<0.1: {int((enrichment.q_binom < .1).sum())}')


##


# Per-mutation summary and genotype tables
summary = pd.DataFrame({
    'mutation_id': true_muts,
    'substitution': substitution.values,
    'n_cuts_3reads': strong[true_muts].sum().values,
    'placenta_AD': placenta['AD_alt'].values.astype(int),
    'placenta_DP': placenta['DP'].values.astype(int),
    'tree_assignment': tree.values,
    'lineage_class': lineage_class,
    'site_background': site_bg,
}).set_index('mutation_id')
summary.to_csv(os.path.join(path_results, 'LINEAGE_SUMMARY.tsv'), sep='\t')

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
genotypes[f'state_af{int(AF_MAIN * 100)}'] = state.ravel()
genotypes.to_csv(os.path.join(path_results, 'GENOTYPES_TRUE.tsv.gz'), sep='\t', index=False)
