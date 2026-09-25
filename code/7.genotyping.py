"""
Binomial genotyping of the final callset, and spatial analyses on genotype calls.

Pericentromeric sites are dropped first. Two mutation sets are then defined
(in_robust, in_backbone), each cut x mutation cell is genotyped against the
sequencing error rate, and the septum/ventricle comparisons are re-run treating
underpowered zeros as missing rather than as absences.
"""

import os
import itertools
import numpy as np
import pandas as pd
from scipy.stats import binom, mannwhitneyu, wilcoxon
from sklearn.metrics import pairwise_distances

import statsmodels.formula.api as smf


##


# Approximate hg38 centromere intervals (Mb), padded by PAD below. Acrocentric
# p-arms (13, 14, 15, 21, 22) are counted from 0. Replace with the UCSC
# centromere/gap track if exact boundaries matter.
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

ERR = 5e-4          # per-base error rate for the binomial genotyper
MIN_DP = 30         # minimum depth for a read count to support a mutation
P_PRESENT = 0.01    # significance for a present call
AF_EXCLUDE = 0.10   # absence calls must exclude a clone at this AF (95% power)
PRIOR_BASES = 1000  # shrinkage weight of the substitution-class rate, in bases
PLACENTA_AF_FRACTION = 0.2   # placenta AF above this share of heart AF = real presence
N_PERM = 2000

REGION_ABBR = {
    'Left_septum': 'LS', 'Centre_septum': 'CS', 'Right_septum': 'RS',
    'Left_Ventricle': 'LV', 'Right_Ventricle': 'RV',
}


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


def region_mean(S, a, b, labels):
    """
    Mean of S over cut pairs between region a and region b (diagonal excluded).
    """
    ia, ib = np.where(labels == a)[0], np.where(labels == b)[0]
    sub = S[np.ix_(ia, ib)]
    if a == b:
        sub = sub[~np.eye(len(ia), dtype=bool)]
    return np.nanmean(sub)


def perm_pvalue(null, obs):
    """
    One-sided permutation p-value.
    """
    return min(1, min((null >= obs).mean()))
    #return min(1, 2 * min((null >= obs).mean(), (null <= obs).mean()))


##


# Paths
path_main = '/Users/cossa/Desktop/projects/manas_heart'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results')

# Read data: final callset annotations, and the full force-called matrix (zeros included)
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


##


# 1. Drop pericentromeric sites first
muts['pericentromeric'] = [is_pericentromeric(m) for m in muts.index]
kept = muts.index[~muts['pericentromeric']]
print(f'Mutations: {len(muts)} -> {len(kept)} after dropping '
      f'{int(muts["pericentromeric"].sum())} pericentromeric sites')

heart = full.query('tissue == "heart" and mutation_id in @kept')
AD = heart.pivot(index='Sample_ID', columns='mutation_id', values='AD_alt')
DP = heart.pivot(index='Sample_ID', columns='mutation_id', values='DP')
region = samples['region'].loc[AD.index].map(REGION_ABBR)
chunk = samples['chunk'].loc[AD.index]

# (AD!=0).sum().describe()
# (DP<10).sum().describe()
# 
# top = (AD!=0).sum().sort_values(ascending=False).head(10).index
# 
# AD[top].values
# DP[top].values
# (AD[top].values / DP[top].values).round(2)




##


# 2. Mutation sets: supporting reads only count in cuts with enough depth
AD_ok = AD.where(DP >= MIN_DP)
in_robust = (AD_ok >= 3).sum() >= 1
in_backbone = (AD_ok >= 5).sum() >= 2
print(f'in_robust: {int(in_robust.sum())} | in_backbone: {int(in_backbone.sum())}')

# Per-site background, estimated from the placenta cuts. Sites whose placenta
# signal looks like real presence (above error and comparable to their heart AF)
# cannot serve as their own background, so they take their substitution class
# rate. Site estimates are shrunk towards that class rate, since ~210x of
# placenta depth only resolves ~5e-3.
cols = in_robust[in_robust].index
placenta = (
    full.query('tissue == "placenta" and mutation_id in @cols')
    .groupby('mutation_id').agg(AD_alt=('AD_alt', 'sum'), DP=('DP', 'sum'))
    .reindex(cols)
)
heart_af = AD[cols].sum() / DP[cols].sum()
placenta_af = placenta['AD_alt'] / placenta['DP']

#zero_in_placenta = placenta.query('AD_alt<=2').index
#AD.sum()[heart_af.sort_values(ascending=False).index]
#
#for mut in heart_af.sort_values(ascending=False).index:
#    df_ = (
#        pd.DataFrame({'AD':AD[mut], 'DP':DP[mut],}).join(samples)
#        .groupby('region')[['AD', 'DP']].sum()
#    )
#    if (df_['AD']==0).any() & (df_['AD']>10).any():
#        print(
#            f'Mut: {mut}', df_,
#            '\n'
#        )
#
#heart_af['chr12_51771708_C_T']
#placenta.loc['chr12_51771708_C_T']
#in_robust.loc['chr12_51771708_C_T']

carries = (
    (binom.sf(placenta['AD_alt']-1, placenta['DP'], ERR) < 0.01)
    & (placenta_af > PLACENTA_AF_FRACTION * heart_af)
)
substitution = muts['SBS6'].reindex(cols).fillna('NA')
class_rate = (
    placenta.loc[~carries, 'AD_alt'].groupby(substitution[~carries]).sum()
    / placenta.loc[~carries, 'DP'].groupby(substitution[~carries]).sum()
)
prior_rate = substitution.map(class_rate).fillna(ERR).values
site_error = np.where(
    carries.values, prior_rate,
    (placenta['AD_alt'].values + PRIOR_BASES * prior_rate) / (placenta['DP'].values + PRIOR_BASES)
).clip(min=1e-5)
print(f'\nPlacenta background: {int(carries.sum())} of {len(cols)} sites excluded as placenta-carrying')
print('  rate by substitution class: ' +
      ', '.join(f'{k} {v:.1e}' for k, v in class_rate.round(6).items()))
print(f'  per-site background: median {np.median(site_error):.1e}, '
      f'90th pct {np.percentile(site_error, 90):.1e}, max {site_error.max():.1e}')


##


# 3. Genotype every cut x mutation cell of in_robust, against its site background
A = AD[cols].values
D = DP[cols].values
p_error = binom.sf(A - 1, D, site_error[None, :])   # P(>= observed reads | background)
lod = detection_limit(D)
state = np.full(A.shape, 'nocall', dtype=object)
state[(A >= 2) & (p_error < 0.05)] = 'present'
state[A == 1] = 'weak'
state[(A == 0) & (lod <= AF_EXCLUDE)] = 'absent'

geno = pd.DataFrame(state, index=AD.index, columns=cols)
present = geno.eq('present').values
absent = geno.eq('absent').values
determinate = present | absent

print('\nGenotype states (in_robust, %d cells):' % A.size)
for s in ['present', 'weak', 'absent', 'nocall']:
    print(f'  {s:8s} {int((state == s).sum()):6d} ({100 * (state == s).mean():4.1f}%)')
print(f'  detection limit (AF, 95% power): median {np.median(lod):.2f}')

# Long-format genotype table
out = (
    pd.DataFrame({
        'Sample_ID': np.repeat(AD.index.values, len(cols)),
        'mutation_id': np.tile(cols.values, len(AD.index)),
        'AD_alt': A.ravel(), 'DP': D.ravel(),
        'state': state.ravel(), 'detection_limit_AF': lod.ravel().round(3),
        'p_error': p_error.ravel(),
        'site_background': np.tile(site_error, len(AD.index)),
    })
    .assign(
        region=lambda x: x['Sample_ID'].map(samples['region']),
        chunk=lambda x: x['Sample_ID'].map(samples['chunk']),
        in_backbone=lambda x: x['mutation_id'].map(in_backbone),
    )
)
out.to_csv(os.path.join(path_results, 'GENOTYPES.tsv.gz'), sep='\t', index=False)

df_ = out[['Sample_ID', 'mutation_id', 'AD_alt', 'DP', 'state', 'detection_limit_AF']].assign(AF=lambda x: x['AD_alt']/x['DP'])

df_ = df_.query('state=="present" and AF > detection_limit_AF')

df_ = df_.pivot(values='AF', columns='mutation_id', index='Sample_ID').fillna(0)

D = pairwise_distances((df_>0).values, metric='jaccard')

D = pd.DataFrame(D, index=df_.index, columns=df_.index)

regions = samples



(
    pd.DataFrame({'pericentromeric': muts['pericentromeric']})
    .join(pd.DataFrame({'in_robust': in_robust, 'in_backbone': in_backbone}))
    .fillna(False)
    # .to_csv(os.path.join(path_results, 'MUTATION_SETS.tsv'), sep='\t')
)


##


# 4. Region-level genotypes: a region is 'absent' only if every cut there is a
# powered absence, i.e. underpowered zeros never count as absences
rows = []
for r in REGION_ABBR.values():
    idx = np.where(region.values == r)[0]
    rows.append(pd.Series(
        np.where(present[idx].any(0), 'present',
                 np.where(absent[idx].all(0), 'absent', 'nocall')),
        index=cols, name=r
    ))
region_geno = pd.concat(rows, axis=1)
print('\nRegion-level genotypes (%d mutations x 5 regions):' % len(cols))
print(region_geno.apply(pd.Series.value_counts).fillna(0).astype(int).to_string())

determinate_regions = region_geno.ne('nocall').sum(1)
restricted = region_geno.apply(
    lambda x: (x == 'present').sum() >= 1 and (x == 'absent').sum() >= 1, axis=1
)
print(f'\nMutations with all 5 regions determinate: {int((determinate_regions == 5).sum())}')
print(f'Mutations present in >=1 region and absent (with power) in >=1 other: {int(restricted.sum())}')
if restricted.any():
    print(region_geno[restricted].to_string())


##


# 5. Rare lineages from present calls only (>=2 reads), and their sharing
n_present = present.sum(0)
rare = (n_present >= 2) & (n_present <= 6)
print(f'\nRare tier on present calls (2-6 cuts): {int(rare.sum())} mutations '
      f'(vs {int((((AD[cols] > 0).sum() >= 2) & ((AD[cols] > 0).sum() <= 6)).sum())} using >=1 read)')

xyz = pd.read_csv(os.path.join(path_data, 'Heart_final_coorindates_135.csv')).set_index('name')
coords = xyz.loc[AD.index, ['x', 'y', 'z']].values
dist = pairwise_distances(coords)
iu = np.triu_indices(len(AD.index), 1)

shared = (present[:, rare].astype(float) @ present[:, rare].astype(float).T)
# pairs are comparable only where both cuts could be called: normalise by the
# number of rare mutations determinate in both cuts
comparable = (determinate[:, rare].astype(float) @ determinate[:, rare].astype(float).T)
rate = np.divide(shared, comparable, out=np.full_like(shared, np.nan), where=comparable > 0)

is_septum = np.isin(region.values, ['LS', 'CS', 'RS'])
crosses = is_septum[iu[0]] ^ is_septum[iu[1]]
touches_rv = (region.values[iu[0]] == 'RV') | (region.values[iu[1]] == 'RV')
pair_type = np.where(
    crosses, np.where(touches_rv, 'septum-RV', 'septum-LV'),
    np.where(chunk.values[iu[0]] == chunk.values[iu[1]], 'same chunk',
             np.where(~is_septum[iu[0]] & ~is_septum[iu[1]], 'LV-RV',
                      'septum-septum/same region')))
pairs = pd.DataFrame({
    'shared': shared[iu], 'comparable': comparable[iu], 'rate': rate[iu],
    'dist': dist[iu], 'type': pair_type,
}).query('comparable > 0')
print('\nRare-lineage sharing by pair type (present calls, power-normalised):')
print(
    pairs.groupby('type')
    .agg(n_pairs=('shared', 'size'), shared=('shared', 'mean'),
         comparable=('comparable', 'mean'), rate=('rate', 'mean'), dist=('dist', 'mean'))
    .round(4).sort_values('rate', ascending=False).to_string()
)

model = smf.poisson(
    'shared ~ np.log(dist + 1) + C(type, Treatment("septum-LV"))',
    data=pairs, offset=np.log(pairs['comparable'])
).fit(disp=0, cov_type='HC0')
print('\nShared rare lineages ~ log(distance) + pair type (offset: comparable mutations):')
print(pd.DataFrame({'ratio_vs_septum_LV': np.exp(model.params), 'p': model.pvalues}).round(3).to_string())

# Cuts differ in how much they share with everything: a cut holding a clonal
# patch shares more with every other cut. Normalise each ventricle cut by its
# own sharing with the other ventricle cuts before comparing LV and RV.
np.fill_diagonal(rate, np.nan)
is_vent = ~is_septum
prop = []
for i in np.where(is_vent)[0]:
    others = is_vent & (np.arange(len(region)) != i)
    with_septum = np.nanmean(rate[i, is_septum])
    with_ventricles = np.nanmean(rate[i, others])
    prop.append({
        'Sample_ID': AD.index[i], 'region': region.values[i],
        'with_septum': with_septum, 'with_ventricles': with_ventricles,
        'ratio': with_septum / with_ventricles if with_ventricles > 0 else np.nan,
    })
prop = pd.DataFrame(prop)
lv, rv = prop.query('region == "LV"'), prop.query('region == "RV"')
print('\nPer-ventricle-cut sharing with the septum:')
print(f'  raw:        LV median {lv["with_septum"].median():.4f} | RV median {rv["with_septum"].median():.4f} | '
      f'MWU p = {mannwhitneyu(lv["with_septum"].dropna(), rv["with_septum"].dropna())[1]:.3f}')
print(f'  normalised: LV median {lv["ratio"].median():.2f} | RV median {rv["ratio"].median():.2f} | '
      f'MWU p = {mannwhitneyu(lv["ratio"].dropna(), rv["ratio"].dropna())[1]:.3f}')


##


# 6. Septum x ventricle comparisons, underpowered zeros treated as missing
rng = np.random.default_rng(0)
sim_present = rate                                   # power-normalised sharing, all tiers below
sims = {'rare (2-6 present calls)': rate}

for name, sel in [('intermediate (7-25)', (n_present >= 7) & (n_present <= 25)),
                  ('common (>25)', n_present > 25)]:
    sh = present[:, sel].astype(float) @ present[:, sel].astype(float).T
    cp = determinate[:, sel].astype(float) @ determinate[:, sel].astype(float).T
    sims[name] = np.divide(sh, cp, out=np.full_like(sh, np.nan), where=cp > 0)

# AF profile similarity uses AF directly (no genotype threshold), missing where DP < MIN_DP
AF = np.where(D >= MIN_DP, A / np.maximum(D, 1), np.nan)
AF_filled = np.nan_to_num(AF)
sims['AF cosine (in_robust)'] = 1 - pairwise_distances(AF_filled, metric='cosine')

labels = region.values
is_vent = np.isin(labels, ['LV', 'RV'])
is_lat = np.isin(labels, ['LS', 'RS'])

print('\nSeptum section x ventricle (cut-level permutation p):')
rows = []
for name, S in sims.items():
    row = {'measure': name}
    for s in ['LS', 'CS', 'RS']:
        for v in ['LV', 'RV']:
            row[f'{s}-{v}'] = round(region_mean(S, s, v, labels), 4)
    for s in ['LS', 'CS', 'RS']:
        obs = region_mean(S, s, 'LV', labels) - region_mean(S, s, 'RV', labels)
        null = []
        for _ in range(N_PERM):
            lab = labels.copy()
            lab[is_vent] = rng.permutation(labels[is_vent])
            null.append(region_mean(S, s, 'LV', lab) - region_mean(S, s, 'RV', lab))
        row[f'p({s}: LV vs RV)'] = round(perm_pvalue(np.array(null), obs), 3)
    for v in ['LV', 'RV']:
        obs = region_mean(S, 'LS', v, labels) - region_mean(S, 'RS', v, labels)
        null = []
        for _ in range(N_PERM):
            lab = labels.copy()
            lab[is_lat] = rng.permutation(labels[is_lat])
            null.append(region_mean(S, 'LS', v, lab) - region_mean(S, 'RS', v, lab))
        row[f'p(LS vs RS: {v})'] = round(perm_pvalue(np.array(null), obs), 3)
    rows.append(row)
print(pd.DataFrame(rows).set_index('measure').to_string())


##

# 7. Sharing between regions, by paired comparison rather than a matrix-wide
# correction. Cuts differ in their overall level of sharing (depth, and whether
# they hold a clonal patch), so each contrast holds one side fixed: the cut on
# the fixed side cancels exactly, and the partner cuts are divided by their own
# overall level so their differences cancel too.
rate_matrix = np.divide(shared, comparable,
                        out=np.full(shared.shape, np.nan), where=comparable > 0)
np.fill_diagonal(rate_matrix, np.nan)

similarity = {'shared rare lineages': rate_matrix}
for name, members in [('AF cosine, in_robust', in_robust), ('AF cosine, in_backbone', in_backbone)]:
    sub = members[members].index
    af = (AD[sub] / DP[sub].clip(lower=1)).values
    S = 1 - pairwise_distances(af, metric='cosine')
    np.fill_diagonal(S, np.nan)
    similarity[name] = S

is_sep = np.isin(region.values, ['LS', 'CS', 'RS'])
labels = region.values


def normalised(S):
    """
    Divide each column by that cut's mean similarity to all others, so a
    partner's overall level cancels.
    """
    level = np.nanmean(S, axis=1)
    return S / level[None, :]


def mean_to(S, source, targets):
    """
    Mean similarity of one cut to a group of cuts, excluding itself.
    """
    sel = targets.copy()
    sel[source] = False
    return np.nanmean(S[source, sel])


print('\nRegion-pair means (raw, no correction):')
for name, S in similarity.items():
    tab = pd.DataFrame(
        {b: {a: round(np.nanmean(S[np.ix_(labels == a, labels == b)]), 3)
             for a in ['LS', 'CS', 'RS']} for b in ['LS', 'CS', 'RS', 'LV', 'RV']}
    )
    print(f'\n  {name}:')
    print(tab.to_string())

print('\nPaired tests (partner-normalised; one row per fixed-side cut):')
for name, S in similarity.items():
    X = normalised(S)
    print(f'\n  {name}')

    # (a) septum cuts: closer to LV or to RV?
    for group, members in [('all septum', is_sep), ('LS', labels == 'LS'),
                           ('CS', labels == 'CS'), ('RS', labels == 'RS')]:
        lv = np.array([mean_to(X, i, labels == 'LV') for i in np.where(members)[0]])
        rv = np.array([mean_to(X, i, labels == 'RV') for i in np.where(members)[0]])
        ok = ~(np.isnan(lv) | np.isnan(rv))
        p = wilcoxon(lv[ok], rv[ok])[1] if ok.sum() >= 6 else np.nan
        print(f'    {group:10s} vs LV {np.nanmean(lv):.3f} | vs RV {np.nanmean(rv):.3f} | '
              f'n={int(ok.sum())} | p = {p:.3f}' if ok.sum() >= 6 else
              f'    {group:10s} vs LV {np.nanmean(lv):.3f} | vs RV {np.nanmean(rv):.3f} | '
              f'n={int(ok.sum())} (too few for a test); per-cut differences '
              f'{np.round(lv[ok] - rv[ok], 3)}')

    # (b) lateral cuts: do LS cuts resemble CS more than RS cuts do? Here the
    # compared cuts are on the free side, so their own level does not cancel:
    # each cut's similarity to CS is expressed relative to its similarity to
    # every other cut, and the CS partners are normalised as before.
    preference = []
    for i in range(len(labels)):
        to_cs = mean_to(X, i, labels == 'CS')
        to_rest = mean_to(X, i, (labels != 'CS'))
        preference.append(to_cs / to_rest if to_rest else np.nan)
    preference = np.array(preference)
    ls_vals = preference[(labels == 'LS') & ~np.isnan(preference)]
    rs_vals = preference[(labels == 'RS') & ~np.isnan(preference)]
    print(f'    preference for CS (own level divided out): LS cuts {ls_vals.mean():.3f} '
          f'(n={len(ls_vals)}) | RS cuts {rs_vals.mean():.3f} (n={len(rs_vals)}) | '
          f'Mann-Whitney p = {mannwhitneyu(ls_vals, rs_vals)[1]:.3f}')

    # (c) septal cuts: own region (different chunk) vs the other septal regions
    own, other = [], []
    for i in np.where(is_sep)[0]:
        same_region = (labels == labels[i]) & (chunk.values != chunk.values[i])
        other_region = is_sep & (labels != labels[i])
        own.append(mean_to(X, i, same_region))
        other.append(mean_to(X, i, other_region))
    own, other = np.array(own), np.array(other)
    ok = ~(np.isnan(own) | np.isnan(other))
    print(f'    own region (diff chunk) {np.nanmean(own):.3f} | other septal regions '
          f'{np.nanmean(other):.3f} | n={int(ok.sum())} | '
          f'p = {wilcoxon(own[ok], other[ok])[1]:.3f}')
