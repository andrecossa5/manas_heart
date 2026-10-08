"""
Lineage inputs of the final figure, from the final callset (results/ALLELIC_TABLE_FINAL.tsv.gz and
results/ALLELIC_TABLE.tsv.gz): TRUE mutation set, placenta counts, per-cell read counts, lineage
classes, and the read-level region-vs-rest enrichment test.

TRUE mutation: >=3 alt reads in >=2 LCM cuts, with >=5 alt reads summed over one chunk holding
such a cut; pericentromeric sites dropped.

Outputs (results/): LINEAGE_SUMMARY.tsv, GENOTYPES_TRUE.tsv.gz, REGION_ENRICHMENT.tsv.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import binom
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
ERROR_FLOOR = 2e-3      # minimum rest-of-heart VAF in the enrichment test: about the sequencing
                        # error rate (pooled placenta alt-read rate 1.99e-3); without it a region
                        # scores p = 0 when the rest of the heart has no alt read

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


# 2. Placenta counts and substitution class, reported in the summary
placenta = (
    full.query('tissue == "placenta" and mutation_id in @true_muts')
    .groupby('mutation_id').agg(AD_alt=('AD_alt', 'sum'), DP=('DP', 'sum'))
    .reindex(true_muts).fillna(0)
)
substitution = muts['SBS6'].reindex(true_muts).fillna('NA')


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


# 4. Region vs rest of the heart, on read counts. For every
# SNV x region, a read-level one-sided binomial: region AD/DP against the rest-of-heart VAF
# (floored at ERROR_FLOOR); BH over all tests. Every read counts as a replicate, so the
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
            p_binom=binom.sf(a_r - 1, d_r, max(af_o, ERROR_FLOOR)),
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
}).set_index('mutation_id')
summary.to_csv(os.path.join(path_results, 'LINEAGE_SUMMARY.tsv'), sep='\t')

genotypes = pd.DataFrame({
    'Sample_ID': np.repeat(AD.index.values, len(true_muts)),
    'mutation_id': np.tile(true_muts.values, len(AD.index)),
    'AD_alt': A.ravel().astype(int), 'DP': D.ravel().astype(int),
    'AF': AF.ravel().round(4),
}).assign(
    region=lambda x: x['Sample_ID'].map(samples['region']),
    chunk=lambda x: x['Sample_ID'].map(samples['chunk']),
)
genotypes.to_csv(os.path.join(path_results, 'GENOTYPES_TRUE.tsv.gz'), sep='\t', index=False)
