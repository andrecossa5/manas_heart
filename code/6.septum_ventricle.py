"""
Septal LCM samples against the ventricles, behind panel f and supp fi-fiii of the final figure.

For every septal sample, the cosine distance of its VAF profile to the pooled LV and to the pooled RV
profile, over all SNVs, the shared SNVs (pre-gastrulation + other shared) and the heart-specific SNVs;
and its mean Euclidean distance in space to the LV and to the RV samples. Per septal group, the number
of samples closer to LV and a two-sided paired Wilcoxon test on d(RV) - d(LV) (samples are treated as
independent; the shared set was defined after seeing the data).

Outputs (results/): SEPTUM_LV_CUTS.tsv, SEPTUM_LV_CUTS_TESTS.tsv.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.metrics import pairwise_distances


##


REGION_ABBR = {
    'Left_septum': 'LS', 'Centre_septum': 'CS', 'Right_septum': 'RS',
    'Left_Ventricle': 'LV', 'Right_Ventricle': 'RV',
}
GROUPS = ['Septum', 'LV', 'RV']
SHOW_GROUPS = ['Septum', 'LS', 'CS', 'RS']
SPACE = 'Physical space'


##


path_main = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))     # repository root
path_results = os.path.join(path_main, 'results')

geno = pd.read_csv(os.path.join(path_results, 'GENOTYPES_TRUE.tsv.gz'), sep='\t')
summary = pd.read_csv(os.path.join(path_results, 'LINEAGE_SUMMARY.tsv'), sep='\t').set_index('mutation_id')
xyz = pd.read_csv(os.path.join(path_main, 'data', 'Heart_final_coorindates_135.csv')).set_index('name')

geno['reg'] = geno['region'].map(REGION_ABBR)
AD_df = geno.pivot(index='Sample_ID', columns='mutation_id', values='AD_alt')
DP_df = geno.pivot(index='Sample_ID', columns='mutation_id', values='DP')
cuts, muts = AD_df.index, AD_df.columns
A, D = AD_df.values.astype(float), DP_df.values.astype(float)
AF = A / D
reg = geno.drop_duplicates('Sample_ID').set_index('Sample_ID').loc[cuts, 'reg'].values
coords = xyz.loc[cuts, ['x', 'y', 'z']]
klass = summary['lineage_class'].reindex(muts).values
group = np.where(np.isin(reg, ['LS', 'CS', 'RS']), 'Septum', reg)
group_idx = {g: np.where(group == g)[0] for g in GROUPS}
print(f'{len(muts)} SNVs x {len(cuts)} cuts | ' +
      ' | '.join(f'{g} {len(i)} cuts' for g, i in group_idx.items()))


def pooled(rows, cols):
    """
    Pooled VAF, sum AD / sum DP over the given cuts, for the given SNVs.
    """
    return A[np.ix_(rows, cols)].sum(0) / D[np.ix_(rows, cols)].sum(0)


##


all_snvs = np.arange(len(muts))
lv_rows, rv_rows = group_idx['LV'], group_idx['RV']
septal = group_idx['Septum']
sreg = reg[septal]
grp_cuts = {'Septum': np.arange(len(septal)), 'LS': np.where(sreg == 'LS')[0],
            'CS': np.where(sreg == 'CS')[0], 'RS': np.where(sreg == 'RS')[0]}
lv_vaf, rv_vaf = pooled(lv_rows, all_snvs), pooled(rv_rows, all_snvs)

SNV_SETS = {'All SNVs': all_snvs,
            'Pre-gastrulation + Other shared': np.where(klass != 'Heart-specific')[0],
            'Heart-specific': np.where(klass == 'Heart-specific')[0]}


def cosine_to_ventricles(cols):
    """
    Cosine distance of each septal cut to the pooled LV and RV profile, on the given SNVs.
    """
    d_ = pairwise_distances(AF[septal][:, cols], np.vstack([lv_vaf, rv_vaf])[:, cols], metric='cosine')
    return np.nan_to_num(d_, nan=1.0)


cut_tab, test_tab, closer_counts = [], [], {}
for set_name, cols in SNV_SETS.items():
    full = cosine_to_ventricles(cols)
    cut_tab.append(pd.DataFrame({'snv_set': set_name, 'n_snvs': len(cols), 'Sample_ID': cuts[septal],
                                 'region': sreg, 'd_LV': full[:, 0], 'd_RV': full[:, 1],
                                 'closer_LV': full[:, 0] < full[:, 1]}))
    closer_counts[set_name] = int((full[:, 0] < full[:, 1]).sum())
    for g in SHOW_GROUPS:
        k = grp_cuts[g]
        diff = full[k, 1] - full[k, 0]
        test_tab.append(dict(snv_set=set_name, septal_group=g, n_cuts=len(k), cuts_closer_LV=int((diff > 0).sum()),
                             p_wilcoxon=wilcoxon(diff)[1]))

# Sanity check on the septal samples closer to LV (of 36)
expected = {'All SNVs': 25, 'Pre-gastrulation + Other shared': 26, 'Heart-specific': 24}
assert closer_counts == expected, f'septal samples closer to LV changed: {closer_counts}'

# Physical space: each septal sample's mean Euclidean distance (coordinate units) to the LV samples
# and to the RV samples; same paired test
xyz_s = coords.values[septal]
d_lv_sp = pairwise_distances(xyz_s, coords.values[lv_rows]).mean(1)
d_rv_sp = pairwise_distances(xyz_s, coords.values[rv_rows]).mean(1)
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

print('\nSeptal samples closer to LV (cosine on raw VAF, pooled LV and RV profiles; physical space: mean distance)')
print(tests.round(4).to_string(index=False))
