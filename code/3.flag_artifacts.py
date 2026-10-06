"""
Trinucleotide context of the filtered force-called SNVs, COSMIC v3.4 signature fit of their pseudobulk
SBS96 spectrum (SigProfilerAssignment), and flag of the SNVs whose context is assigned to a sequencing /
library artefact signature. Keeps the positive SNVs only.

Reads results/ALLELIC_TABLE_FILTERED.tsv.gz and resources/GRCh38.d1.vd1.fa; writes results/sigprofiler/ and
results/ALLELIC_TABLE_FILTERED_ANNOTATED.tsv.gz.
"""

import os
import pysam
import pandas as pd
from SigProfilerAssignment import Analyzer


##


MUT_ORDER = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
COMP = str.maketrans('ACGT', 'TGCA')

# COSMIC v3.x signatures flagged as possible sequencing/library artefacts
ARTEFACT_SIGS = [
    'SBS27', 'SBS43',
    'SBS45', 'SBS46', 'SBS47', 'SBS48', 'SBS49', 'SBS50',
    'SBS51', 'SBS52', 'SBS53', 'SBS54', 'SBS55', 'SBS56', 'SBS57',
    'SBS58', 'SBS59', 'SBS60', 'SBS95',
]


##


def revcomp(s):
    return s.translate(COMP)[::-1]


def annotate_ctx(df, path_ref):
    """
    Annotate each row of `df` with SBS6 and SBS96 trinucleotide context
    (pyrimidine notation). Indels / missing contexts get None.

    Requires columns: CHROM, POS (1-based), REF, ALT.
    """
    fasta = pysam.FastaFile(path_ref)
    sbs6, sbs96 = [], []

    for chrom, pos, ref, alt in zip(df['CHROM'], df['POS'].astype(int), df['REF'], df['ALT']):
        ref, alt = ref.upper(), alt.upper()

        if len(ref) != 1 or len(alt) != 1:
            sbs6.append(None); sbs96.append(None); continue

        try:
            tri = fasta.fetch(chrom, pos - 2, pos + 1).upper()
        except Exception:
            sbs6.append(None); sbs96.append(None); continue

        if len(tri) != 3 or 'N' in tri:
            sbs6.append(None); sbs96.append(None); continue

        if ref in ('C', 'T'):
            ctx_ref, ctx_alt, ctx_tri = ref, alt, tri
        else:
            ctx_ref, ctx_alt, ctx_tri = revcomp(ref), revcomp(alt), revcomp(tri)

        sbs6.append(f"{ctx_ref}>{ctx_alt}")
        sbs96.append(f"{ctx_tri[0]}[{ctx_ref}>{ctx_alt}]{ctx_tri[2]}")

    df = df.copy()
    df['SBS6'] = sbs6
    df['SBS96'] = sbs96
    return df


##


# Paths
path_main = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))     # repository root
path_results = os.path.join(path_main, 'results')
path_ref = os.path.join(path_main, 'resources', 'GRCh38.d1.vd1.fa')
path_sig = os.path.join(path_results, 'sigprofiler')


def main():

    df = pd.read_csv(os.path.join(path_results, 'ALLELIC_TABLE_FILTERED.tsv.gz'), sep='\t')
    df = annotate_ctx(df, path_ref)

    # Pseudobulk SBS96 matrix of the positive SNVs, with all 96 canonical contexts
    bases = ['A', 'C', 'G', 'T']
    SBS96_CTX = [f'{p}[{ref}>{alt}]{n}'
                 for mut in MUT_ORDER
                 for ref, alt in [mut.split('>')]
                 for p in bases for n in bases]
    os.makedirs(path_sig, exist_ok=True)
    (
        df[df['AF']>0]
        .drop_duplicates('mutation_id')
        .assign(cohort='cohort', SBS96=lambda x: pd.Categorical(x['SBS96'], categories=SBS96_CTX, ordered=True))
        .groupby('cohort')['SBS96']
        .value_counts()
        .reset_index()
        .rename(columns={'SBS96': 'MutationType'})
        .pivot(index='MutationType', columns='cohort', values='count')
        .fillna(0)
        .to_csv(os.path.join(path_sig, 'SBS96_matrix.tsv'), sep='\t')
    )

    # COSMIC v3.4 fit, artefact signatures allowed
    Analyzer.cosmic_fit(
        samples=os.path.join(path_sig, 'SBS96_matrix.tsv'),
        output=os.path.join(path_sig, 'assignment_pseudobulk'),
        input_type='matrix',
        context_type='96',
        genome_build='GRCh38',
        cosmic_version=3.4,
        collapse_to_SBS96=True,
        export_probabilities=True,
        export_probabilities_per_mutation=True,
        make_plots=False,
        cpu=-1,
    )

    # Crisp signature assignment of each context (most likely signature), and artefact flag
    ctx_probabilities = pd.read_csv(
        os.path.join(path_sig, 'assignment_pseudobulk', 'Assignment_Solution', 'Activities',
                     'Decomposed_MutationType_Probabilities.txt'), sep='\t'
    )
    assert not ctx_probabilities['SBS96'].any()
    ctx_probabilities.drop(columns=['SBS96', 'Sample Names'], inplace=True)
    signatures = ctx_probabilities.iloc[:,1:].columns
    ctx_probabilities['assignment'] = (
        ctx_probabilities.iloc[:,1:]
        .apply(lambda x: signatures[x.argmax()], axis=1)
    )
    ctx_probabilities.rename(columns={'MutationType': 'SBS96'}, inplace=True)

    df = df.merge(ctx_probabilities[['assignment', 'SBS96']], on='SBS96', how='left')
    df['artifact_flag'] = (
        df['assignment'].apply(lambda x: 'Artefact' if x in ARTEFACT_SIGS else 'No artefact')
    )

    # Only positive SNVs
    df = df.query('AF>0').copy()
    n_clean = df.loc[df['artifact_flag'] == 'No artefact', 'mutation_id'].nunique()
    print(f'ALLELIC_TABLE_FILTERED_ANNOTATED: {df["mutation_id"].nunique()} SNVs, {n_clean} not flagged as artefacts')
    df.to_csv(os.path.join(path_results, 'ALLELIC_TABLE_FILTERED_ANNOTATED.tsv.gz'), sep='\t', index=False)


if __name__ == '__main__':
    main()
