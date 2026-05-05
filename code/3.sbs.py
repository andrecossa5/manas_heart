"""
Mutational spectrum assignment to COSMIC v3 signatures, 
using SigProfilerAssignment. This is done on a pseudobulk matrix of SBS96 contexts,
built from the filtered table and annotated with trinucleotide contexts.
"""

import os
import pysam
import pandas as pd
from SigProfilerExtractor import sigpro as sig
from SigProfilerAssignment import Analyzer


##


MUT_ORDER = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
COMP = str.maketrans('ACGT', 'TGCA')


##


def revcomp(s):
    return s.translate(COMP)[::-1]


##


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


# Paths
path_main = '/Users/cossa/Desktop/projects/manas_heart'
path_data = os.path.join(path_main, 'data')
path_resources = os.path.join(path_main, 'resources')
path_filtered = os.path.join(path_main, 'results')
path_figures = os.path.join(path_main, 'figures')


##


# Read forcecall results
def main():

    # Read filtered table
    df = pd.read_csv(os.path.join(path_filtered, 'ALLELIC_TABLE_FILTERED.tsv.gz'), sep='\t')

    # Annotate contexts
    path_ref = os.path.join(path_resources, 'GRCh38.d1.vd1.fa')
    df = annotate_ctx(df, path_ref)

    # Force SBS96 column to a Categorical with all 96 canonical contexts
    bases = ['A', 'C', 'G', 'T']
    SBS96_CTX = [f'{p}[{ref}>{alt}]{n}'
                 for mut in MUT_ORDER
                 for ref, alt in [mut.split('>')]
                 for p in bases for n in bases]
    df['SBS96'] = pd.Categorical(df['SBS96'], categories=SBS96_CTX, ordered=True)


    ##

    # Build SigProfiler matrix
    path_sig = os.path.join(path_main, 'results', 'sigprofiler')
    os.makedirs(path_sig, exist_ok=True)
    (
        df[df['AF']>0]
        .drop_duplicates('mutation_id')
        .assign(cohort='cohort')
        .groupby('cohort')['SBS96']
        .value_counts()
        .reset_index()
        .rename(columns={'SBS96': 'MutationType'})
        .pivot(index='MutationType', columns='cohort', values='count')
        .fillna(0)
        .to_csv(os.path.join(path_sig, 'SBS96_matrix.tsv'), sep='\t')
    )

    # Per-chunk COSMIC v3 fit, unrestricted (with artefact signatures allowed)
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
        make_plots=True,
        cpu=-1,
    )

    # Per-chunk COSMIC v3 fit, restricted (artefact signatures excluded)
    Analyzer.cosmic_fit(
        samples=os.path.join(path_sig, 'SBS96_matrix.tsv'),
        output=os.path.join(path_sig, 'assignment_pseudobulk_no_artefacts'),
        input_type='matrix',
        context_type='96',
        genome_build='GRCh38',
        cosmic_version=3.4,
        collapse_to_SBS96=True,
        exclude_signature_subgroups=['Artifact_signatures'],
        export_probabilities=True,
        export_probabilities_per_mutation=True,
        make_plots=True,
        cpu=-1,
    )


    ##


# Run
if __name__ == '__main__':
    main()
