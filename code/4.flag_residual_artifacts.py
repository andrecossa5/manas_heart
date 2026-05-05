"""
Flag artifacts still present in filtered forcecalled results, using COSMIC SBS96 signature extraction.
"""

import os
import numpy as np
import pysam
import pandas as pd
import matplotlib
import plotting_utils as plu
import matplotlib.pyplot as plt
import scanpy as sc
matplotlib.use('macOSX')
plu.set_rcParams()


##


sb6_colors = {
    "C>A": "#03BDEF",
    "C>G": "#010101",
    "C>T": "#E42926",
    "T>A": "#CBCACA",
    "T>C": "#A2CF63",
    "T>G": "#ECC7C5",
}
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


def calculate_sbs96(df, context='SBS96'):

    bases = ['A', 'C', 'G', 'T']
    total = len(df)
    groups = ['SBS6', context]
    counts = (
        df.groupby(groups)
        .size()
        .div(total)
        .reset_index(name='fraction')
    )

    # Build complete index of all 96 contexts
    ctx_per_mut = {}
    for mut in MUT_ORDER:
        ref, alt = mut[0], mut[2]
        if context == 'SBS96':
            ctx_per_mut[mut] = sorted([f"{p}[{ref}>{alt}]{n}" for p in bases for n in bases])
        else:
            ctx_per_mut[mut] = sorted([f"{p}{ref}{n}" for p in bases for n in bases])

    full_idx = pd.MultiIndex.from_tuples(
        [(mut, ctx) for mut in MUT_ORDER for ctx in ctx_per_mut[mut]],
        names=['SBS6', context]
    )
    counts = (
        counts.set_index(groups)['fraction']
        .reindex(full_idx, fill_value=0)
        .reset_index()
    )

    return counts


##


def mut_profile(df=None, counts=None, context='SBS96', figsize=(12, 3), legend_kwargs={}) -> matplotlib.figure.Figure:
    """
    Plot raw fraction of MT-SNVs across SBS96 (or 3nt) contexts,
    stratified by mutation type (one axis each).
    """

    df = df.drop_duplicates('mutation_id') if df is not None else None

    if counts is None:
        total = len(df)
        counts = calculate_sbs96(df, context=context)
    else:
        total = None

    fig, axs = plt.subplots(
        1, len(MUT_ORDER), figsize=figsize, sharey=True,
        constrained_layout=True
    )

    for i, mut in enumerate(MUT_ORDER):
        ax = axs[i]
        df_ = counts.query('SBS6 == @mut')
        x_order = sorted(df_[context].unique())
        plu.bar(
            df_, x=context, y='fraction',
            color=sb6_colors[mut],
            x_order=x_order,
            width=0.8, alpha=1.0, edgecolor=None,
            with_label=False, ax=ax
        )
        n_mut = int(round(df_['fraction'].sum() * total)) if total is not None else None
        plu.format_ax(
            ax, xlabel=mut, rotx=90,
            title=f'n: {n_mut}' if n_mut is not None else '',
            ylabel='Fraction of total SBSs' if i == 0 else '',
            reduced_spines=True, xticks_size=6
        )

    plu.add_legend(
        ax=axs[-1],
        colors={'H': '#444444', 'L': '#bbbbbb'},
        label='Strand', ncols=1,
        loc='upper left', bbox_to_anchor=(1, 1),
        **legend_kwargs
    )

    return fig


##


# Paths
path_main = '/Users/cossa/Desktop/projects/manas_heart'
path_data = os.path.join(path_main, 'data')
path_filtered = os.path.join(path_main, 'results')
path_resources = os.path.join(path_main, 'resources')
path_assignmet = os.path.join(path_main, 'results/sigprofiler/assignment_pseudobulk/Assignment_Solution/Activities')
path_figures = os.path.join(path_main, 'figures')


##


# Read filtered forceallcelled results
df = pd.read_csv(os.path.join(path_filtered, 'ALLELIC_TABLE_FILTERED.tsv.gz'), sep='\t')

# Annotate context
df = annotate_ctx(df, os.path.join(path_resources, 'GRCh38.d1.vd1.fa'))

# COSMIC v3.x signatures flagged as possible sequencing/library artefacts
ARTEFACT_SIGS = [
    'SBS27', 'SBS43',
    'SBS45', 'SBS46', 'SBS47', 'SBS48', 'SBS49', 'SBS50',
    'SBS51', 'SBS52', 'SBS53', 'SBS54', 'SBS55', 'SBS56', 'SBS57',
    'SBS58', 'SBS59', 'SBS60', 'SBS95',
]

# Per-sample COSMIC activities
ctx_probabilities = pd.read_csv(
    os.path.join(path_assignmet, 'Decomposed_MutationType_Probabilities.txt'), sep='\t'
)

# Get crisp assignment of each context
assert not ctx_probabilities['SBS96'].any()
ctx_probabilities.drop(columns=['SBS96', 'Sample Names'], inplace=True)
signatures = ctx_probabilities.iloc[:,1:].columns
ctx_probabilities['assignment'] = (
    ctx_probabilities.iloc[:,1:]
    .apply(lambda x: signatures[x.argmax()], axis=1)
)
ctx_probabilities.rename(columns={'MutationType': 'SBS96'}, inplace=True)

# Merge with allelic table
df = df.merge(ctx_probabilities[['assignment', 'SBS96']], on='SBS96', how='left')
df['artifact_flag'] = (
    df['assignment'].apply(lambda x: 'Artefact' if x in ARTEFACT_SIGS else 'No artefact')
)

# (
#     df.loc[
#     (df['AF']>0) & \
#     (~df['assignment'].isin(ARTEFACT_SIGS))]
#     ['mutation_id'].nunique()
# )


##


# Only positive muts
df = df.query('AF>0').copy()

fig, ax = plt.subplots(figsize=(5, 5))

colors = plu.create_palette(df, 'assignment', palette=sc.pl.palettes.default_20)
colors['Unassigned'] = '#cccccc'
plu.bb_plot(df, 'Sample_ID', 'assignment', ax=ax, categorical_cmap=colors)
plu.format_ax(ax, reduced_spines=True, yticks_size=3)
plu.add_legend(
    ax=ax, colors=colors, label='COSMIC SBS',
    loc='upper left', bbox_to_anchor=(1, 1), ncols=1,
    artists_size=7, label_size=8, ticks_size=6
)
fig.subplots_adjust(right=0.7, left=.2, top=.8, bottom=.2)
fig.savefig(os.path.join(path_figures, 'forcecall_COSMIC_assignment.pdf'))


##


fig, ax = plt.subplots(figsize=(5, 5))

colors = plu.create_palette(df, 'artifact_flag', palette='Set1')
plu.bb_plot(df, 'Sample_ID', 'artifact_flag', ax=ax, categorical_cmap=colors)
plu.format_ax(ax, reduced_spines=True, yticks_size=3)
plu.add_legend(
    ax=ax, colors=colors, label='Artifact',
    loc='upper left', bbox_to_anchor=(1, 1), ncols=1,
    artists_size=7, label_size=8, ticks_size=6
)
fig.subplots_adjust(right=0.7, left=.2, top=.8, bottom=.2)
fig.savefig(os.path.join(path_figures, 'forcecall_artifact_flag.pdf'))


##


# Viz spectra
fig = mut_profile(df.query('artifact_flag == "Artefact"'), context='SBS96', figsize=(12, 3))
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'forcecall_artifact_spectrum.pdf'))

# Write
(
    df.query('artifact_flag == "No artefact"')
    .to_csv(os.path.join(path_filtered, 'ALLELIC_TABLE_NO_ARTIFACTS.tsv.gz'), sep='\t', index=False)
)


##


