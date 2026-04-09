#!/usr/bin/env python3
"""
Add SBS6 and SBS96 trinucleotide context (pyrimidine) to ALL_FILTERED.tsv.gz.

Indels receive NA for both columns. SNVs on purine REF strands are
reverse-complemented so all contexts are reported in pyrimidine notation.

Usage:
    python add_context.py --input ALL_FILTERED.tsv.gz --ref GRCh38.d1.vd1.fa \\
                          --output ALL_FILTERED_ctx.tsv.gz [--chunksize N] [--n-workers N]
"""

import sys
import gzip
import argparse
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import pysam

COMP = str.maketrans('ACGT', 'TGCA')

def revcomp(s):
    return s.translate(COMP)[::-1]


# ---------------------------------------------------------------------------
# Per-record annotation (called inside worker processes)
# ---------------------------------------------------------------------------

_fasta = None  # process-local pysam handle

def _init_worker(ref_path):
    global _fasta
    _fasta = pysam.FastaFile(ref_path)


def _annotate_records(records):
    """
    records : list of (chrom, pos, ref, alt)   pos is 1-based
    returns : list of (sbs6, sbs96)             None for indels / missing context
    """
    results = []
    for chrom, pos, ref, alt in records:
        ref = ref.upper()
        alt = alt.upper()

        if len(ref) != 1 or len(alt) != 1:
            results.append((None, None))
            continue

        try:
            # pysam fetch is 0-based half-open; pos-2 gives [pos-1, pos, pos+1] in 1-based
            tri = _fasta.fetch(chrom, pos - 2, pos + 1).upper()
        except Exception:
            results.append((None, None))
            continue

        if len(tri) != 3 or 'N' in tri:
            results.append((None, None))
            continue

        if ref in ('C', 'T'):
            ctx_ref, ctx_alt, ctx_tri = ref, alt, tri
        else:
            ctx_ref = revcomp(ref)
            ctx_alt = revcomp(alt)
            ctx_tri = revcomp(tri)

        results.append((
            f"{ctx_ref}>{ctx_alt}",
            f"{ctx_tri[0]}[{ctx_ref}>{ctx_alt}]{ctx_tri[2]}",
        ))

    return results


# ---------------------------------------------------------------------------
# Chunk-level helper (single-process fallback and pool target)
# ---------------------------------------------------------------------------

def _annotate_chunk_records(args):
    """Pool-safe wrapper: receives (records,) and returns annotations."""
    (records,) = args
    return _annotate_records(records)


def annotate_chunk(df, ref_path, executor, n_workers):
    records = list(zip(df['CHROM'], df['POS'].astype(int), df['REF'], df['ALT']))

    if executor is None or n_workers == 1:
        annotations = _annotate_records(records)
    else:
        # Split records evenly across workers and submit
        size = max(1, len(records) // n_workers)
        batches = [records[i:i + size] for i in range(0, len(records), size)]
        annotations = []
        for batch_result in executor.map(_annotate_chunk_records, [(b,) for b in batches]):
            annotations.extend(batch_result)

    sbs6, sbs96 = zip(*annotations) if annotations else ([], [])
    df = df.copy()
    df['SBS6']  = list(sbs6)
    df['SBS96'] = list(sbs96)
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input',     '-i', required=True,  help='Input TSV.gz (ALL_FILTERED.tsv.gz)')
    p.add_argument('--ref',       '-r', required=True,  help='Indexed reference FASTA (.fa + .fai)')
    p.add_argument('--output',    '-o', required=True,  help='Output TSV.gz path')
    p.add_argument('--chunksize', type=int, default=500_000, metavar='N',
                   help='Rows per pandas chunk [%(default)s]')
    p.add_argument('--n-workers', type=int, default=1,  metavar='N',
                   help='Parallel worker processes [%(default)s]')
    return p.parse_args()


def main():
    args = parse_args()

    executor = None
    if args.n_workers > 1:
        executor = ProcessPoolExecutor(
            max_workers=args.n_workers,
            initializer=_init_worker,
            initargs=(args.ref,),
        )
    else:
        _init_worker(args.ref)   # init in main process

    first  = True
    n_rows = 0

    try:
        with gzip.open(args.output, 'wt') as out_fh:
            for chunk in pd.read_csv(args.input, sep='\t', chunksize=args.chunksize):
                chunk = annotate_chunk(chunk, args.ref, executor, args.n_workers)
                chunk.to_csv(out_fh, sep='\t', index=False, header=first)
                first   = False
                n_rows += len(chunk)
                print(f'  {n_rows:>12,} rows processed', file=sys.stderr, end='\r')
    finally:
        if executor is not None:
            executor.shutdown(wait=False)

    print(f'\nDone. Wrote {n_rows:,} rows to {args.output}', file=sys.stderr)


if __name__ == '__main__':
    main()
