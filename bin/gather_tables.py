#!/usr/bin/env python3

"""
Gather per-chunk filter_mutect2 output into a single table.

Modes:
  muts   — concatenate *_filtered.tsv files into ALL_FILTERED.tsv
  stats  — concatenate *_filtered.stats files into ALL_FILTERED.stats

Usage:
    gather_tables.py --mode muts  --output ALL_FILTERED.tsv   *.tsv
    gather_tables.py --mode stats --output ALL_FILTERED.stats *.stats
"""

import sys
import argparse
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('files', nargs='+', help='Input files to gather')
    p.add_argument('--mode',   required=True, choices=['muts', 'stats'],
                   help='Gathering mode: muts or stats')
    p.add_argument('--output', '-o', required=True, help='Output file path')
    return p.parse_args()


def gather_muts(files):
    L = []
    for f in files:
        L.append(pd.read_csv(f, sep='\t'))
    out = pd.concat(L, ignore_index=True)
    out['region'] = out['chunk'].map(lambda x: x.split('.')[0])
    return out


def gather_stats(files):
    L = []
    for f in files:
        chunk = f.replace('_filtered.stats', '').split('/')[-1]
        df = pd.read_csv(f, sep='\t').assign(chunk=chunk)
        L.append(df)
    out = pd.concat(L, ignore_index=True)
    out['n_records'] = out['n_records'].astype(int)
    return out


def main():
    args = parse_args()

    if args.mode == 'muts':
        out = gather_muts(args.files)
    else:
        out = gather_stats(args.files)

    compression = 'gzip' if args.output.endswith('.gz') else None
    out.to_csv(args.output, sep='\t', index=False, compression=compression)
    print(f'Wrote {len(out):,} rows to {args.output}', file=sys.stderr)


if __name__ == '__main__':
    main()
