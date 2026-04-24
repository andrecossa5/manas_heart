#!/usr/bin/env python3

"""
Gather per-chunk filter_mutect2 output into a single table.

Modes:
  muts     — concatenate *_filtered.tsv files into ALL_FILTERED.tsv (adds 'region' col)
  stats    — concatenate *_filtered.stats files into ALL_FILTERED.stats
  allelic  — concatenate *.allelic.tsv.gz forcecall files into ALLELIC_TABLE.tsv.gz

Usage:
    gather_tables.py --mode muts    --output ALL_FILTERED.tsv.gz  *.tsv
    gather_tables.py --mode stats   --output ALL_FILTERED.stats   *.stats
    gather_tables.py --mode allelic --output ALLELIC_TABLE.tsv.gz *.allelic.tsv.gz
"""

import sys
import csv
import gzip
import argparse
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('files', nargs='+', help='Input files to gather')
    p.add_argument('--mode',   required=True, choices=['muts', 'stats', 'allelic'],
                   help='Gathering mode: muts, stats, or allelic')
    p.add_argument('--output', '-o', required=True, help='Output file path')
    return p.parse_args()


def gather_muts(files, output):
    """Stream all TSVs into a single gzip TSV — O(1) memory."""
    open_fn = gzip.open if output.endswith('.gz') else open
    n_rows = 0
    writer = None

    with open_fn(output, 'wt') as out_fh:
        for f in files:
            in_open = gzip.open if f.endswith('.gz') else open
            with in_open(f, 'rt', newline='') as in_fh:
                reader = csv.DictReader(in_fh, delimiter='\t')
                if writer is None:
                    fieldnames = reader.fieldnames + ['region']
                    writer = csv.DictWriter(out_fh, fieldnames=fieldnames,
                                            delimiter='\t')
                    writer.writeheader()
                for row in reader:
                    row['region'] = row['chunk'].split('.')[0]
                    writer.writerow(row)
                    n_rows += 1

    print(f'Wrote {n_rows:,} rows to {output}', file=sys.stderr)


def gather_stats(files, output):
    L = []
    for f in files:
        chunk = f.replace('_filtered.stats', '').split('/')[-1]
        df = pd.read_csv(f, sep='\t').assign(chunk=chunk)
        L.append(df)
    out = pd.concat(L, ignore_index=True)
    out['n_records'] = out['n_records'].astype(int)
    out.to_csv(output, sep='\t', index=False)
    print(f'Wrote {len(out):,} rows to {output}', file=sys.stderr)


def gather_allelic(files, output):
    """Stream all per-(sample, shard) allelic TSVs into one gzip TSV."""
    open_fn = gzip.open if output.endswith('.gz') else open
    n_rows = 0
    writer = None

    with open_fn(output, 'wt') as out_fh:
        for f in files:
            in_open = gzip.open if f.endswith('.gz') else open
            with in_open(f, 'rt', newline='') as in_fh:
                reader = csv.DictReader(in_fh, delimiter='\t')
                if writer is None:
                    writer = csv.DictWriter(out_fh, fieldnames=reader.fieldnames,
                                            delimiter='\t')
                    writer.writeheader()
                for row in reader:
                    writer.writerow(row)
                    n_rows += 1

    print(f'Wrote {n_rows:,} rows to {output}', file=sys.stderr)


def main():
    args = parse_args()
    if args.mode == 'muts':
        gather_muts(args.files, args.output)
    elif args.mode == 'allelic':
        gather_allelic(args.files, args.output)
    else:
        gather_stats(args.files, args.output)


if __name__ == '__main__':
    main()
