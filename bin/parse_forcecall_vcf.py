#!/usr/bin/env python3

"""

Parse a forcecall (bcftools mpileup | bcftools call -C alleles) VCF into a
long-format allelic TSV.

For every site in --sites (CHROM POS REF ALT), emit one row with AD_ref,
AD_alt, DP and AF. Sites absent from the VCF (true zero coverage) are
emitted with AD_ref=0, AD_alt=0, DP=0, AF=NaN.

Usage:
    parse_forcecall_vcf.py --vcf <vcf.gz> --sites <shard.tsv> \\
        --sample-id <id> --tissue <heart|placenta> --chunk <chunk> \\
        --output <out.tsv.gz>

"""

import sys
import csv
import gzip
import argparse
from cyvcf2 import VCF


COLS = ['Sample_ID', 'tissue', 'chunk',
        'CHROM', 'POS', 'REF', 'ALT',
        'AD_ref', 'AD_alt', 'DP', 'AF']


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--vcf',       required=True, help='Per-sample forcecall VCF (.vcf.gz)')
    p.add_argument('--sites',     required=True, help='Shard sites TSV (no header): CHROM POS REF ALT')
    p.add_argument('--sample-id', required=True)
    p.add_argument('--tissue',    required=True)
    p.add_argument('--chunk',     required=True)
    p.add_argument('--output',    required=True, help='Output TSV (.tsv.gz)')
    return p.parse_args()


def load_sites(path):
    sites = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            chrom, pos, ref, alt = parts[0], int(parts[1]), parts[2], parts[3]
            sites.append((chrom, pos, ref, alt))
    return sites


def index_vcf(vcf_path):
    """Build a dict keyed by (CHROM, POS, REF, ALT) → (AD_ref, AD_alt, DP)."""
    idx = {}
    vcf = VCF(vcf_path)
    for v in vcf:
        if v.ALT is None:
            continue
        # AD FORMAT shape: (1, n_alleles); DP FORMAT shape: (1, 1)
        ad = v.format('AD')
        dp_fmt = v.format('DP')
        dp = int(dp_fmt[0, 0]) if dp_fmt is not None else 0
        ad_ref = int(ad[0, 0]) if ad is not None else 0

        alts = v.ALT if isinstance(v.ALT, list) else [v.ALT]
        for i, alt in enumerate(alts):
            ad_alt = int(ad[0, i + 1]) if ad is not None and ad.shape[1] > i + 1 else 0
            idx[(v.CHROM, v.POS, v.REF, alt)] = (ad_ref, ad_alt, dp)

    vcf.close()
    return idx


def main():
    args = parse_args()

    sites = load_sites(args.sites)
    idx   = index_vcf(args.vcf)

    n_hit = 0
    with gzip.open(args.output, 'wt', newline='') as out_fh:
        w = csv.DictWriter(out_fh, fieldnames=COLS, delimiter='\t')
        w.writeheader()

        for chrom, pos, ref, alt in sites:
            rec = idx.get((chrom, pos, ref, alt))
            if rec is not None:
                ad_ref, ad_alt, dp = rec
                n_hit += 1
            else:
                ad_ref, ad_alt, dp = 0, 0, 0

            af = (ad_alt / dp) if dp > 0 else ''

            w.writerow({
                'Sample_ID': args.sample_id,
                'tissue':    args.tissue,
                'chunk':     args.chunk,
                'CHROM':     chrom,
                'POS':       pos,
                'REF':       ref,
                'ALT':       alt,
                'AD_ref':    ad_ref,
                'AD_alt':    ad_alt,
                'DP':        dp,
                'AF':        af,
            })

    print(f'Wrote {len(sites):,} rows ({n_hit:,} with coverage) to {args.output}',
          file=sys.stderr)


if __name__ == '__main__':
    main()
