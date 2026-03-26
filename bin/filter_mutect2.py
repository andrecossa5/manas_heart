#!/usr/bin/env python3

"""

Filter MUTECT2 unfiltered VCF for embryonic somatic mutations.
Samples in VCF: heart (col 10, tumor), placenta (col 11, normal).
Usage:
    python filter_mutect2.py <vcf> [--chunk NAME] [options]

"""

import sys
import os
import argparse
import pandas as pd
from cyvcf2 import VCF
from scipy.stats import fisher_exact


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('vcf', help='Input unfiltered VCF (.vcf.gz)')
    p.add_argument('--chunk', help='Chunk name (inferred from filename if omitted)')
    p.add_argument('--ad-placenta-max', type=int,   default=100,     metavar='N',
                   help='Max alt AD allowed in placenta [%(default)s]')
    p.add_argument('--ad-heart-min',    type=int,   default=3,     metavar='N',
                   help='Min alt AD in heart [%(default)s]')
    p.add_argument('--dp-min',          type=int,   default=10,    metavar='N',
                   help='Min DP required in both heart and placenta [%(default)s]')
    p.add_argument('--af-heart-max',    type=float, default=0.75,  metavar='F',
                   help='Max alt AF in heart [%(default)s]')
    p.add_argument('--mmq-min',         type=int,   default=20,    metavar='Q',
                   help='Min alt allele median mapping quality [%(default)s]')
    p.add_argument('--mbq-min',         type=int,   default=20,    metavar='Q',
                   help='Min alt allele median base quality [%(default)s]')
    p.add_argument('--sb-pval-min',     type=float, default=0.001, metavar='P',
                   help='Min Fisher strand bias p-value — discard if < this [%(default)s]')
    p.add_argument('--alt2-ad-ratio-max', type=float, default=0.2,  metavar='F',
                   help='Max ALT2/ALT1 AD ratio in heart for multiallelic records [%(default)s]')
    p.add_argument('--mpos-min',        type=int,   default=3,     metavar='N',
                   help='Min median distance from read end [%(default)s]')
    p.add_argument('--output', '-o', default=None,
                   help='Output TSV path (default: <chunk>_filtered.tsv)')
    return p.parse_args()


def infer_chunk(vcf_path):
    """Derive chunk name from filename, e.g. Myocardium.0_unfiltered.vcf.gz -> Myocardium.0."""
    name = os.path.basename(vcf_path)
    if '_unfiltered' in name:
        return name.split('_unfiltered')[0]
    return name.split('.vcf')[0]


def main():
    args = parse_args()
    chunk = args.chunk or infer_chunk(args.vcf)

    counters = dict(total=0, multiallelic=0, mmq=0, mbq=0, mpos=0,
                    ad_placenta=0, ad_heart=0, af_heart=0, dp=0, sb=0, passed=0)
    records = []
    vcf = VCF(args.vcf)

    for v in vcf:

        counters['total'] += 1

        # --- discard multiallelic only if ALT2 AD > 1/5 of ALT1 AD in heart ---
        if len(v.ALT) > 1:
            _ad = v.format('AD')
            if int(_ad[0, 2]) > int(_ad[0, 1]) * args.alt2_ad_ratio_max:
                counters['multiallelic'] += 1
                continue

        # --- INFO filters (cheap, no FORMAT parsing yet) ---
        mmq  = v.INFO.get('MMQ')   # Number=R: (ref_MMQ, alt_MMQ)
        mbq  = v.INFO.get('MBQ')   # Number=R: (ref_MBQ, alt_MBQ)
        mpos = v.INFO.get('MPOS')  # Number=A: scalar or (alt_MPOS,)

        if mmq is None or mbq is None or mpos is None:
            continue

        # Number=R → tuple of 2; Number=A with 1 alt → scalar or tuple of 1
        alt_mmq  = mmq[1]  if hasattr(mmq,  '__len__') else mmq
        alt_mbq  = mbq[1]  if hasattr(mbq,  '__len__') else mbq
        alt_mpos = mpos[0] if hasattr(mpos, '__len__') else mpos

        if alt_mmq < args.mmq_min:
            counters['mmq'] += 1
            continue
        if alt_mbq < args.mbq_min:
            counters['mbq'] += 1
            continue
        if alt_mpos < args.mpos_min:
            counters['mpos'] += 1
            continue

        # --- FORMAT filters ---
        # AD  Number=R → shape (n_samples, n_alleles): heart=0, placenta=1; ref=0, alt=1
        ad = v.format('AD')
        ad_heart_alt    = int(ad[0, 1])
        ad_placenta_alt = int(ad[1, 1])

        if ad_placenta_alt > args.ad_placenta_max:
            counters['ad_placenta'] += 1
            continue
        if ad_heart_alt <= args.ad_heart_min:
            counters['ad_heart'] += 1
            continue

        # AF  Number=A → shape (n_samples, 1) for biallelic
        af = v.format('AF')
        af_heart    = float(af[0, 0])
        af_placenta = float(af[1, 0])

        if af_heart >= args.af_heart_max:
            counters['af_heart'] += 1
            continue

        # DP  Number=1 → shape (n_samples, 1)
        dp = v.format('DP')
        dp_heart     = int(dp[0, 0])
        dp_placenta  = int(dp[1, 0])

        if dp_heart < args.dp_min or dp_placenta < args.dp_min:
            counters['dp'] += 1
            continue

        # --- Strand bias: Fisher's exact on heart FORMAT SB ---
        # SB  Number=4 → shape (n_samples, 4): ref_fwd, ref_rev, alt_fwd, alt_rev
        sb = v.format('SB')[0]
        _, sb_pval = fisher_exact([[int(sb[0]), int(sb[1])],
                                   [int(sb[2]), int(sb[3])]])

        if sb_pval < args.sb_pval_min:
            counters['sb'] += 1
            continue

        counters['passed'] += 1
        records.append({
            'chunk':        chunk,
            'CHROM':        v.CHROM,
            'POS':          v.POS,
            'REF':          v.REF,
            'ALT':          v.ALT[0],
            'AD_placenta':  ad_placenta_alt,
            'AD_heart':     ad_heart_alt,
            'DP_heart':     dp_heart,
            'DP_placenta':  dp_placenta,
            'AF_placenta':  af_placenta,
            'AF_heart':     af_heart,
            'median_BQ':    int(alt_mbq),
            'median_MQ':    int(alt_mmq),
            'SB_pval':      sb_pval,
            'MPOS':         int(alt_mpos),
        })

    vcf.close()

    # --- build dataframe ---
    cols = ['chunk', 'CHROM', 'POS', 'REF', 'ALT',
            'AD_placenta', 'AD_heart',
            'DP_heart', 'DP_placenta', 'AF_placenta', 'AF_heart',
            'median_MQ', 'median_BQ', 'SB_pval', 'MPOS']
    df = pd.DataFrame(records, columns=cols)

    out_tsv   = args.output or f'{chunk}_filtered.tsv'
    out_stats = out_tsv.replace('_filtered.tsv', '_filtered.stats')

    df.to_csv(out_tsv, sep='\t', index=False)

    # --- stats file ---
    stats = pd.DataFrame([
        ('total',          counters['total']),
        ('multiallelic',   counters['multiallelic']),
        ('failed_mmq',     counters['mmq']),
        ('failed_mbq',     counters['mbq']),
        ('failed_mpos',    counters['mpos']),
        ('failed_ad_placenta', counters['ad_placenta']),
        ('failed_ad_heart',    counters['ad_heart']),
        ('failed_af_heart',    counters['af_heart']),
        ('failed_dp',      counters['dp']),
        ('failed_sb',      counters['sb']),
        ('passed',         counters['passed']),
    ], columns=['filter', 'n_records'])
    stats.to_csv(out_stats, sep='\t', index=False)

    print(f'Wrote {len(df):,} records to {out_tsv}', file=sys.stderr)
    print(f'Wrote stats to {out_stats}', file=sys.stderr)


if __name__ == '__main__':
    main()


##