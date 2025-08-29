#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import warnings

warnings.filterwarnings("ignore")

# ----------------------------------------------------------------------
# Logging
# ----------------------------------------------------------------------
def loginfo(s):
    print(f"[INFO] {s}")

# ----------------------------------------------------------------------
# Argument Parsing
# ----------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Calculate per-sample feature means from patstats and MAF.")
    parser.add_argument("--patstats", "-p", required=True, help="Path to masked patstats CSV")
    parser.add_argument("--maf", "-m", required=True, help="Path to MAF CSV")
    parser.add_argument("--ref", required=True, help="Path to HXB2 reference data CSV")
    parser.add_argument("--output", "-o", required=True, help="Filename for saving output feature table (CSV)")
    return parser.parse_args()

# ----------------------------------------------------------------------
# Load Reference Data (HXB2)
# ----------------------------------------------------------------------
def load_reference_data(ref_path):
    hxb2 = pd.read_csv(ref_path)
    hxb2['position'] = hxb2['HXB2 base position']
    hxb2.set_index('position', inplace=True)

    # Get 3rd codon positions from RF1, RF2, RF3
    rf3_3cp = hxb2.groupby(['RF3 protein', 'RF3 aa position'])['HXB2 base position'].max()
    rf2_3cp = hxb2.groupby(['RF2 protein', 'RF2 aa position'])['HXB2 base position'].max()
    rf1_3cp = hxb2.groupby(['RF1 protein', 'RF1 aa position'])['HXB2 base position'].max()

    # Get 1st + 2nd codon positions
    rf1_12 = set(hxb2.groupby(['RF1 protein', 'RF1 aa position'])['HXB2 base position'].nsmallest(2).reset_index()['HXB2 base position'])
    rf2_12 = set(hxb2.groupby(['RF2 protein', 'RF2 aa position'])['HXB2 base position'].nsmallest(2).reset_index()['HXB2 base position'])
    rf3_12 = set(hxb2.groupby(['RF3 protein', 'RF3 aa position'])['HXB2 base position'].nsmallest(2).reset_index()['HXB2 base position'])

    first_second_codon_pos = rf1_12 | rf2_12 | rf3_12
    third_codon_pos = (set(rf1_3cp) | set(rf2_3cp) | set(rf3_3cp)) - first_second_codon_pos

    genes = {
        'gag': set(hxb2[hxb2['RF1 protein'] == 'gag'].index),
        'pol': set(hxb2[hxb2['RF3 protein'] == 'pol'].index),
        'gp120': set(hxb2[hxb2['RF3 protein'] == 'gp120'].index),
        'gp41': set(hxb2[hxb2['RF3 protein'] == 'gp41'].index),
    }

    return first_second_codon_pos, third_codon_pos, genes

# ----------------------------------------------------------------------
# Load PatStats
# ----------------------------------------------------------------------
def load_patstats(fpath):
    df = pd.read_csv(fpath)
    # Aggregate by host.id and xcoord, then unstack
    Xlrtt = df.groupby(['host.id', 'xcoord'])['normalised.largest.rtt'].mean().unstack()
    Xtips = df.groupby(['host.id', 'xcoord'])['tips'].mean().unstack()
    Xdual = df.groupby(['host.id', 'xcoord'])['solo.dual.count'].mean().unstack()
    return Xlrtt, Xtips, Xdual

# ----------------------------------------------------------------------
# Load MAF
# ----------------------------------------------------------------------
def load_maf(fpath):
    maf = pd.read_csv(fpath, index_col=0)
    maf.columns = [int(float(c)) for c in maf.columns]
    return maf

# ----------------------------------------------------------------------
# Calculate Feature Means (Refactored)
# ----------------------------------------------------------------------
def calculate_means(Xlrtt, Xtips, Xdual, maf, first_second, third, genes):
    df = pd.DataFrame(index=Xlrtt.index)

    # Genome-wide features
    df['genome_lrtt'] = Xlrtt.mean(axis=1)
    df['genome_tips'] = Xtips.mean(axis=1)
    df['genome_dual'] = Xdual.mean(axis=1)
    df['genome_maf12c'] = maf.loc[:, maf.columns.intersection(first_second)].mean(axis=1)
    df['genome_maf3c'] = maf.loc[:, maf.columns.intersection(third)].mean(axis=1)

    # Per-gene features
    for gene, positions in genes.items():
        gene_pos = maf.columns.intersection(positions)
        maf12 = gene_pos.intersection(first_second)
        maf3 = gene_pos.intersection(third)

        df[f'{gene}_lrtt'] = Xlrtt.loc[:, Xlrtt.columns.intersection(gene_pos)].mean(axis=1)
        df[f'{gene}_tips'] = Xtips.loc[:, Xtips.columns.intersection(gene_pos)].mean(axis=1)
        df[f'{gene}_dual'] = Xdual.loc[:, Xdual.columns.intersection(gene_pos)].mean(axis=1)
        df[f'{gene}_maf12c'] = maf.loc[:, maf12].mean(axis=1)
        df[f'{gene}_maf3c'] = maf.loc[:, maf3].mean(axis=1)

    return df

# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    args = parse_args()

    loginfo("Loading reference data...")
    first_second, third, genes = load_reference_data(args.ref)

    loginfo("Loading patstats...")
    Xlrtt, Xtips, Xdual = load_patstats(args.patstats)

    loginfo("Loading MAF...")
    maf = load_maf(args.maf)

    loginfo("Calculating means...")
    df_means = calculate_means(Xlrtt, Xtips, Xdual, maf, first_second, third, genes)

    loginfo(f"Saving output to {args.output}")
    df_means.to_csv(args.output)

# ----------------------------------------------------------------------
# Entry Point
# ----------------------------------------------------------------------
if __name__ == "__main__":
    main()