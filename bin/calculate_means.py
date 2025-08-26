#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import warnings

warnings.filterwarnings("ignore")

def loginfo(s):
    print(f"[INFO] {s}")

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate per-sample feature means from patstats and MAF.")
    parser.add_argument("--patstats", "-p", required=True, help="Path to masked patstats CSV")
    parser.add_argument("--maf", "-m", required=True, help="Path to MAF CSV")
    parser.add_argument("--ref", required=True, help="Path to HXB2 reference data CSV")
    parser.add_argument("--output", "-o", required=True, help="Filename for saving output feature table (CSV)")
    return parser.parse_args()

def load_reference_data(ref_path):
    hxb2 = pd.read_csv(ref_path)
    hxb2['position'] = hxb2['HXB2 base position']
    hxb2.set_index('position', inplace=True)

    rf3_3cp = hxb2.groupby(['RF3 protein', 'RF3 aa position'])['HXB2 base position'].max()
    rf2_3cp = hxb2.groupby(['RF2 protein', 'RF2 aa position'])['HXB2 base position'].max()
    rf1_3cp = hxb2.groupby(['RF1 protein', 'RF1 aa position'])['HXB2 base position'].max()

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

def load_patstats(fpath):
    df = pd.read_csv(fpath)
    # Pivot to have host.id as index, xcoord as columns, values are features
    Xlrtt = df.pivot(index='host.id', columns='xcoord', values='normalised.largest.rtt')
    Xtips = df.pivot(index='host.id', columns='xcoord', values='tips')
    Xdual = df.pivot(index='host.id', columns='xcoord', values='solo.dual.count')
    return Xlrtt, Xtips, Xdual

def load_maf(fpath):
    maf = pd.read_csv(fpath, index_col=0)
    maf.columns = [int(float(c)) for c in maf.columns]
    return maf

def calculate_means(Xlrtt, Xtips, Xdual, maf, first_second, third, genes):
    df = pd.DataFrame(index=Xlrtt.index)

    # Whole genome means
    df['genome_lrtt'] = Xlrtt.mean(axis=1)
    df['genome_tips'] = Xtips.mean(axis=1)
    df['genome_dual'] = Xdual.mean(axis=1)
    df['genome_maf12c'] = maf.loc[:, maf.columns.intersection(first_second)].mean(axis=1)
    df['genome_maf3c'] = maf.loc[:, maf.columns.intersection(third)].mean(axis=1)

    # Per gene means
    for gene in ['gag', 'pol', 'gp120', 'gp41']:
        gene_pos = maf.columns.intersection(genes[gene])
        df[f'{gene}_lrtt'] = Xlrtt.loc[:, Xlrtt.columns.intersection(gene_pos)].mean(axis=1)
        df[f'{gene}_tips'] = Xtips.loc[:, Xtips.columns.intersection(gene_pos)].mean(axis=1)
        df[f'{gene}_dual'] = Xdual.loc[:, Xdual.columns.intersection(gene_pos)].mean(axis=1)
        df[f'{gene}_maf12c'] = maf.loc[:, gene_pos.intersection(first_second)].mean(axis=1)
        df[f'{gene}_maf3c'] = maf.loc[:, gene_pos.intersection(third)].mean(axis=1)

    return df

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

if __name__ == "__main__":
    main()