#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import re

def parse_coords(coord_str):
    start, end = map(int, coord_str.split('-'))
    return start, end

def find_non_nan_ranges(valid_cols):
    ranges = []
    if not valid_cols:
        return ranges

    sorted_cols = sorted(valid_cols)
    start = prev = sorted_cols[0]

    for col in sorted_cols[1:]:
        if col == prev + 1:
            prev = col
        else:
            ranges.append((start, prev))
            start = prev = col
    ranges.append((start, prev))
    return ranges

def natural_sort_key(s):
    # Split string into numeric and non-numeric parts for natural sort
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', str(s))]

def main():
    parser = argparse.ArgumentParser(description="Mask primer and non-product positions in a table and print valid data regions.")
    parser.add_argument("-r", "--regions", required=True, help="CSV file with primer and product coords")
    parser.add_argument("-m", "--maf", required=True, help="Input CSV table to mask")
    parser.add_argument("-o", "--output", required=True, help="Output masked CSV file")
    args = parser.parse_args()

    # Load main data
    df = pd.read_csv(args.maf, index_col=0)
    positions = df.columns.astype(int)
    max_pos = positions.max()

    # Load regions table
    regions_df = pd.read_csv(args.regions)
    regions_df.columns = regions_df.columns.str.strip()

    primer_positions = set()
    product_positions = set()

    for _, row in regions_df.iterrows():
        sp_start, sp_end = parse_coords(row["sense_primer_coords"])
        asp_start, asp_end = parse_coords(row["antisense_primer_coords"])
        pr_start, pr_end = parse_coords(row["product_coords_incl_primers"])

        primer_positions.update(range(sp_start, sp_end + 1))
        primer_positions.update(range(asp_start, asp_end + 1))
        product_positions.update(range(pr_start, pr_end + 1))

    # Apply masking
    primer_mask = positions.isin(primer_positions)
    df.loc[:, primer_mask] = np.nan

    keep_mask = positions.isin(product_positions)
    df.loc[:, ~keep_mask] = np.nan

    # Sort index naturally but preserve index
    df = df.loc[sorted(df.index, key=natural_sort_key)]

    # Save masked table, keeping index
    df.to_csv(args.output)

    # Identify columns with valid (non-NaN) values
    non_nan_cols = df.columns[df.notna().any(axis=0)].astype(int).tolist()
    valid_ranges = find_non_nan_ranges(non_nan_cols)

    print("\nValid regions (non-NaN columns after masking):\n")
    for start, end in valid_ranges:
        print(f"  Valid region: {start}-{end}")

if __name__ == "__main__":
    main()