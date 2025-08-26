#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import pandas.api.types as ptypes

def parse_coords(coord_str):
    """Parse coordinate string of format 'start-end' into a tuple of integers."""
    start, end = map(int, coord_str.split('-'))
    return start, end

def main():
    parser = argparse.ArgumentParser(description="Mask patstats rows based on product and primer regions.")
    parser.add_argument("-r", "--regions", required=True, help="CSV with primer and product coordinates")
    parser.add_argument("-p", "--patstats", required=True, help="Input patstats CSV file")
    parser.add_argument("-o", "--output", required=True, help="Output masked patstats CSV file")
    args = parser.parse_args()

    # Load patstats
    df = pd.read_csv(args.patstats)
    if "xcoord" not in df.columns:
        raise ValueError("Input patstats file must contain 'xcoord' column.")
    df['xcoord'] = df['xcoord'].astype(int)

    # Load regions table
    regions_df = pd.read_csv(args.regions)
    regions_df.columns = regions_df.columns.str.strip()  # remove any accidental whitespace

    # Collect all primer and product positions
    primer_positions = set()
    product_positions = set()

    for _, row in regions_df.iterrows():
        sp_start, sp_end = parse_coords(row["sense_primer_coords"])
        asp_start, asp_end = parse_coords(row["antisense_primer_coords"])
        pr_start, pr_end = parse_coords(row["product_coords_incl_primers"])

        primer_positions.update(range(sp_start, sp_end + 1))
        primer_positions.update(range(asp_start, asp_end + 1))
        product_positions.update(range(pr_start, pr_end + 1))

    # Function to determine if a 250bp window is valid (fully within product, no primer overlap)
    def is_valid_window(xcoord):
        window_start = xcoord - 125
        window_end = xcoord + 124
        window_range = set(range(window_start, window_end + 1))

        # Must be fully inside product region
        if not window_range.issubset(product_positions):
            return False

        # Must not overlap any primer region
        if window_range.intersection(primer_positions):
            return False

        return True

    # Apply the window check
    df['keep'] = df['xcoord'].apply(is_valid_window)

    # Mask numeric columns except host.id, tree.id, xcoord
    exclude_cols = ['host.id', 'tree.id', 'xcoord', 'keep']

    numeric_cols = [col for col in df.columns if col not in exclude_cols and ptypes.is_numeric_dtype(df[col])]
    
    df.loc[~df['keep'], numeric_cols] = np.nan

    # Drop helper column
    df.drop(columns=['keep'], inplace=True)

    # Sort for easier inspection and stable output
    df.sort_values(by=["host.id", "xcoord"], inplace=True)

    # Save masked data
    df.to_csv(args.output, index=False)
    print(f"\n✅ Masked patstats saved to: {args.output}")

if __name__ == "__main__":
    main()