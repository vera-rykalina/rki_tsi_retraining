#!/usr/bin/env python3

import pandas as pd
import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Merge metadata with features table for retraining step 2, keep only means, viral_load, and known_tsi_years")
    parser.add_argument("--metadata", "-m", required=True, help="CSV file with metadata (includes 'scount', 'viral_load', 'known_tsi_years')")
    parser.add_argument("--features", "-f", required=True, help="CSV file with features (includes 'host.id' and mean features)")
    parser.add_argument("--output", "-o", required=True, help="Output CSV filename for merged and filtered table")
    return parser.parse_args()

def main():
    args = parse_args()

    # Load metadata
    meta_df = pd.read_csv(args.metadata)
    if 'scount' not in meta_df.columns:
        sys.exit("Error: 'scount' column not found in metadata file.")

    # Rename 'scount' to 'host.id'
    meta_df.rename(columns={'scount': 'host.id'}, inplace=True)

    if 'host.id' not in meta_df.columns:
        sys.exit("Error: 'host.id' column missing after renaming in metadata.")

    # Load features
    feat_df = pd.read_csv(args.features)
    if 'host.id' not in feat_df.columns:
        sys.exit("Error: 'host.id' column not found in features file.")

    # Merge on 'host.id'
    merged = pd.merge(meta_df, feat_df, on='host.id', how='inner')

    required_cols = ['host.id', 'viral_load', 'known_tsi_years']
    for col in required_cols:
        if col not in merged.columns:
            sys.exit(f"Error: '{col}' column missing in merged dataframe.")

    # Select feature columns from feat_df except 'host.id'
    feature_cols = [col for col in feat_df.columns if col != 'host.id']
    cols_to_keep = ['host.id', 'viral_load', 'known_tsi_years'] + feature_cols

    # Filter merged dataframe
    merged_filtered = merged[cols_to_keep]

    # Add 'sample_id' as the first column, copying 'host.id'
    merged_filtered.insert(0, 'sample_id', merged_filtered['host.id'])

    # Save output
    merged_filtered.to_csv(args.output, index=False)
    print(f"Merged and filtered table saved to {args.output}")

if __name__ == "__main__":
    main()
