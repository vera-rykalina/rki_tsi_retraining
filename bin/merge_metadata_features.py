#!/usr/bin/env python3

import pandas as pd
import argparse

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

    # Load features
    feat_df = pd.read_csv(args.features)

    # Merge on scount == host.id
    merged = pd.merge(meta_df, feat_df, left_on='scount', right_on='host.id', how='inner')

    # Select columns:
    # - IDs: 'scount' (to rename to sample_id) and 'host.id'
    # - metadata cols viral_load and known_tsi_years
    # - all feature columns except 'host.id'
    feature_cols = [col for col in feat_df.columns if col != 'host.id']

    cols_to_keep = ['scount', 'host.id', 'viral_load', 'known_tsi_years'] + feature_cols

    # Filter merged dataframe
    merged_filtered = merged.loc[:, cols_to_keep]

    # Rename scount to sample_id
    merged_filtered.rename(columns={'scount': 'sample_id'}, inplace=True)

    # Save
    merged_filtered.to_csv(args.output, index=False)
    print(f"Merged and filtered table saved to {args.output}")

if __name__ == "__main__":
    main()