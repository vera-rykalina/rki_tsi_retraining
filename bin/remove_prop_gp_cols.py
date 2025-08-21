#!/usr/bin/env python3

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Remove columns starting with 'prop.gp.' from CSV")
    parser.add_argument('--input', '-i', required=True, help="Input CSV file")
    parser.add_argument('--output', '-o', required=True, help="Output CSV file")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    # Select columns that do NOT start with 'prop.gp.'
    cols_to_keep = [col for col in df.columns if not col.startswith('prop.gp.')]
    df_filtered = df[cols_to_keep]

    df_filtered.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()