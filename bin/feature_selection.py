#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import KNNImputer
from itertools import combinations

def knn_impute(df, n_neighbors=3):
    """Impute missing data using KNN."""
    if df.isna().sum().sum() == 0:
        print("No missing data to impute, continuing.")
        return df
    imputer = KNNImputer(n_neighbors=n_neighbors)
    imputed_array = imputer.fit_transform(df)
    return pd.DataFrame(imputed_array, columns=df.columns, index=df.index)

def evaluate_combos(df, features, target_col, include_vl=False):
    """Evaluate feature combinations and return results."""
    results = []
    # If include viral_load and viral_load not already in features, add it
    if include_vl and 'viral_load' not in features and 'viral_load' in df.columns:
        features = features + ['viral_load']

    # Limit max combination size for speed (e.g., max 5 features)
    max_feats = min(5, len(features))
    
    for size in range(1, max_feats + 1):
        for combo in combinations(features, size):
            X_sub = df[list(combo)]
            X_imputed = knn_impute(X_sub)
            y = df[target_col].values
            rf = RandomForestRegressor(n_estimators=100, oob_score=True, random_state=42)
            try:
                rf.fit(X_imputed, y)
            except ValueError as e:
                print(f"Skipping combo {combo} due to error: {e}")
                continue
            oob_score = rf.oob_score_
            results.append({
                'features': combo,
                'num_features': len(combo),
                'oob_score': oob_score
            })
    return results

def save_csv_report(results, output_csv):
    """Save top 10 combos with viral_load included as CSV."""
    # Sort descending by oob_score and keep top 10
    sorted_res = sorted(results, key=lambda x: x['oob_score'], reverse=True)[:10]
    rows = []
    for r in sorted_res:
        rows.append({
            'num_features': r['num_features'],
            'oob_score': r['oob_score'],
            'feature_names': ','.join(r['features'])
        })
    df_out = pd.DataFrame(rows)
    df_out.to_csv(output_csv, index=False)

def save_txt_report(results, output_txt):
    """Save best combo WITHOUT viral_load as txt (only features names)."""
    # Filter combos without viral_load
    filtered = [r for r in results if 'viral_load' not in r['features']]
    if not filtered:
        print("No feature combos found without viral_load for TXT report.")
        return
    # Pick best by oob_score
    best = max(filtered, key=lambda x: x['oob_score'])
    with open(output_txt, 'w') as f:
        f.write(f"Best feature combo WITHOUT viral_load (OOB score: {best['oob_score']:.4f}):\n")
        f.write('\n'.join(best['features']))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Input CSV file')
    parser.add_argument('--output-csv', required=True, help='Output CSV file for top combos with viral_load')
    parser.add_argument('--output-txt', required=True, help='Output TXT file for best combo without viral_load')
    args = parser.parse_args()

    df = pd.read_csv(args.input)

    target_col = 'known_tsi_years'
    if target_col not in df.columns:
        raise ValueError(f"Input CSV must contain '{target_col}' column")

    # Exclude columns with these prefixes because they have no values or are irrelevant
    exclude_prefixes = ['gp41_', 'gp120_', 'genome_']

    # Exclude identifiers and target
    exclude_cols = ['sample_id', 'host.id', target_col]
    
    # Filter features, exclude prefixes
    features_all = [c for c in df.columns 
                    if c not in exclude_cols 
                    and not any(c.startswith(pref) for pref in exclude_prefixes)]

    # Separate features with viral_load and without viral_load
    features_without_vl = [f for f in features_all if f != 'viral_load']

    # Evaluate combos without viral_load for TXT report
    results_without_vl = evaluate_combos(df, features_without_vl, target_col, include_vl=False)
    save_txt_report(results_without_vl, args.output_txt)

    # Evaluate combos with viral_load for CSV report
    results_with_vl = evaluate_combos(df, features_without_vl, target_col, include_vl=True)
    save_csv_report(results_with_vl, args.output_csv)

if __name__ == '__main__':
    main()