#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import KNNImputer
from itertools import combinations

def loginfo(msg):
    print(msg)

def _impute_knn(X, k=3):
    ''' Use k nearest rows which have a feature to fill in each row's missing features. '''
    if np.isnan(X).sum() == 0:
        loginfo('No missing data to impute, continuing.')
        return X

    nm = (np.isnan(X).sum(axis=1) != 0).sum()
    loginfo(f'Using KNN ({k} nearest neighbours) to impute {np.isnan(X).sum()} missing values ({nm} sample{"s" if nm > 1 else ""})')

    # Show where missing values are (row and column indices)
    rows, cols = np.where(np.isnan(X))
    unique_rows = np.unique(rows)
    unique_cols = np.unique(cols)
    loginfo(f'Missing data found at rows: {unique_rows.tolist()} and columns indexes: {unique_cols.tolist()}')

    imputer = KNNImputer(n_neighbors=k)
    X_knn = imputer.fit_transform(X)
    return X_knn

def evaluate_combos(df, features, target_col, include_vl=False, n_estimators_list=[100, 200, 300]):
    results = []

    if include_vl and 'viral_load' not in features and 'viral_load' in df.columns:
        features = features + ['viral_load']

    max_feats = min(5, len(features))
    mandatory_feats = ['is_mrc']

    if 'is_mrc' not in df.columns:
        raise ValueError("'is_mrc' column must be present in dataframe before evaluating")

    optional_feats = [f for f in features if f not in mandatory_feats]

    for size in range(1, max_feats + 1):
        for combo in combinations(optional_feats, size):
            combo_with_mand = tuple(mandatory_feats + list(combo))
            X_sub = df[list(combo_with_mand)].values

            X_imputed = _impute_knn(X_sub)
            y = np.sqrt(df[target_col].values)

            for n_estimators in n_estimators_list:
                rf = RandomForestRegressor(n_estimators=n_estimators, oob_score=True, random_state=42)
                try:
                    rf.fit(X_imputed, y)
                except ValueError as e:
                    loginfo(f"Skipping combo {combo_with_mand} with n_estimators={n_estimators} due to error: {e}")
                    continue

                # OOB prediction and outlier detection
                oob_preds = rf.oob_prediction_
                errors = np.abs(y - oob_preds)

                # Detecting outliers based on the errors
                top_n = 5
                outlier_indices = np.argsort(errors)[-top_n:][::-1]
                outlier_info = [(idx, errors[idx]) for idx in outlier_indices]

                # Track outlier categories: mild, moderate, extreme
                mild_threshold = np.percentile(errors, 75)
                moderate_threshold = np.percentile(errors, 90)
                extreme_threshold = np.percentile(errors, 95)

                outlier_categories = {
                    'mild': [(idx, err) for idx, err in zip(outlier_indices, errors[outlier_indices]) if err < mild_threshold],
                    'moderate': [(idx, err) for idx, err in zip(outlier_indices, errors[outlier_indices]) if mild_threshold <= err < moderate_threshold],
                    'extreme': [(idx, err) for idx, err in zip(outlier_indices, errors[outlier_indices]) if err >= extreme_threshold]
                }

                loginfo(f"Outliers for combo {combo_with_mand} with n_estimators={n_estimators}:")
                loginfo(f"Mild: {outlier_categories['mild']}")
                loginfo(f"Moderate: {outlier_categories['moderate']}")
                loginfo(f"Extreme: {outlier_categories['extreme']}")

                # Feature importance
                feature_importance = rf.feature_importances_

                results.append({
                    'features': combo_with_mand,
                    'num_features': len(combo_with_mand),
                    'n_estimators': n_estimators,
                    'oob_score': rf.oob_score_,
                    'outliers': outlier_info,
                    'feature_importance': feature_importance.tolist()
                })
    return results

def save_csv_report(results, output_csv):
    """Save top 10 combos sorted by oob_score descending to CSV."""
    sorted_res = sorted(results, key=lambda x: x['oob_score'], reverse=True)[:10]
    rows = []
    for r in sorted_res:
        rows.append({
            'num_features': r['num_features'],
            'n_estimators': r['n_estimators'],
            'oob_score': r['oob_score'],
            'feature_names': ','.join(r['features']),
            'feature_importance': ','.join(map(str, r['feature_importance']))
        })
    df_out = pd.DataFrame(rows)
    df_out.to_csv(output_csv, index=False)

def save_txt_report(results, output_txt):
    """Save best combo WITHOUT viral_load as txt (only feature names, no header)."""
    filtered = [r for r in results if 'viral_load' not in r['features']]
    if not filtered:
        loginfo("No feature combos found without viral_load for TXT report.")
        return
    best = max(filtered, key=lambda x: x['oob_score'])
    with open(output_txt, 'w') as f:
        f.write('\n'.join(best['features']))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Input CSV file')
    parser.add_argument('--output-csv-with-vl', required=True, help='Output CSV file for top combos WITH viral_load')
    parser.add_argument('--output-csv-without-vl', required=True, help='Output CSV file for top combos WITHOUT viral_load')
    parser.add_argument('--output-txt', required=True, help='Output TXT file for best combo without viral_load')
    parser.add_argument('--amplicons', action='store_true', help='Set is_mrc=1 if data are amplicons, else 0')
    args = parser.parse_args()

    df = pd.read_csv(args.input)

    # Preprocess viral_load to log10(viral_load + 1)
    if 'viral_load' in df.columns:
        df['viral_load'] = np.log10(df['viral_load'] + 1)

    target_col = 'known_tsi_years'
    if target_col not in df.columns:
        raise ValueError(f"Input CSV must contain '{target_col}' column")

    exclude_prefixes = ['gp41_', 'gp120_', 'genome_']
    exclude_cols = ['sample_id', 'host.id', target_col]

    features_all = [c for c in df.columns
                    if c not in exclude_cols
                    and not any(c.startswith(pref) for pref in exclude_prefixes)]

    df['is_mrc'] = int(args.amplicons)

    features_without_vl = [f for f in features_all if f != 'viral_load']

    results_without_vl = evaluate_combos(df, features_without_vl, target_col, include_vl=False,
                                        n_estimators_list=[100, 200, 300])
    save_csv_report(results_without_vl, args.output_csv_without_vl)
    save_txt_report(results_without_vl, args.output_txt)

    results_with_vl = evaluate_combos(df, features_without_vl, target_col, include_vl=True,
                                     n_estimators_list=[100, 200, 300])
    save_csv_report(results_with_vl, args.output_csv_with_vl)

if __name__ == '__main__':
    main()