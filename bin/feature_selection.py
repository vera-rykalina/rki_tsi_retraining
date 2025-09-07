#!/usr/bin/env python3

import argparse, os
import pandas as pd
import numpy as np
from itertools import combinations
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import KNNImputer
from sklearn.model_selection import GridSearchCV, train_test_split

def loginfo(msg):
    print(msg)

def _impute_knn(X, k=3):
    if np.isnan(X).sum() == 0:
        loginfo("No missing data; skipping imputation.")
        return X
    loginfo(f"Imputing {np.isnan(X).sum()} missing values with {k}-NN.")
    return KNNImputer(n_neighbors=k).fit_transform(X)

def evaluate_combos(df, base_features, include_flags=None, target_col='known_tsi_years'):
    results = []
    features = base_features.copy()
    if include_flags:
        for flag in include_flags:
            if flag in df.columns and flag not in features:
                features.append(flag)

    # Impute once on all features combined
    X_full_imputed = _impute_knn(df[features].values)

    y_full = np.sqrt(df[target_col].values)  # Keep target in sqrt space

    for size in range(1, min(8, len(features)) + 1):
        for combo in combinations(features, size):
            idxs = [features.index(f) for f in combo]
            X_full = X_full_imputed[:, idxs]

            if X_full.shape[1] == 0:
                loginfo(f"Skipping combo: {combo} (empty after imputation)")
                continue

            # Split once per combo into train/validation sets (in sqrt space)
            X_train, X_valid, y_train, y_valid = train_test_split(
                X_full, y_full, test_size=0.2, random_state=42
            )

            rf = RandomForestRegressor(random_state=42, n_jobs=8)
            param_grid = {
                'n_estimators': [100, 200],
                'max_depth': [5, 10, 20, None],
                'min_samples_split': [2, 5, 10]
            }

            grid = GridSearchCV(rf, param_grid, cv=10, scoring='neg_mean_absolute_error', n_jobs=8)
            try:
                grid.fit(X_train, y_train)
            except Exception as e:
                loginfo(f"Error fitting combo {combo}: {e}")
                continue

            best_params = grid.best_params_

            # Retrain best model on train with OOB enabled
            rf_best = RandomForestRegressor(
                random_state=42,
                n_jobs=8,
                oob_score=True,
                **best_params
            )
            rf_best.fit(X_train, y_train)

            # Evaluate validation MAE in original linear scale by squaring predictions and targets
            val_pred = rf_best.predict(X_valid)
            val_mae = np.mean(np.abs((y_valid**2) - (val_pred**2)))  # MAE on linear scale

            results.append({
                'features': combo,
                'n_estimators': best_params['n_estimators'],
                'oob_score': rf_best.oob_score_,
                'val_mae': val_mae,
                'best_params': best_params
            })

    return results

def save_csv(results, outdir, label):
    os.makedirs(outdir, exist_ok=True)
    df_out = pd.DataFrame({
        'num_features': [len(r['features']) for r in results],
        'n_estimators': [r['n_estimators'] for r in results],
        'oob_score': [r['oob_score'] for r in results],
        'val_mae': [r['val_mae'] for r in results],
        'features': [','.join(r['features']) for r in results],
        'best_params': [r['best_params'] for r in results]
    })
    # Sort by highest OOB score, then lowest val_mae
    df_top = df_out.sort_values(['oob_score', 'val_mae'], ascending=[False, True]).head(10)
    path = os.path.join(outdir, f"top10_{label}.csv")
    df_top.to_csv(path, index=False)
    return path

def save_txt_best(results, outdir, filename):
    if not results:
        loginfo(f"No combos available for TXT output: {filename}")
        return None
    # Pick best combo by OOB score first, then val_mae if tie
    best = max(results, key=lambda r: (r['oob_score'], -r['val_mae']))
    path = os.path.join(outdir, filename)
    with open(path, 'w') as f:
        for feat in best['features']:
            f.write(f"{feat}\n")
    return path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('--amplicons', action='store_true')  # Still accepted, but no longer used
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    if 'known_tsi_years' not in df:
        raise ValueError("Missing 'known_tsi_years' column.")

    if 'viral_load' in df.columns:
        df['viral_load'] = np.log10(df['viral_load'] + 1)

    exclude_cols = {'sample_id', 'host.id', 'known_tsi_years', 'viral_load', 'is_mrc'}
    # Removed 'dual_' from blacklist to include those features
    blacklist = ['genome_', 'gp120_', 'gp41_']
    base_features = [
        c for c in df.columns
        if c not in exclude_cols and not any(c.startswith(pref) for pref in blacklist)
    ]

    outdir = args.output

    # 1. Base only
    res_base = evaluate_combos(df, base_features)
    csv_base = save_csv(res_base, outdir, 'base_only')
    txt_base = save_txt_best(res_base, outdir, 'best_features_base_only.txt')

    # 2. Base + viral_load
    res_vl = evaluate_combos(df, base_features, include_flags=['viral_load'])
    csv_vl = save_csv(res_vl, outdir, 'with_viral_load')

    loginfo("Saved files:")
    for path in [csv_base, txt_base, csv_vl]:
        if path:
            loginfo(f" • {path}")

if __name__ == '__main__':
    main()
