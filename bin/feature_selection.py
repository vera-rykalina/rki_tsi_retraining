#!/usr/bin/env python3

import argparse, os
import pandas as pd
import numpy as np
from itertools import combinations
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import KNNImputer
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.metrics import r2_score

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

    # Impute once for efficiency
    X_full_imputed = _impute_knn(df[features].values)
    y_full = np.sqrt(df[target_col].values)

    max_features = min(8, len(features))

    for size in range(1, max_features + 1):  # Up to 8 features max
        for combo in combinations(features, size):
            # Ensure is_mrc is always in the first position if it exists
            if 'is_mrc' in features:
                if 'is_mrc' not in combo:
                    # skip combos that do not include is_mrc
                    continue
                # reorder so that is_mrc is first
                combo = ('is_mrc',) + tuple(f for f in combo if f != 'is_mrc')

            idxs = [features.index(f) for f in combo]
            X_full = X_full_imputed[:, idxs]

            if X_full.shape[1] == 0:
                loginfo(f"Skipping combo: {combo} (empty after imputation)")
                continue

            # 75% train, 25% test split
            X_train, X_valid, y_train, y_valid = train_test_split(
                X_full, y_full, test_size=0.25, random_state=42
            )

            rf = RandomForestRegressor(random_state=42, n_jobs=-1)
            param_grid = {
               'n_estimators': [100, 200, 1000],
               'max_depth': [3, 5, 7],
               'min_samples_split': [2, 5, 10],
               'min_samples_leaf': [1, 2, 4]
            }

            grid = GridSearchCV(rf, param_grid, cv=10, scoring='neg_mean_absolute_error', n_jobs=-1)
            try:
                grid.fit(X_train, y_train)
            except Exception as e:
                loginfo(f"Error fitting combo {combo}: {e}")
                continue

            best_params = grid.best_params_

            rf_best = RandomForestRegressor(
                random_state=42,
                n_jobs=-1,
                oob_score=True,
                **best_params
            )
            rf_best.fit(X_train, y_train)

            val_pred = rf_best.predict(X_valid)
            val_mae = np.mean(np.abs((y_valid**2) - (val_pred**2)))  # MAE in original scale
            val_r2 = r2_score(y_valid**2, val_pred**2)               # R2 in original scale

            results.append({
                'features': combo,
                'n_estimators': best_params['n_estimators'],
                'oob_score': rf_best.oob_score_,
                'val_mae': val_mae,
                'val_r2': val_r2,
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
        'val_r2': [r['val_r2'] for r in results],
        'features': [','.join(r['features']) for r in results],
        'best_params': [r['best_params'] for r in results]
    })

    df_top = df_out.sort_values(['oob_score', 'val_mae'], ascending=[False, True]).head(10)
    path = os.path.join(outdir, f"top10_{label}.csv")
    df_top.to_csv(path, index=False)
    return path

def save_txt_best(results, outdir, filename):
    if not results:
        loginfo(f"No combos available for TXT output: {filename}")
        return None
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
    parser.add_argument('--amplicons', action='store_true', help="Set is_mrc=1 (for compatibility, not used)")
    parser.add_argument('--include_viral_load', action='store_true', help="Include viral_load as a feature")
    parser.add_argument('--include_is_mrc', action='store_true', help="Include is_mrc as the first feature")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    if 'known_tsi_years' not in df:
        raise ValueError("Missing 'known_tsi_years' column.")

    # Add is_mrc column dynamically based on --amplicons flag
    if args.amplicons:
        df['is_mrc'] = 1
    else:
        df['is_mrc'] = 0


    exclude_cols = {'sample_id', 'host.id', 'known_tsi_years', 'viral_load', 'is_mrc'}
    blacklist = ['genome_', 'gp120_', 'gp41_']

    base_features = [
        c for c in df.columns
        if c not in exclude_cols and not any(c.startswith(pref) for pref in blacklist)
    ]

    include_flags = []
    if args.include_viral_load:
        include_flags.append('viral_load')

    if args.include_is_mrc:
        include_flags.append('is_mrc')
        # Put is_mrc as first feature in base_features
        if 'is_mrc' in base_features:
            base_features.remove('is_mrc')
        base_features = ['is_mrc'] + base_features

    results = evaluate_combos(df, base_features, include_flags=include_flags)

    label_parts = ['base']
    if args.include_is_mrc:
        label_parts.append('mrc')
    if args.include_viral_load:
        label_parts.append('vl')
    label = "_".join(label_parts)

    csv_path = save_csv(results, args.output, label)
    txt_path = save_txt_best(results, args.output, f'best_features_{label}.txt')

    print("Saved outputs:")
    print(f" - {csv_path}")
    print(f" - {txt_path}")

if __name__ == '__main__':
    main()
