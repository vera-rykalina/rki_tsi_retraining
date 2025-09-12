#!/usr/bin/env python3

import argparse
import os
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.impute import KNNImputer
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import logging

logging.basicConfig(level=logging.INFO)
loginfo = logging.info

def _impute_knn(X, k=3):
    '''Use k nearest rows which have a feature to fill in each row's missing features.'''
    if np.isnan(X).sum() == 0:
        loginfo('No missing data to impute, continuing.')
        return X
    nm = (np.isnan(X).sum(axis=1) != 0).sum()
    loginfo('Using KNN (3 nearest neighbours) to impute {} missing values ({} sample{})'.format(
        np.isnan(X).sum(), nm, ('', 's')[nm > 1]))
    return KNNImputer(n_neighbors=k).fit_transform(X)

def prepare_data(df, features_list, is_amplicons=False):
    '''Prepare feature matrix for modeling.'''
    seqtype = pd.DataFrame(index=df.index)
    seqtype['is_mrc'] = int(is_amplicons)
    data_for_model = pd.concat((seqtype, df), axis=1)[features_list]
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data_for_model)
    data_imputed = _impute_knn(data_scaled)
    return data_imputed, scaler

def train_model(input_file, features_file, output_folder, report_folder,
                amplicons=False, model_name="IGS"):
    loginfo(f"Loading data from {input_file}")
    df = pd.read_csv(input_file, index_col=0)

    loginfo(f"Loading feature list from {features_file}")
    with open(features_file) as f:
        features_list = [line.strip() for line in f if line.strip()]

    if 'known_tsi_years' not in df.columns:
        raise ValueError("Target column 'known_tsi_years' not found in input data.")

    y_raw = df['known_tsi_years'].values
    y = np.sqrt(y_raw)

    data_for_model_imp, scaler = prepare_data(df, features_list, is_amplicons=amplicons)

    outer_cv = KFold(n_splits=5, shuffle=True, random_state=42)
    inner_cv = KFold(n_splits=5, shuffle=True, random_state=42)

    param_grid = {
        'n_estimators': [100, 200, 1000],
        'max_depth': [3, 5, 7],
        'min_samples_split': [2, 5, 10],
        'min_samples_leaf': [1, 2, 4]
    }

    maes = []
    rmses = []
    r2s = []

    preds_all = []

    best_params_list = []
    best_maes_list = []

    fold_num = 1
    for train_idx, test_idx in outer_cv.split(data_for_model_imp):
        loginfo(f"Outer CV fold {fold_num} starting...")

        X_train, X_test = data_for_model_imp[train_idx], data_for_model_imp[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        rf = RandomForestRegressor(random_state=42, n_jobs=-1)
        grid = GridSearchCV(rf, param_grid, cv=inner_cv, scoring='neg_mean_absolute_error', n_jobs=-1)
        grid.fit(X_train, y_train)

        best_params = grid.best_params_
        best_score = -grid.best_score_  # neg MAE to positive MAE
        best_params_list.append(best_params)
        best_maes_list.append(best_score)

        best_rf = grid.best_estimator_

        # Per-tree predictions for test samples
        all_tree_preds = np.array([tree.predict(X_test) for tree in best_rf.estimators_])  # shape: (n_trees, n_samples_test)

        # Mean prediction per sample (sqrt scale)
        y_test_pred = np.mean(all_tree_preds, axis=0)

        # 95% prediction interval in sqrt scale
        lower_bounds = np.percentile(all_tree_preds, 2.5, axis=0)
        upper_bounds = np.percentile(all_tree_preds, 97.5, axis=0)

        # Convert back to original scale
        y_test_pred_linear = y_test_pred ** 2
        y_test_linear = y_test ** 2
        lower_linear = lower_bounds ** 2
        upper_linear = upper_bounds ** 2

        mae = mean_absolute_error(y_test_linear, y_test_pred_linear)
        rmse = np.sqrt(mean_squared_error(y_test_linear, y_test_pred_linear))
        r2 = r2_score(y_test_linear, y_test_pred_linear)

        maes.append(mae)
        rmses.append(rmse)
        r2s.append(r2)

        # Save results with intervals for reporting
        for i, idx in enumerate(test_idx):
            preds_all.append({
                'sample_index': idx,
                'true': y_test_linear[i],
                'pred': y_test_pred_linear[i],
                'lower_95': lower_linear[i],
                'upper_95': upper_linear[i],
                'abs_error': abs(y_test_linear[i] - y_test_pred_linear[i])
            })

        loginfo(f"Fold {fold_num} results: MAE={mae:.4f}, RMSE={rmse:.4f}, R2={r2:.4f}")
        fold_num += 1

    # Aggregate CV results
    mae_mean, mae_std = np.mean(maes), np.std(maes)
    rmse_mean, rmse_std = np.mean(rmses), np.std(rmses)
    r2_mean, r2_std = np.mean(r2s), np.std(r2s)

    loginfo(f"Nested CV results: MAE={mae_mean:.4f}±{mae_std:.4f}, RMSE={rmse_mean:.4f}±{rmse_std:.4f}, R2={r2_mean:.4f}±{r2_std:.4f}")

    # Pick best params based on aggregated CV MAE
    best_index = np.argmin(best_maes_list)
    best_params_agg = best_params_list[best_index]
    loginfo(f"Best hyperparameters from aggregated CV: {best_params_agg}")

    # Retrain final model on full dataset using aggregated best hyperparameters
    final_rf = RandomForestRegressor(random_state=42, n_jobs=-1, **best_params_agg)
    final_rf.fit(data_for_model_imp, y)

    y_pred_full = final_rf.predict(data_for_model_imp)
    y_pred_full_linear = y_pred_full ** 2
    y_linear = y ** 2
    abs_errors_full = np.abs(y_linear - y_pred_full_linear)

    # Train residuals error model (absolute residuals in sqrt scale)
    train_residuals_full = np.abs(y - y_pred_full)
    err_model = RandomForestRegressor(random_state=42, **best_params_agg)
    err_model.fit(data_for_model_imp, train_residuals_full)

    # Save model components
    os.makedirs(output_folder, exist_ok=True)
    with open(os.path.join(output_folder, "full_scaler.selfscale.pkl"), "wb") as f:
        pickle.dump(scaler, f)
    with open(os.path.join(output_folder, f"full_model.selfscale.{model_name}.pkl"), "wb") as f:
        pickle.dump(final_rf, f)
    with open(os.path.join(output_folder, f"err_model_mae.{model_name}.pkl"), "wb") as f:
        pickle.dump(err_model, f)
    with open(os.path.join(output_folder, f"feature_list_{model_name}.txt"), "w") as f:
        for feat in features_list:
            f.write(f"{feat}\n")

    # Prepare report folder and save metrics + plot
    os.makedirs(report_folder, exist_ok=True)
    metrics_path = os.path.join(report_folder, "metrics.csv")
    plot_path = os.path.join(report_folder, "performance_plot.png")

    # Convert preds_all list of dicts to DataFrame sorted by sample index
    preds_df = pd.DataFrame(preds_all).sort_values('sample_index').reset_index(drop=True)

    # Summary metrics dataframe as header lines
    summary_metrics = pd.DataFrame({
        'sample_index': [np.nan]*3,
        'true': ['MAE', 'RMSE', 'R2'],
        'pred': [f"{mae_mean:.4f} ± {mae_std:.4f}",
                 f"{rmse_mean:.4f} ± {rmse_std:.4f}",
                 f"{r2_mean:.4f} ± {r2_std:.4f}"],
        'lower_95': [np.nan]*3,
        'upper_95': [np.nan]*3,
        'abs_error': [np.nan]*3
    })

    full_report_df = pd.concat([summary_metrics, preds_df], ignore_index=True)
    full_report_df.to_csv(metrics_path, index=False)

    # Plot true vs predicted with identity line
    plt.figure(figsize=(8, 6))
    plt.scatter(y_linear, y_pred_full_linear, alpha=0.7, edgecolors='k')
    plt.plot([y_linear.min(), y_linear.max()], [y_linear.min(), y_linear.max()], 'r--')
    plt.xlabel('True known_tsi_years')
    plt.ylabel('Predicted known_tsi_years')
    plt.title('Performance: True vs Predicted (full data)')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(plot_path)
    plt.close()

    loginfo(f"Metrics saved to {metrics_path}")
    loginfo(f"Performance plot saved to {plot_path}")
    loginfo(f"Model components saved in {output_folder}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Train Random Forest model for known_tsi_years prediction."
    )
    parser.add_argument("--input", required=True,
                        help="Input CSV file with features and target (known_tsi_years).")
    parser.add_argument("--features", required=True,
                        help="TXT file with list of features to use for training.")
    parser.add_argument("--modeldir", required=True,
                        help="Folder name to save trained model components.")
    parser.add_argument("--report", required=True,
                        help="Folder name to save metrics CSV and performance plot.")
    parser.add_argument("--amplicons", action="store_true",
                        help="Set is_mrc=1 (amplicons data). Default is False.")
    parser.add_argument("--modelname", default="IGS",
                        help="Model name used in saved file names (default: IGS).")

    args = parser.parse_args()

    train_model(
        input_file=args.input,
        features_file=args.features,
        output_folder=args.modeldir,
        report_folder=args.report,
        amplicons=args.amplicons,
        model_name=args.modelname
    )
