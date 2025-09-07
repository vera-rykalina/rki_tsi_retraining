#!/usr/bin/env python3

import argparse
import os
import pickle
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.impute import KNNImputer
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV, train_test_split, KFold
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

def train_model(input_file, features_file, output_folder, report_file,
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
    X_train, X_test, y_train, y_test = train_test_split(
        data_for_model_imp, y, test_size=0.2, random_state=42
    )

    param_grid = {
        'n_estimators': [100, 200],
        'max_depth': [5, 10, None],
        'min_samples_split': [2, 5],
        'min_samples_leaf': [1, 2]
    }

    rf = RandomForestRegressor(random_state=42, n_jobs=8)
    grid = GridSearchCV(rf, param_grid, cv=10, scoring='neg_mean_absolute_error', n_jobs=8)
    loginfo("Starting hyperparameter tuning with GridSearchCV...")
    grid.fit(X_train, y_train)
    best_rf = grid.best_estimator_
    loginfo(f"Best params: {grid.best_params_}")
    loginfo(f"Best CV MAE (sqrt space): {-grid.best_score_:.4f}")

    # Train error model using best model's parameters
    y_train_pred = best_rf.predict(X_train)
    train_residuals = np.abs(y_train - y_train_pred)
    err_model = RandomForestRegressor(random_state=42, **grid.best_params_)
    err_model.fit(X_train, train_residuals)

    os.makedirs(output_folder, exist_ok=True)

    # Save components
    with open(os.path.join(output_folder, "full_scaler.selfscale.pkl"), "wb") as f:
        pickle.dump(scaler, f)
    with open(os.path.join(output_folder, f"full_model.selfscale.{model_name}.pkl"), "wb") as f:
        pickle.dump(best_rf, f)
    with open(os.path.join(output_folder, f"err_model_mae.{model_name}.pkl"), "wb") as f:
        pickle.dump(err_model, f)
    with open(os.path.join(output_folder, f"feature_list_{model_name}.txt"), "w") as f:
        for feat in features_list:
            f.write(f"{feat}\n")

    loginfo(f"Model training complete. Components saved in {output_folder}")

    # === Evaluate and write report ===
    y_test_pred = best_rf.predict(X_test)
    y_test_pred_linear = y_test_pred ** 2
    y_test_linear = y_test ** 2

    mae = mean_absolute_error(y_test_linear, y_test_pred_linear)
    rmse = np.sqrt(mean_squared_error(y_test_linear, y_test_pred_linear))
    r2 = r2_score(y_test_linear, y_test_pred_linear)

    test_report = pd.DataFrame({
        'true_known_tsi_years': y_test_linear,
        'predicted_known_tsi_years': y_test_pred_linear,
        'abs_error': np.abs(y_test_linear - y_test_pred_linear)
    })

    # Append performance metrics at the end
    test_report.loc[len(test_report.index)] = [np.nan, 'MAE', mae]
    test_report.loc[len(test_report.index)] = [np.nan, 'RMSE', rmse]
    test_report.loc[len(test_report.index)] = [np.nan, 'R2', r2]

    test_report.to_csv(report_file, index=False)
    loginfo(f"Test predictions + performance metrics saved to {report_file}")
    loginfo(f"MAE = {mae:.4f}, RMSE = {rmse:.4f}, R² = {r2:.4f}")

    # === Retrain final models on full data if thresholds met ===
    if mae < 1.0 and r2 >= 0.5:
        loginfo(f"Performance thresholds met (MAE={mae:.3f} < 1.0 and R²={r2:.3f} ≥ 0.5).")
        loginfo("Retraining final model and error model on full dataset with best hyperparameters...")

        X_all, _ = prepare_data(df, features_list, is_amplicons=amplicons)

        final_rf = RandomForestRegressor(random_state=42, n_jobs=8, **grid.best_params_)
        final_rf.fit(X_all, y)

        y_all_pred = final_rf.predict(X_all)
        train_residuals_all = np.abs(y - y_all_pred)
        final_err_model = RandomForestRegressor(random_state=42, **grid.best_params_)
        final_err_model.fit(X_all, train_residuals_all)

        # Overwrite saved models
        with open(os.path.join(output_folder, f"full_model.selfscale.{model_name}.pkl"), "wb") as f:
            pickle.dump(final_rf, f)
        with open(os.path.join(output_folder, f"err_model_mae.{model_name}.pkl"), "wb") as f:
            pickle.dump(final_err_model, f)

        loginfo("Final models retrained and saved on full dataset.")
    else:
        loginfo(f"Performance thresholds NOT met (MAE={mae:.3f}, R²={r2:.3f}). Skipping final retraining.")

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
                        help="CSV filename to save test set predictions report (mandatory).")
    parser.add_argument("--amplicons", action="store_true",
                        help="Set is_mrc=1 (amplicons data). Default is False.")
    parser.add_argument("--modelname", default="IGS",
                        help="Model name used in saved file names (default: IGS).")

    args = parser.parse_args()

    train_model(
        input_file=args.input,
        features_file=args.features,
        output_folder=args.modeldir,
        report_file=args.report,
        amplicons=args.amplicons,
        model_name=args.modelname
    )