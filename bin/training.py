#!/usr/bin/env python3

import argparse
import os
import pickle
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.impute import KNNImputer
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.metrics import mean_absolute_error
import logging

logging.basicConfig(level=logging.INFO)
loginfo = logging.info


def _impute_knn(X, k=3):
    ''' Use k nearest rows which have a feature to fill in each row's missing features. '''
    if np.isnan(X).sum() == 0:
        loginfo('No missing data to impute, continuing.')
        return X
    nm = (np.isnan(X).sum(axis=1) != 0).sum()
    loginfo('Using KNN (3 nearest neighbours) to impute {} missing values ({} sample{})'.format(
        np.isnan(X).sum(), nm, ('', 's')[nm > 1]))
    X_knn = KNNImputer(n_neighbors=k).fit_transform(X)
    return X_knn


def prepare_data(df, features_list, is_amplicons=False):
    ''' Prepare feature matrix for modeling. '''
    seqtype = pd.DataFrame(index=df.index)
    seqtype['is_mrc'] = int(is_amplicons)
    data_for_model = pd.concat((seqtype, df), axis=1)[features_list]
    # Scale features
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data_for_model)
    # Impute missing values
    data_imputed = _impute_knn(data_scaled)
    return data_imputed, scaler


def train_model(input_file, features_file, output_folder, report_file=None, amplicons=False, model_name="IGS"):
    loginfo(f"Loading data from {input_file}")
    df = pd.read_csv(input_file, index_col=0)
    
    loginfo(f"Loading feature list from {features_file}")
    with open(features_file) as f:
        features_list = [line.strip() for line in f if line.strip()]
    
    # Check that target is present
    if 'known_tsi_years' not in df.columns:
        raise ValueError("Target column 'known_tsi_years' not found in input data.")
    y = df['known_tsi_years'].values
    
    # Prepare features (exclude target and viral_load even if present)
    data_for_model_imp, scaler = prepare_data(df, features_list, is_amplicons=amplicons)

    # Split data into train/test (80/20)
    X_train, X_test, y_train, y_test = train_test_split(
        data_for_model_imp, y, test_size=0.2, random_state=42)

    # Hyperparameter tuning
    param_grid = {
        'n_estimators': [100, 200],
        'max_depth': [5, 10, None],
        'min_samples_split': [2, 5]
    }
    rf = RandomForestRegressor(random_state=42)
    grid = GridSearchCV(rf, param_grid, cv=3, scoring='neg_mean_absolute_error', n_jobs=-1)
    loginfo("Starting hyperparameter tuning with GridSearchCV...")
    grid.fit(X_train, y_train)
    best_rf = grid.best_estimator_
    loginfo(f"Best params: {grid.best_params_}")
    loginfo(f"Best CV MAE: {-grid.best_score_:.4f}")

    # Train error model for MAE prediction
    y_train_pred = best_rf.predict(X_train)
    train_residuals = np.abs(y_train - y_train_pred)
    err_model = RandomForestRegressor(random_state=42, n_estimators=100)
    err_model.fit(X_train, train_residuals)

    # Save model components
    os.makedirs(output_folder, exist_ok=True)

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

    # Generate test set predictions report if requested
    if report_file:
        y_test_pred = best_rf.predict(X_test)
        test_report = pd.DataFrame({
            'true_known_tsi_years': y_test,
            'predicted_known_tsi_years': y_test_pred,
            'abs_error': np.abs(y_test - y_test_pred)
        })
        test_report.to_csv(report_file, index=False)
        loginfo(f"Test predictions report saved to {report_file}")

        # Also print final test MAE
        test_mae = mean_absolute_error(y_test, y_test_pred)
        loginfo(f"Test set MAE: {test_mae:.4f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train Random Forest model for known_tsi_years prediction.")
    parser.add_argument("--input", required=True, help="Input CSV file with features and target (known_tsi_years).")
    parser.add_argument("--features", required=True, help="TXT file with list of features to use for training.")
    parser.add_argument("--modeldir", required=True, help="Folder name to save trained model components.")
    parser.add_argument("--report", help="CSV filename to save test set predictions report (optional).")
    parser.add_argument("--amplicons", action="store_true", default=False, help="Set is_mrc=1 (amplicons data). Default is False.")
    parser.add_argument("--modelname", default="IGS", help="Model name used in saved file names (default: IGS).")
    args = parser.parse_args()

    train_model(args.input, args.features, args.modeldir, args.report, amplicons=args.amplicons, model_name=args.modelname)