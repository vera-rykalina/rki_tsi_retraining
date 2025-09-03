#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

def evaluate_model(input_csv, output_csv):
    # Load and validate input
    df = pd.read_csv(input_csv)
    required_cols = {"true_known_tsi_years", "predicted_known_tsi_years", "abs_error"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"CSV must contain columns: {required_cols}")

    # Extract values
    y_true = df["true_known_tsi_years"].values
    y_pred = df["predicted_known_tsi_years"].values
    errors = df["abs_error"].values

    # Calculate metrics
    mae = mean_absolute_error(y_true, y_pred)
    rmse = mean_squared_error(y_true, y_pred, squared=False)
    r2 = r2_score(y_true, y_pred)
    large_error_count = np.sum(errors > 2.0)

    # Create and save summary table
    results = pd.DataFrame({
        "Metric": ["MAE", "RMSE", "R2", "N_Large_Errors_>2y"],
        "Value": [mae, rmse, r2, large_error_count]
    })

    results.to_csv(output_csv, index=False)
    print(f"✅ Performance report saved to {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate model performance and save metrics to CSV.")
    parser.add_argument("--input", required=True, help="Input CSV file (metrics.csv)")
    parser.add_argument("--output", required=True, help="Output CSV file (performance.csv)")
    args = parser.parse_args()

    evaluate_model(args.input, args.output)