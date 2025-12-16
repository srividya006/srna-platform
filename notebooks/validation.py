import json
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
from joblib import load
from sklearn.metrics import (
    r2_score,
    mean_squared_error,
    mean_absolute_error,
    roc_auc_score,
    average_precision_score,
)

# -------- Paths --------
ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "data" / "srna_augmented.csv"
RF_MODEL_PATH = ROOT / "models" / "interaction_rf.joblib"
XGB_MODEL_PATH = ROOT / "models" / "interaction_xgb.joblib"  # optional
METRICS_PATH = ROOT / "data" / "validation_metrics.json"


def load_data():
    print(f"📂 Loading dataset from: {DATA_PATH}")
    df = pd.read_csv(DATA_PATH)

    feature_cols = [
        "deltaG_predicted",
        "gc_content_srna",
        "gc_content_mrna",
        "seed_length",
        "has_seed",
        "hybrid_length",
    ]
    target_col = "label_strength"

    df_clean = df.dropna(subset=feature_cols + [target_col]).copy()
    X = df_clean[feature_cols]
    y = df_clean[target_col].astype(float)

    print(f"✅ Clean dataset: {len(df_clean)} rows usable for validation")
    return df_clean, X, y


def compute_regression_metrics(y_true, y_pred):
    mse = mean_squared_error(y_true, y_pred)
    rmse = np.sqrt(mse)
    mae = mean_absolute_error(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)

    return {
        "r2": float(r2),
        "rmse": float(rmse),
        "mae": float(mae),
    }


def compute_binary_metrics(y_true_cont, scores, threshold=0.8):
    """
    Treat >= threshold as "strong" interaction (1), else 0
    """
    y_bin = (y_true_cont >= threshold).astype(int)
    try:
        roc_auc = roc_auc_score(y_bin, scores)
        ap = average_precision_score(y_bin, scores)
        return {
            "roc_auc": float(roc_auc),
            "average_precision": float(ap),
        }
    except ValueError:
        # Not enough positive/negative cases
        return {
            "roc_auc": None,
            "average_precision": None,
        }


def main():
    df_clean, X, y = load_data()

    results = {
        "dataset": {
            "n_samples": int(len(df_clean)),
            "feature_columns": list(X.columns),
            "target": "label_strength",
        },
        "models": {},
        "generated_at": datetime.utcnow().isoformat() + "Z",
    }

    # ---------- 1) RandomForest ----------
    rf_available = RF_MODEL_PATH.exists()
    if rf_available:
        print(f" Loading RandomForest model from {RF_MODEL_PATH}")
        rf = load(RF_MODEL_PATH)
        y_pred_rf = rf.predict(X)
        reg_rf = compute_regression_metrics(y, y_pred_rf)
        bin_rf = compute_binary_metrics(y, y_pred_rf)

        results["models"]["rf"] = {
            "name": "RandomForest",
            "available": True,
            "regression": reg_rf,
            "binary": bin_rf,
        }
    else:
        print("⚠ RF model not found")
        results["models"]["rf"] = {
            "name": "RandomForest",
            "available": False,
        }

    # ---------- 2) XGBoost (optional) ----------
    xgb_available = XGB_MODEL_PATH.exists()
    if xgb_available:
        try:
            print(f" Loading XGBoost model from {XGB_MODEL_PATH}")
            xgb = load(XGB_MODEL_PATH)
            y_pred_xgb = xgb.predict(X)
            reg_xgb = compute_regression_metrics(y, y_pred_xgb)
            bin_xgb = compute_binary_metrics(y, y_pred_xgb)

            results["models"]["xgb"] = {
                "name": "XGBoost",
                "available": True,
                "regression": reg_xgb,
                "binary": bin_xgb,
            }
        except Exception as e:
            print(f"Failed to load/use XGBoost model: {e}")
            results["models"]["xgb"] = {
                "name": "XGBoost",
                "available": False,
                "error": str(e),
            }
    else:
        print("ℹ No XGBoost model file found; skipping.")
        results["models"]["xgb"] = {
            "name": "XGBoost",
            "available": False,
        }

    # ---------- 3) IntaRNA ΔG baseline ----------
    # Use -ΔG as score, min-max scaled to [0,1]
    if "deltaG_predicted" in df_clean.columns:
        dg = df_clean["deltaG_predicted"].astype(float).to_numpy()
        dg = np.nan_to_num(dg)
        raw_score = -dg  # more negative ΔG → stronger → higher score

        if raw_score.max() == raw_score.min():
            baseline_pred = np.full_like(raw_score, y.mean())
        else:
            # scale to [0,1] for comparison
            baseline_pred = (raw_score - raw_score.min()) / (
                raw_score.max() - raw_score.min()
            )

        reg_base = compute_regression_metrics(y, baseline_pred)
        bin_base = compute_binary_metrics(y, baseline_pred)

        results["models"]["intarna_baseline"] = {
            "name": "IntaRNA ΔG baseline",
            "available": True,
            "regression": reg_base,
            "binary": bin_base,
        }
    else:
        results["models"]["intarna_baseline"] = {
            "name": "IntaRNA ΔG baseline",
            "available": False,
            "error": "deltaG_predicted column missing",
        }

    # ---------- Save metrics JSON ----------
    METRICS_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(METRICS_PATH, "w") as f:
        json.dump(results, f, indent=2)

    print(f"✅ Saved validation metrics to: {METRICS_PATH}")

    # ---------- Save per-interaction predictions CSV ----------
    df_out = df_clean.copy()

    if rf_available:
        df_out["rf_score"] = y_pred_rf
    if xgb_available and "y_pred_xgb" in locals():
        df_out["xgb_score"] = y_pred_xgb

    if "deltaG_predicted" in df_clean.columns:
        try:
            df_out["intarna_baseline_score"] = baseline_pred
        except NameError:
            pass

    out_path = ROOT / "data" / "validation_with_predictions.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(out_path, index=False)
    print(f"✅ Saved per-interaction predictions to: {out_path}")

if __name__ == "__main__":
    main()
