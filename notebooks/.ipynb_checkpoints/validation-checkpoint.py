import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from joblib import load

from sklearn.metrics import (
    r2_score,
    mean_squared_error,
    mean_absolute_error,
    roc_auc_score,
    average_precision_score,
)
from sklearn.metrics import roc_curve, precision_recall_curve

DATA_PATH = Path("data/srna_augmented.csv")
MODEL_PATH = Path("models/interaction_rf.joblib")

df = pd.read_csv(DATA_PATH)
print(df.shape)
df.head()

feature_cols = [
    "deltaG_predicted",
    "gc_content_srna",
    "gc_content_mrna",
    "seed_length",
    "has_seed",
    "hybrid_length",
]

target_col = "label_strength"

# Keep only rows with all needed values
df_clean = df.dropna(subset=feature_cols + [target_col]).copy()
df_clean.shape

model = load(MODEL_PATH)

X = df_clean[feature_cols]
y = df_clean[target_col].astype(float)

X.head(), y.head()

y_pred = model.predict(X)

print("First 5 predictions:", y_pred[:5])
r2 = r2_score(y, y_pred)
mse = mean_squared_error(y, y_pred)
rmse = np.sqrt(mse)
mae = mean_absolute_error(y, y_pred)

print(f"R²        : {r2:.3f}")
print(f"RMSE      : {rmse:.3f}")
print(f"MAE       : {mae:.3f}")

plt.figure(figsize=(6, 6))
plt.scatter(y, y_pred, alpha=0.7)
plt.xlabel("True label_strength")
plt.ylabel("Predicted label_strength")
plt.title("True vs Predicted label_strength")
plt.plot([0, 1], [0, 1], "r--")
plt.grid(True)
plt.show()

y_bin = (y >= 0.8).astype(int)

# Some models are regressors; we treat y_pred as "score"
scores = y_pred

roc_auc = roc_auc_score(y_bin, scores)
ap = average_precision_score(y_bin, scores)

print(f"ROC AUC            : {roc_auc:.3f}")
print(f"Average Precision  : {ap:.3f}")

fpr, tpr, thr = roc_curve(y_bin, scores)

plt.figure(figsize=(6, 6))
plt.plot(fpr, tpr, label=f"ROC AUC={roc_auc:.3f}")
plt.plot([0, 1], [0, 1], "k--")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve (strong vs non-strong)")
plt.legend()
plt.grid(True)
plt.show()

prec, rec, thr_pr = precision_recall_curve(y_bin, scores)

plt.figure(figsize=(6, 6))
plt.plot(rec, prec, label=f"AP={ap:.3f}")
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title("Precision–Recall Curve")
plt.legend()
plt.grid(True)
plt.show()


