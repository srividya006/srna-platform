import os
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error
import joblib

DATA_PATH = "../data/final_training_dataset_with_accessibility.csv"
MODEL_DIR = "/tmp/models"
MODEL_PATH = os.path.join(MODEL_DIR, "interaction_rf.joblib")

# Load data
df = pd.read_csv(DATA_PATH)

feature_cols = [
    "deltaG_predicted",
    "gc_content_srna",
    "gc_content_mrna",
    "seed_length",
    "has_seed",
    "hybrid_length",
    "target_accessibility",
]
target_col = "label_strength"

# Fill missing feature values
df["deltaG_predicted"] = df["deltaG_predicted"].fillna(0)
df["gc_content_srna"] = df["gc_content_srna"].fillna(0)
df["gc_content_mrna"] = df["gc_content_mrna"].fillna(0)
df["seed_length"] = df["seed_length"].fillna(0)
df["has_seed"] = df["has_seed"].fillna(0)
df["hybrid_length"] = df["hybrid_length"].fillna(0)
df["target_accessibility"] = df["target_accessibility"].fillna(0)

# Keep only rows with all features + label present
df_clean = df.dropna(subset=[target_col])

X = df_clean[feature_cols].astype(float)
y = df_clean[target_col].astype(float)

# Simple train / test split
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.25, random_state=42
)

model = RandomForestRegressor(
    n_estimators=200,
    max_depth=6,
    random_state=42
)

model.fit(X_train, y_train)

y_pred = model.predict(X_test)
print("R²:", r2_score(y_test, y_pred))
print("MSE:", mean_squared_error(y_test, y_pred))

os.makedirs(MODEL_DIR, exist_ok=True)
joblib.dump(model, MODEL_PATH)
print(f" Saved model to {MODEL_PATH}")
