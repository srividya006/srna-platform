import os
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error
import joblib

DATA_PATH = "data/srna_augmented.csv"
MODEL_DIR = "models"
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
]
target_col = "label_strength"

# Keep only rows with all features + label present
df_clean = df.dropna(subset=feature_cols + [target_col])

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
print(f"✅ Saved model to {MODEL_PATH}")
