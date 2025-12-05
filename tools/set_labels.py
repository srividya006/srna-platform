import pandas as pd

PATH = "data/srna_target_validated.csv"

df = pd.read_csv(PATH, encoding="utf-8", encoding_errors="replace")
def map_validation(v):
    if not isinstance(v, str):
        return ""

    v = v.strip().lower()

    if "high" in v or "strong" in v:
        if "medium" in v:        # catches "medium-high"
            return 0.8
        return 1.0

    if "medium-high" in v:
        return 0.8

    if "medium" in v:
        return 0.7

    if "medium-low" in v:
        return 0.5

    if "low" in v:
        return 0.4

    if "very low" in v or "weak" in v:
        return 0.1

    return ""  # unknown label

df["label_strength"] = df["validation_level"].apply(map_validation)

df.to_csv(PATH, index=False)

print("✓ label_strength column updated.")
