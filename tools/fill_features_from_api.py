import pandas as pd
import requests

API_URL = "http://localhost:8080/predict_intarna"
INPUT_CSV = "data/srna_target_validated.csv"
OUTPUT_CSV = "data/srna_augmented.csv"

# IMPORTANT:
# Your CSV must have columns 'srna_sequence' and 'mrna_sequence'
# along with srna_name and mrna_name

df = pd.read_csv(INPUT_CSV)

updated_rows = []

for idx, row in df.iterrows():
    raw_srna = row.get("srna_sequence", "")
    raw_mrna = row.get("mrna_sequence", "")

    srna_seq = "" if (pd.isna(raw_srna) or str(raw_srna).strip().lower() in ("", "nan")) else str(raw_srna).strip()
    mrna_seq = "" if (pd.isna(raw_mrna) or str(raw_mrna).strip().lower() in ("", "nan")) else str(raw_mrna).strip()

    if not srna_seq or not mrna_seq:
        print(f"[WARN] Missing sequences for row {idx} ({row.get('srna_name')} → {row.get('mrna_name')}). Skipping API call.")
        updated_rows.append(row)
        continue

    payload = {
        "srna": srna_seq,
        "mrna": mrna_seq,
        "srna_name": row.get("srna_name", f"srna_{idx}"),
        "mrna_name": row.get("mrna_name", f"mrna_{idx}"),
        "max_hits": 1,
        "demo": False,
    }

    try:
        r = requests.post(API_URL, json=payload, timeout=60)
        res = r.json()
    except Exception as e:
        print(f"[ERROR] API error for row {idx}: {e}")
        updated_rows.append(row)
        continue

    if res.get("interactions"):
        best = res["interactions"][0]

        row["deltaG_predicted"] = best.get("deltaG", "")
        row["gc_content_srna"] = res.get("gc_content_srna", "")
        row["gc_content_mrna"] = res.get("gc_content_mrna", "")

        seed = res.get("seed_features", {}) or {}
        row["seed_length"] = seed.get("length", "")
        row["has_seed"] = int(seed.get("has_seed", False))
        row["hybrid_length"] = best.get("hybrid_length", "")
    else:
        print(f"[INFO] No interactions returned for row {idx} ({row.get('srna_name')} → {row.get('mrna_name')}).")

    updated_rows.append(row)

out_df = pd.DataFrame(updated_rows)
out_df.to_csv(OUTPUT_CSV, index=False)
print(f"✅ Saved augmented dataset with features to: {OUTPUT_CSV}")
