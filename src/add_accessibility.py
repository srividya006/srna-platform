import pandas as pd

# Load dataset
df = pd.read_csv("../data/final_training_dataset.csv", encoding="latin1")

# Function to compute accessibility (AU-based approximation)
def compute_accessibility(seq):
    if pd.isna(seq) or len(seq) == 0:
        return 0
    seq = seq.upper().replace("T", "U")  # normalize DNA → RNA
    au_count = seq.count("A") + seq.count("U")
    return au_count / len(seq)

# Apply to mRNA sequences
df["target_accessibility"] = df["mrna_sequence"].apply(compute_accessibility)

# Save updated dataset
df.to_csv("../data/final_training_dataset_with_accessibility.csv", index=False)

print(" target_accessibility added successfully!")
