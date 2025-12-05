#!/usr/bin/env python3

"""
Convert sRNATarBase-like TSV/CSV to srna_validated_interactions.csv schema.

Usage:
    python convert_srna_db_to_csv.py input.tsv output.csv
"""

import sys
import pandas as pd

def main(input_path, output_path):
    # Auto-detect delimiter (TSV/CSV)
    if input_path.endswith(".tsv"):
        df = pd.read_csv(input_path, sep="\t")
    else:
        df = pd.read_csv(input_path)

    # Try to guess column names from common patterns
    # Adjust these mappings based on your actual file
    colmap_candidates = {
        "organism": ["Organism", "Species", "organism"],
        "srna_name": ["sRNA", "srna", "sRNA_name"],
        "mrna_name": ["Target", "Target_gene", "mRNA", "target_name"],
        "evidence_type": ["Evidence", "Evidence_type", "Validation"],
        "pmid": ["PMID", "PubMed", "pmid"],
        "context": ["Condition", "Experiment_condition", "Context"],
    }

    def find_col(possible):
        for name in possible:
            if name in df.columns:
                return name
        return None

    org_col = find_col(colmap_candidates["organism"])
    srna_col = find_col(colmap_candidates["srna_name"])
    mrna_col = find_col(colmap_candidates["mrna_name"])
    ev_col = find_col(colmap_candidates["evidence_type"])
    pmid_col = find_col(colmap_candidates["pmid"])
    ctx_col = find_col(colmap_candidates["context"])

    if not (org_col and srna_col and mrna_col):
        raise ValueError("Could not infer essential columns (organism/sRNA/target). Adjust mappings in script.")

    out = {}

    out["organism"] = df[org_col]
    out["srna_name"] = df[srna_col]
    out["mrna_name"] = df[mrna_col]

    out["evidence_type"] = df[ev_col] if ev_col else ""
    out["reference_doi_or_pmid"] = df[pmid_col].apply(lambda x: f"PMID:{x}" if pd.notna(x) else "") if pmid_col else ""
    out["interaction_context"] = df[ctx_col] if ctx_col else ""

    # Initialize model-related columns as empty; you'll fill them later with your pipeline
    out["deltaG_predicted"] = ""
    out["gc_content_srna"] = ""
    out["gc_content_mrna"] = ""
    out["seed_length"] = ""
    out["has_seed"] = ""
    out["hybrid_length"] = ""

    # Simple label heuristic: treat all imported as strong positives for now (you can refine later)
    out["label_strength"] = 1.0

    # Optional notes column
    out["notes"] = df.get("Notes", "")

    out_df = pd.DataFrame(out)
    out_df.to_csv(output_path, index=False)
    print(f"Saved converted dataset to: {output_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_srna_db_to_csv.py input.tsv output.csv")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
