#!/usr/bin/env python3

"""
Parse an HTML table (e.g. from sRNATarBase) and convert to srna_validated_interactions.csv.
"""

from bs4 import BeautifulSoup
import pandas as pd
import sys

def html_table_to_csv(html_path, output_csv):
    with open(html_path, "r", encoding="utf-8") as f:
        soup = BeautifulSoup(f, "html.parser")

    table = soup.find("table")
    if table is None:
        raise ValueError("No <table> element found in HTML.")

    # Extract header
    headers = [th.get_text(strip=True) for th in table.find_all("th")]
    rows = []
    for tr in table.find_all("tr")[1:]:
        cells = [td.get_text(strip=True) for td in tr.find_all(["td", "th"])]
        if len(cells) != len(headers):
            continue
        rows.append(cells)

    raw_df = pd.DataFrame(rows, columns=headers)

    # Map as before (reuse logic from previous script; simplified here)
    colmap = {
        "organism": next((c for c in raw_df.columns if "Organism" in c or "Species" in c), None),
        "srna_name": next((c for c in raw_df.columns if "sRNA" in c), None),
        "mrna_name": next((c for c in raw_df.columns if "Target" in c or "mRNA" in c), None),
        "evidence_type": next((c for c in raw_df.columns if "Evidence" in c), None),
        "pmid": next((c for c in raw_df.columns if "PMID" in c or "PubMed" in c), None),
    }

    if not (colmap["organism"] and colmap["srna_name"] and colmap["mrna_name"]):
        raise ValueError("Could not identify core columns in HTML table")

    out = {}
    out["organism"] = raw_df[colmap["organism"]]
    out["srna_name"] = raw_df[colmap["srna_name"]]
    out["mrna_name"] = raw_df[colmap["mrna_name"]]
    out["evidence_type"] = raw_df[colmap["evidence_type"]] if colmap["evidence_type"] else ""
    out["reference_doi_or_pmid"] = raw_df[colmap["pmid"]].apply(lambda x: f"PMID:{x}") if colmap["pmid"] else ""

    out["interaction_context"] = ""
    out["deltaG_predicted"] = ""
    out["gc_content_srna"] = ""
    out["gc_content_mrna"] = ""
    out["seed_length"] = ""
    out["has_seed"] = ""
    out["hybrid_length"] = ""
    out["label_strength"] = 1.0
    out["notes"] = ""

    df = pd.DataFrame(out)
    df.to_csv(output_csv, index=False)
    print(f"Saved parsed CSV to {output_csv}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python html_table_to_csv.py input.html output.csv")
        sys.exit(1)

    html_table_to_csv(sys.argv[1], sys.argv[2])
