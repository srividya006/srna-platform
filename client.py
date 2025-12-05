"""
Simple Python client for the sRNA–mRNA Interaction Platform.
Requires: requests
"""

import requests
from typing import Any, Dict, List, Optional

BASE_URL = "http://localhost:8080"


def _post(path: str, json: Dict[str, Any]) -> Dict[str, Any]:
    url = f"{BASE_URL}{path}"
    resp = requests.post(url, json=json)
    resp.raise_for_status()
    return resp.json()


def rnafold(sequence: str, name: str = "seq1") -> Dict[str, Any]:
    """
    Call /rnafold endpoint.
    """
    payload = {
        "sequence": sequence,
        "name": name
    }
    return _post("/rnafold", payload)


def predict_intarna(
    srna: str,
    mrna: str,
    srna_name: str = "srna",
    mrna_name: str = "mrna",
    max_hits: int = 5,
    demo: bool = False,
) -> Dict[str, Any]:
    """
    Call /predict_intarna endpoint.
    """
    payload = {
        "srna": srna,
        "mrna": mrna,
        "srna_name": srna_name,
        "mrna_name": mrna_name,
        "max_hits": max_hits,
        "demo": demo,
    }
    return _post("/predict_intarna", payload)


def batch_predict(requests_list: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Call /batch_predict endpoint.

    requests_list: list of IntaRNARequest-like dicts, e.g.
        {
          "srna": "...",
          "mrna": "...",
          "srna_name": "...",
          "mrna_name": "...",
          "max_hits": 5,
          "demo": false
        }
    """
    payload = {
        "requests": requests_list
    }
    return _post("/batch_predict", payload)


def explain(
    srna: str,
    mrna: str,
    srna_name: str = "srna",
    mrna_name: str = "mrna",
    max_hits: int = 5,
    demo: bool = False,
) -> Dict[str, Any]:
    """
    Call /explain endpoint.
    Returns explanation for the top-ranked interaction.
    """
    payload = {
        "srna": srna,
        "mrna": mrna,
        "srna_name": srna_name,
        "mrna_name": mrna_name,
        "max_hits": max_hits,
        "demo": demo,
    }
    return _post("/explain", payload)


if __name__ == "__main__":
    # Example sequences
    srna_seq = "AUGCUACGUGAAGGCU"
    mrna_seq = "CUCCGCUUUCACGCGGAUUACG"

    print("=== /rnafold ===")
    rf = rnafold(srna_seq, name="example_srna")
    print("MFE:", rf.get("mfe"), "kcal/mol")
    print("GC-content:", rf.get("gc_content"))

    print("\n=== /predict_intarna ===")
    pi = predict_intarna(srna_seq, mrna_seq, "my_srna", "my_target", max_hits=3)
    print("Interactions:", len(pi.get("interactions", [])))
    if pi.get("interactions"):
        print("Top ΔG:", pi["interactions"][0].get("deltaG"))

    print("\n=== /batch_predict ===")
    bp = batch_predict([
        {
            "srna": srna_seq,
            "mrna": mrna_seq,
            "srna_name": "job1_srna",
            "mrna_name": "job1_target",
            "max_hits": 3,
            "demo": False
        },
        {
            "srna": "GGGACCUUAGCUAUGCC",
            "mrna": "AGGCUCCAUGGCAAGCGAC",
            "srna_name": "job2_srna",
            "mrna_name": "job2_target",
            "max_hits": 3,
            "demo": True
        }
    ])
    print("Batch count:", bp.get("count"))

    print("\n=== /explain ===")
    ex = explain(srna_seq, mrna_seq, "my_srna", "my_target", max_hits=5)
    print("Why rank 1?:", ex.get("explain", {}).get("why_rank1"))
    print("Relative affinity:", ex.get("explain", {}).get("relative_affinity"))
    print("Seed features:", ex.get("seed_features"))
