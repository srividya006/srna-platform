from fastapi import FastAPI
from pydantic import BaseModel
from typing import List, Optional
from .ml_model import (
    is_available as ml_is_available, 
    score_interaction as ml_score_interaction
)
import subprocess
import shutil
import os
import uuid
import csv

app = FastAPI(title="sRNA–mRNA Interaction Platform", version="0.5.0")


# ==========
# Utility functions
# ==========

def gc_content(seq: str) -> float:
    """
    Return GC fraction (0–1) for an RNA/DNA sequence.
    """
    if not seq:
        return 0.0
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq)


def find_seed_match(srna_seq: str, mrna_seq: str,
                    min_len: int = 6, max_len: int = 8) -> dict:
    """
    Simple seed detection heuristic.

    Looks for a perfect contiguous match of length 6–8 nt
    between sRNA and mRNA sequences (exact substring match).

    Returns a dict:
    {
        "has_seed": bool,
        "length": int,
        "srna_start": int or None,  # 1-based
        "mrna_start": int or None   # 1-based
    }
    """
    if not srna_seq or not mrna_seq:
        return {
            "has_seed": False,
            "length": 0,
            "srna_start": None,
            "mrna_start": None,
        }

    srna = srna_seq.upper()
    mrna = mrna_seq.upper()

    best = None

    # Prefer longer seeds first (8 -> 6), then earliest mRNA position
    for L in range(max_len, min_len - 1, -1):
        for i in range(len(srna) - L + 1):
            window = srna[i:i + L]
            j = mrna.find(window)
            if j != -1:
                candidate = {
                    "length": L,
                    "srna_start": i + 1,  # convert to 1-based
                    "mrna_start": j + 1,  # convert to 1-based
                }
                if best is None:
                    best = candidate
                else:
                    if candidate["length"] > best["length"]:
                        best = candidate
                    elif candidate["length"] == best["length"] and candidate["mrna_start"] < best["mrna_start"]:
                        best = candidate

    if best is None:
        return {
            "has_seed": False,
            "length": 0,
            "srna_start": None,
            "mrna_start": None,
        }

    best["has_seed"] = True
    return best


# ==========
# Data models
# ==========

class RNAFoldRequest(BaseModel):
    sequence: str
    name: Optional[str] = "seq1"


class RNAFoldResult(BaseModel):
    name: str
    sequence: str
    structure: str
    mfe: float
    raw_stdout: str


class IntaRNARequest(BaseModel):
    srna: str
    mrna: str
    srna_name: Optional[str] = "srna"
    mrna_name: Optional[str] = "mrna"
    max_hits: int = 5  # limit interactions returned
    demo: bool = False  # relaxed demo mode


class IntaRNAInteraction(BaseModel):
    rank: int
    srna_name: str
    mrna_name: str
    deltaG: float
    target_start: int
    target_end: int
    srna_start: int
    srna_end: int
    hybrid_length: Optional[int] = None
    raw_row: str


class IntaRNAResult(BaseModel):
    srna_name: str
    mrna_name: str
    interactions: List[IntaRNAInteraction]
    raw_stdout: str
    raw_stderr: str


class BatchIntaRNARequest(BaseModel):
    requests: List[IntaRNARequest]


# ==========
# Root / health
# ==========

@app.get("/")
def root():
    return {
        "status": "ok",
        "message": "sRNA–mRNA platform backend is running",
        "endpoints": [
            "/rnafold",
            "/predict_simple",
            "/predict_intarna",
            "/batch_predict",
            "/explain",
        ],
    }


# ==========
# RNAfold endpoint
# ==========

@app.post("/rnafold")
def rnafold_endpoint(req: RNAFoldRequest):
    """
    Run RNAfold on a single RNA sequence.
    Returns structure and minimum free energy (MFE), or an error message.
    """
    if not shutil.which("RNAfold"):
        return {"error": "RNAfold binary not found on PATH"}

    fasta_input = f">{req.name}\n{req.sequence.strip()}\n"

    try:
        proc = subprocess.run(
            ["RNAfold", "--noPS"],
            input=fasta_input,
            text=True,
            capture_output=True,
            timeout=30
        )
    except Exception as e:
        return {"error": f"Failed to run RNAfold: {e}"}

    stdout = proc.stdout.strip()
    stderr = proc.stderr.strip()

    if proc.returncode != 0:
        return {
            "error": "RNAfold exited with non-zero status",
            "stdout": stdout,
            "stderr": stderr,
        }

    # RNAfold output lines: >name, sequence, "structure (  -x.xx)"
    lines = stdout.splitlines()
    if len(lines) < 3:
        return {
            "error": "Unexpected RNAfold output format",
            "stdout": stdout,
            "stderr": stderr,
        }

    name_line = lines[0].lstrip(">").strip()
    seq_line = lines[1].strip()
    struct_line = lines[2].strip()

    # Example struct_line: "........... (  0.00)"
    try:
        if "(" in struct_line and ")" in struct_line:
            paren_start = struct_line.rfind("(")
            paren_end = struct_line.rfind(")")
            structure = struct_line[:paren_start].strip()
            mfe_str = struct_line[paren_start + 1:paren_end].strip()
            mfe_val = float(mfe_str)
        else:
            return {
                "error": "No parentheses with MFE found in RNAfold output line",
                "line": struct_line,
                "stdout": stdout,
                "stderr": stderr,
            }
    except Exception as e:
        return {
            "error": "Could not parse RNAfold output line",
            "line": struct_line,
            "exception": str(e),
            "stdout": stdout,
            "stderr": stderr,
        }

    return {
        "name": name_line,
        "sequence": seq_line,
        "structure": structure.strip(),
        "mfe": mfe_val,
        "gc_content": gc_content(seq_line),
        "raw_stdout": stdout,
    }


# ==========
# Old raw IntaRNA endpoint (debug)
# ==========

@app.post("/predict_simple")
def predict_simple(payload: dict):
    """
    Simple/debug endpoint:
    POST JSON: {"srna": "...", "mrna": "..."}
    Returns raw IntaRNA stdout/stderr.
    """
    srna = payload.get("srna", "").strip()
    mrna = payload.get("mrna", "").strip()
    if not srna or not mrna:
        return {"error": "Provide 'srna' and 'mrna' in JSON body."}

    if not shutil.which("IntaRNA"):
        return {"error": "IntaRNA not found on PATH"}

    qf = f"/tmp/q_{uuid.uuid4().hex}.fa"
    tf = f"/tmp/t_{uuid.uuid4().hex}.fa"

    try:
        with open(qf, "w") as fq:
            fq.write(">srna\n" + srna + "\n")
        with open(tf, "w") as ft:
            ft.write(">mrna\n" + mrna + "\n")

        proc = subprocess.run(
            ["IntaRNA", "-q", qf, "-t", tf],
            capture_output=True,
            text=True,
            timeout=60,
        )
        out = proc.stdout
        err = proc.stderr
    finally:
        for path in (qf, tf):
            try:
                os.remove(path)
            except FileNotFoundError:
                pass

    return {"stdout": out, "stderr": err}


# ==========
# Core IntaRNA runner (reused by multiple endpoints)
# ==========

def run_intarna(req: IntaRNARequest) -> dict:
    """
    Core logic for running IntaRNA and parsing structured output.
    Reused by /predict_intarna, /batch_predict, and /explain.
    """
    if not shutil.which("IntaRNA"):
        return {"error": "IntaRNA binary not found on PATH"}

    srna_seq = req.srna.strip()
    mrna_seq = req.mrna.strip()
    if not srna_seq or not mrna_seq:
        return {"error": "Both 'srna' and 'mrna' sequences must be non-empty."}

    # Compute GC-content for both sequences (for ML/features)
    srna_gc = gc_content(srna_seq)
    mrna_gc = gc_content(mrna_seq)

    # Compute simple seed heuristic (6–8 bp exact match)
    seed_info = find_seed_match(srna_seq, mrna_seq)

    qf = f"/tmp/q_{uuid.uuid4().hex}.fa"
    tf = f"/tmp/t_{uuid.uuid4().hex}.fa"

    stdout = ""
    stderr = ""

    try:
        with open(qf, "w") as fq:
            fq.write(f">{req.srna_name}\n{srna_seq}\n")
        with open(tf, "w") as ft:
            ft.write(f">{req.mrna_name}\n{mrna_seq}\n")

        # Use CSV mode; let IntaRNA choose default columns.
        cmd = [
            "IntaRNA",
            "-q", qf,
            "-t", tf,
            "--outMode=C",
        ]

        # Relax constraints in demo mode to force more hits
        if req.demo:
            cmd.extend([
                "--noSeed",      # no strict seed needed
                "--acc=N",       # disable accessibility penalty
                "--outMaxE=10",  # allow interactions up to +10 kcal/mol
                "--outNumber=5", # request up to 5 ranked interaction predictions
            ])

        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120,
        )

        stdout = proc.stdout.strip()
        stderr = proc.stderr.strip()

    except Exception as e:
        return {
            "error": f"Exception while running IntaRNA: {e}",
            "stdout": stdout,
            "stderr": stderr,
        }
    finally:
        for path in (qf, tf):
            try:
                os.remove(path)
            except FileNotFoundError:
                pass

    # If IntaRNA exited with error, show details
    if proc.returncode != 0:
        return {
            "error": "IntaRNA exited with non-zero status",
            "exit_code": proc.returncode,
            "stdout": stdout,
            "stderr": stderr,
            "cmd": " ".join(cmd),
        }

    # Split lines, ignore empty ones and comments
    lines = [ln for ln in stdout.splitlines() if ln.strip() and not ln.startswith("#")]

    # Base response scaffold
    base_response = {
        "srna_name": req.srna_name,
        "mrna_name": req.mrna_name,
        "interactions": [],
        "raw_stdout": stdout,
        "raw_stderr": stderr,
        "gc_content_srna": srna_gc,
        "gc_content_mrna": mrna_gc,
        "seed_features": seed_info,
    }

    if not lines:
        # Absolutely no CSV lines at all
        base_response["note"] = "No CSV output from IntaRNA (no interactions or header)."
        return base_response

    # If only a header line like:
    # id1;start1;end1;id2;start2;end2;subseqDP;hybridDP;E
    if len(lines) == 1:
        base_response["note"] = "No interactions reported by IntaRNA (CSV header only, no data rows)."
        return base_response

    # Otherwise, parse CSV with ';' delimiter
    reader = csv.DictReader(lines, delimiter=';')
    interactions = []

    def safe_int(value):
        try:
            return int(value)
        except Exception:
            return None

    def safe_float(value):
        try:
            return float(value)
        except Exception:
            return 0.0

    for idx, row in enumerate(reader, start=1):
        try:
            rank = idx

            deltaG = safe_float(
                row.get("E")
                or row.get("E_total")
                or row.get("energy")
                or "0"
            )

            # IntaRNA 3.x typically uses: id1;start1;end1;id2;start2;end2;subseqDP;hybridDP;E
            q_start = safe_int(row.get("start1"))
            q_end = safe_int(row.get("end1"))
            tgt_start = safe_int(row.get("start2"))
            tgt_end = safe_int(row.get("end2"))

            hybrid = row.get("hybridDP", "") or row.get("hybrid", "")
            mrna_name = row.get("id1") or row.get("seq1") or req.mrna_name
            srna_name = row.get("id2") or row.get("seq2") or req.srna_name

            hybrid_len = len(hybrid.replace("&", "").replace(" ", "")) if hybrid else None

            # --- ML score (uses global seed_info) ---
            ml_score = None
            if ml_is_available():
                seed_len = seed_info.get("length") if seed_info else None
                seed_flag = 1.0 if (seed_info.get("has_seed") if seed_info else False) else 0.0

                ml_score = ml_score_interaction(
                    deltaG=deltaG,
                    gc_srna=srna_gc,
                    gc_mrna=mrna_gc,
                    seed_length=seed_len,
                    has_seed=seed_flag,
                    hybrid_length=hybrid_len,
                )

            interaction = {
                "rank": rank,
                "srna_name": srna_name,
                "mrna_name": mrna_name,
                "deltaG": deltaG,
                "target_start": tgt_start,
                "target_end": tgt_end,
                "srna_start": q_start,
                "srna_end": q_end,
                "hybrid_length": hybrid_len,
                "ml_score": ml_score,
                "raw_row": ";".join([row.get(k, "") for k in reader.fieldnames]),
            }
            interactions.append(interaction)
        except Exception:
            # Skip malformed row but continue
            continue

        if len(interactions) >= req.max_hits:
            break

    # --- Re-ranking by ML score if model is available ---
    if ml_is_available() and interactions:
        scored = [it for it in interactions if it.get("ml_score") is not None]
        unscored = [it for it in interactions if it.get("ml_score") is None]

        # High score = better
        scored.sort(key=lambda x: x["ml_score"], reverse=True)

        new_interactions = []
        for new_rank, it in enumerate(scored + unscored, start=1):
            it["rank"] = new_rank
            new_interactions.append(it)

        interactions = new_interactions

    base_response["interactions"] = interactions
    return base_response

# ==========
# Structured IntaRNA endpoint
# ==========

@app.post("/predict_intarna")
def predict_intarna(req: IntaRNARequest):
    """
    Structured IntaRNA wrapper.
    Uses IntaRNA 3.x CSV output (outMode=C) with ';' as delimiter.
    Returns parsed top interactions plus GC-content and seed features.
    """
    return run_intarna(req)


# ==========
# Batch IntaRNA endpoint
# ==========

@app.post("/batch_predict")
def batch_predict(req: BatchIntaRNARequest):
    """
    Run IntaRNA predictions in batch.

    Accepts:
    {
      "requests": [
        { "srna": "...", "mrna": "...", "srna_name": "...", "mrna_name": "...", "max_hits": 5, "demo": false },
        ...
      ]
    }

    Returns:
    {
      "count": N,
      "results": [ <run_intarna result>, ... ]
    }
    """
    results = []
    for item in req.requests:
        res = run_intarna(item)
        results.append(res)

    return {
        "count": len(results),
        "results": results,
    }


# ==========
# Explainability endpoint
# ==========

@app.post("/explain")
def explain_endpoint(req: IntaRNARequest):
    """
    Explain the top-ranked IntaRNA interaction using ΔG, seed features, and rank.

    Returns:
    - top_interaction
    - explain: {
        "why_rank1": "...",
        "seed_match": bool,
        "relative_affinity": 0–1,
        "deltaG_rank1": float,
        "rank1_index": int
      }
    - GC and seed features
    - raw_stdout, raw_stderr
    """
    base = run_intarna(req)

    # If run_intarna returned an error, bubble it up directly
    if "error" in base:
        return base

    interactions = base.get("interactions", [])

    # If no interactions, still return feature context
    if not interactions:
        return {
            "srna_name": base.get("srna_name"),
            "mrna_name": base.get("mrna_name"),
            "top_interaction": None,
            "explain": {
                "why_rank1": "No interactions available to explain.",
                "seed_match": False,
                "relative_affinity": 0.0,
                "deltaG_rank1": None,
                "rank1_index": None,
            },
            "gc_content_srna": base.get("gc_content_srna"),
            "gc_content_mrna": base.get("gc_content_mrna"),
            "seed_features": base.get("seed_features"),
            "raw_stdout": base.get("raw_stdout"),
            "raw_stderr": base.get("raw_stderr"),
            "note": base.get("note", "No interactions returned by IntaRNA."),
        }

    # Compute relative affinity for rank1 based on ΔG distribution
    energies = [it.get("deltaG", 0.0) for it in interactions]
    best = min(energies)  # most negative / strongest
    worst = max(energies)  # weakest

    top = interactions[0]
    top_energy = top.get("deltaG", 0.0)

    if worst == best:
        rel_aff = 1.0
    else:
        # Map such that best → 1, worst → 0
        rel_aff = (worst - top_energy) / (worst - best)
        if rel_aff < 0:
            rel_aff = 0.0
        if rel_aff > 1:
            rel_aff = 1.0

    seed_features = base.get("seed_features", {}) or {}
    seed_match = bool(seed_features.get("has_seed"))

    # Build human-readable explanation string
    reasons = []
    if top_energy == best:
        reasons.append("lowest ΔG among predicted interactions")
    if seed_match:
        reasons.append("presence of a contiguous seed-like match")
    hybrid_length = top.get("hybrid_length")
    if hybrid_length is not None and hybrid_length <= 25:
        reasons.append("compact hybridization window")

    if not reasons:
        reasons.append("top IntaRNA rank based on energy and interaction score")

    why_rank1 = "; ".join(reasons)

    explain = {
        "why_rank1": why_rank1,
        "seed_match": seed_match,
        "relative_affinity": round(rel_aff, 3),
        "deltaG_rank1": top_energy,
        "rank1_index": top.get("rank", 1),
    }

    return {
        "srna_name": base.get("srna_name"),
        "mrna_name": base.get("mrna_name"),
        "top_interaction": top,
        "explain": explain,
        "gc_content_srna": base.get("gc_content_srna"),
        "gc_content_mrna": base.get("gc_content_mrna"),
        "seed_features": seed_features,
        "raw_stdout": base.get("raw_stdout"),
        "raw_stderr": base.get("raw_stderr"),
    }

