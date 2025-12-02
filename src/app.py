from fastapi import FastAPI
import subprocess, shutil, os, uuid, json

app = FastAPI(title="sRNA Predictor Demo")

@app.get("/")
def root():
    return {"status": "ready", "note": "sRNA predictor demo API"}

@app.post("/predict_simple")
def predict_simple(payload: dict):
    """
    POST JSON: {"srna":"AUGCUAUGCUA", "mrna":"GGGGAUACGAUAGCUAGCUA"}
    Runs IntaRNA on temp files and returns raw output.
    """

    # Extract input
    srna = payload.get("srna", "").strip()
    mrna = payload.get("mrna", "").strip()

    if not srna or not mrna:
        return {"error": "Provide both 'srna' and 'mrna' sequences as strings."}

    # Ensure IntaRNA exists
    if not shutil.which("IntaRNA"):
        return {"error": "IntaRNA not found on PATH. Check installation."}

    # Create temporary FASTA files
    qf = f"/tmp/q_{uuid.uuid4().hex}.fa"
    tf = f"/tmp/t_{uuid.uuid4().hex}.fa"

    try:
        with open(qf, "w") as f:
            f.write(">srna\n" + srna + "\n")

        with open(tf, "w") as f:
            f.write(">mrna\n" + mrna + "\n")

        # Run IntaRNA
        proc = subprocess.run(
            ["IntaRNA", "-q", qf, "-t", tf],
            capture_output=True,
            text=True,
            timeout=30
        )

        out = proc.stdout
        err = proc.stderr

    finally:
        # Clean up
        if os.path.exists(qf):
            os.remove(qf)
        if os.path.exists(tf):
            os.remove(tf)

    return {"stdout": out, "stderr": err}
