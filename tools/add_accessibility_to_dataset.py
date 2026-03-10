import pandas as pd
import subprocess
import tempfile
import os

DATA_PATH = "data/srna_augmented.csv"


def run_rnaplfold(sequence):
    """
    Run RNAplfold and return accessibility score
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta = os.path.join(tmpdir, "seq.fa")

        with open(fasta, "w") as f:
            f.write(">seq\n")
            f.write(sequence + "\n")

        subprocess.run(
            ["RNAplfold", "-W", "80", "-L", "40", "-u", "30"],
            cwd=tmpdir,
            input=sequence,
            text=True,
            capture_output=True
        )

        lunp = os.path.join(tmpdir, "plfold_lunp")

        if not os.path.exists(lunp):
            return None

        probs = []

        with open(lunp) as f:
            for line in f:
                if line.startswith("#"):
                    continue

                cols = line.strip().split()

                # probability that region length=1 is unpaired
                p_unpaired = float(cols[1])
                probs.append(p_unpaired)

        if len(probs) == 0:
            return None

        return sum(probs) / len(probs)


df = pd.read_csv(DATA_PATH)

accessibility_scores = []

for seq in df["mrna_sequence"]:
    if pd.isna(seq) or seq == "Same as above":
        accessibility_scores.append(None)
        continue

    try:
        score = run_rnaplfold(seq)
    except:
        score = None

    accessibility_scores.append(score)


df["target_accessibility"] = accessibility_scores

df.to_csv(DATA_PATH, index=False)

print(" Added RNAplfold accessibility column")
