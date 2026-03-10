import subprocess
import tempfile
import os


def compute_accessibility(sequence):

    with tempfile.TemporaryDirectory() as tmpdir:

        # adapt window size for short sequences
        W = min(80, len(sequence))

        cmd = [
            "RNAplfold",
            "-W", str(W),
            "-L", "40",
            "-u", "1"
        ]

        subprocess.run(
            cmd,
            input=(sequence + "\n").encode(),
            cwd=tmpdir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        lunp_file = os.path.join(tmpdir, "plfold_lunp")

        accessibility = []

        if not os.path.exists(lunp_file):
            return accessibility

        with open(lunp_file) as f:

            for line in f:

                if line.startswith("#"):
                    continue

                cols = line.strip().split()

                if len(cols) != 2:
                    continue

                try:
                    prob = float(cols[1])
                    accessibility.append(prob)
                except ValueError:
                    continue

        return accessibility


def region_accessibility(accessibility, start, end):

    region = accessibility[start:end]

    if not region:
        return 0

    return sum(region) / len(region)
