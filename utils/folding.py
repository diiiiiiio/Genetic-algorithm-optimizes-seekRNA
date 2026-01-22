import os
import subprocess
import tempfile
import re
from statistics import mean

def fold_rnafold(C):
    """RNAfold -> (dot-bracket, dG)"""
    p = subprocess.run(
        ["RNAfold", "--noPS"],
        input=C + "\n",
        capture_output=True,
        text=True,
    )
    out = p.stdout.strip().splitlines()
    if len(out) < 2:
        return None, None
    line = out[1]
    struct = line.split()[0]
    try:
        e = float(line[line.rfind("(")+1:line.rfind(")")])
    except ValueError:
        e = 0.0
    return struct, e


def fold_mxfold2(C):
    """mxfold2 -> (dot-bracket, score)"""
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as f:
        f.write(">C\n")
        f.write(C + "\n")
        fa_path = f.name

    try:
        p = subprocess.run(
            ["mxfold2", "predict", fa_path],
            capture_output=True,
            text=True,
        )

        if p.returncode != 0:
            return None, None

        out = [line.strip() for line in p.stdout.splitlines() if line.strip()]
        if len(out) < 3:
            return None, None

        line = out[-1]
        struct = line.split()[0]

        score = 0.0
        if "(" in line and ")" in line:
            try:
                score = float(line[line.rfind("(")+1:line.rfind(")")])
            except ValueError:
                score = 0.0

        return struct, score
    finally:
        try:
            os.remove(fa_path)
        except OSError:
            pass


def fold_centroid(C: str):
    """centroid_fold -e McCaskill -> (dot-bracket, pseudo_energy)"""
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as f:
        f.write(">C\n")
        f.write(C + "\n")
        fa_path = f.name

    try:
        p = subprocess.run(
            ["centroid_fold", "-e", "McCaskill", fa_path],
            capture_output=True,
            text=True,
        )
        out = p.stdout.strip().splitlines()
        if len(out) < 2:
            return None, None

        line = out[-1].strip()
        struct = line.split()[0]

        # line is now similar to:
        # "....((...)).... (g=1,th=0.5,e=-3.57)"
        e = 0.0
        m = re.search(r"e=([\-0-9\.]+)", line)
        if m:
            try:
                e = float(m.group(1))
            except ValueError:
                e = 0.0

        return struct, e
    finally:
        os.remove(fa_path)


def Fold(C, mode="full"):
    """
    Unified wrapper:
      mode="cheap": RNAfold only (fastest)
      mode="medium": RNAfold + mxfold2
      mode="full": all three (RNAfold + centroid + mxfold2)
    Returns:
      - stru_C: [str, ...]          dot-bracket predicted by each model
      - E: [float, ...]             normalized "stability score"
    """
    stru_C_raw = []
    L_C = max(len(C), 1)

    # RNAfold always runs
    s, e = fold_rnafold(C)
    if s is not None and e is not None:
        stru_C_raw.append(("rnafold", s, e))

    if mode in ["medium", "full"]:
        s, e = fold_mxfold2(C)
        if s is not None and e is not None:
            stru_C_raw.append(("mxfold2", s, e))

    if mode == "full":
        s, e = fold_centroid(C)
        if s is not None and e is not None:
            stru_C_raw.append(("centroid", s, e))

    # Unify direction + normalize
    stru_C = []
    E = []

    for model, s, e in stru_C_raw:
        if model == "rnafold":
            score = -e
        elif model == "centroid":
            score = -e
        elif model == "mxfold2":
            score = e
        else:
            continue

        score = score / L_C
        score_norm = score / (1.0 + abs(score))

        stru_C.append(s)
        E.append(score_norm)

    return stru_C, E

