"""
Generate the initial candidate population for the seekRNA GA.
For each valid insertion point i in sequence A, construct chimeric
sequence C = A[:i] + B + A[i:], score it with Phy_mark, and collect
all passing candidates into a CSV file.
"""

import pandas as pd
from config_loader import SEQUENCES, CONFIG
from utils.folding import Fold
from core.scoring import Phy_mark


def mk(A: str, B: str, i: int) -> str:
    """Create chimeric sequence C = A[:i] + B + A[i:]"""
    return A[:i] + B + A[i:]


def generate(fold_mode="cheap"):
    A = SEQUENCES["A"]
    B = SEQUENCES["B"]
    LA = len(A)

    # Fold reference sequences (full mode for references, only done once)
    print("Folding reference sequences A and B...")
    struct_ref_A, _ = Fold(A, mode="full")
    struct_ref_B, _ = Fold(B, mode="full")

    print(f"Scanning {LA - 3} insertion points (i = 2 .. {LA - 2})...")
    rows = []

    for i in range(2, LA - 1):
        C = mk(A, B, i)

        res = Phy_mark(C, struct_ref_A, struct_ref_B, i,
                       coords=None, fold_mode=fold_mode)

        if res is None:
            continue

        res["i"] = i
        rows.append(res)

    df = pd.DataFrame(rows)

    if df.empty:
        print("Warning: no candidates passed Phy_mark. CSV not written.")
        return df

    df = df.sort_values("Score", ascending=False).reset_index(drop=True)

    out_path = CONFIG.get("initial_population_csv", "initial_candidates.csv")
    df.to_csv(out_path, index=False)
    print(f"Generated {len(df)} initial candidates -> {out_path}")

    return df


if __name__ == "__main__":
    generate(fold_mode="cheap")
