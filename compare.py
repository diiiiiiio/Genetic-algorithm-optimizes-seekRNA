"""compare.py — Run GA vs Random Search and produce comparison metrics and plots.

Usage:
    python compare.py

Outputs:
    ga_history.csv         — Per-round best_score / average_score for GA
    rs_history.csv         — Per-round best_score / average_score for Random Search
    comparison_plot.png    — Side-by-side convergence curves
    ga_results.csv         — Final GA population (top candidates)
    rs_results.csv         — Final RS best pool (top candidates)
"""

import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # non-interactive backend; safe in headless environments
import matplotlib.pyplot as plt

from config_loader import SEQUENCES, GA_PARAMS, MOTIF_LENGTHS, CONFIG
from utils.folding import Fold
from core.ga_engine import GAEngine
from core.rs_engine import RandomSearchEngine

A_MOTIF_LEN = MOTIF_LENGTHS["A_MOTIF_LEN"]
B_MOTIF_LEN = MOTIF_LENGTHS["B_MOTIF_LEN"]


def add_motif_columns_if_missing(df):
    """Ensure all motif coordinate columns exist in the DataFrame."""
    motif_cols = ["a1m", "a2m", "A1m", "A2m", "b1m", "b2m", "B1m", "B2m"]
    if all(col in df.columns for col in motif_cols):
        return df

    df = df.copy()
    df["a1m"] = df["a2"] - A_MOTIF_LEN
    df["a2m"] = df["a2"]
    df["A1m"] = df["A2"] - A_MOTIF_LEN
    df["A2m"] = df["A2"]
    df["b1m"] = df["b2"] - B_MOTIF_LEN
    df["b2m"] = df["b2"]
    df["B1m"] = df["B2"] - B_MOTIF_LEN
    df["B2m"] = df["B2"]
    print("Added missing motif columns to DataFrame")
    return df


def load_initial_population():
    """Load and prepare the initial candidate population."""
    input_csv = CONFIG.get("initial_population_csv", "initial_candidates.csv")
    if not os.path.exists(input_csv):
        raise FileNotFoundError(
            f"{input_csv} not found. "
            "Run 'python generate_initial_population.py' first."
        )

    df = pd.read_csv(input_csv)
    df = add_motif_columns_if_missing(df)

    int_cols = [
        "i", "a1", "a2", "A1", "A2", "b1", "b2", "B1", "B2",
        "a1m", "a2m", "A1m", "A2m", "b1m", "b2m", "B1m", "B2m",
    ]
    for col in int_cols:
        if col in df.columns:
            df[col] = df[col].astype(int)
    return df


def plot_comparison(ga_history, rs_history, output_path="comparison_plot.png"):
    """Generate side-by-side convergence plots for GA vs Random Search."""
    ga_df = pd.DataFrame(ga_history)
    rs_df = pd.DataFrame(rs_history)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Genetic Algorithm vs. Random Search — seekRNA Optimization", fontsize=13)

    # --- Best Score Convergence ---
    ax = axes[0]
    ax.plot(ga_df["round"], ga_df["best_score"], marker="o", markersize=3,
            color="steelblue", label="GA")
    ax.plot(rs_df["round"], rs_df["best_score"], marker="s", markersize=3,
            color="tomato", label="Random Search")
    ax.set_title("Best Score Convergence")
    ax.set_xlabel("Round")
    ax.set_ylabel("Best Score")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.5)

    # --- Average Population Score ---
    ax = axes[1]
    ax.plot(ga_df["round"], ga_df["average_score"], marker="o", markersize=3,
            color="steelblue", label="GA")
    ax.plot(rs_df["round"], rs_df["average_score"], marker="s", markersize=3,
            color="tomato", label="Random Search")
    ax.set_title("Average Score per Round")
    ax.set_xlabel("Round")
    ax.set_ylabel("Average Score")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"Comparison plot saved to {output_path}")


def main():
    # 1) Load shared initial population
    print("Loading initial population...")
    initial_pop = load_initial_population()
    print(f"  {len(initial_pop)} candidates loaded.\n")

    # 2) Pre-fold reference sequences (shared between GA and RS)
    print("Folding reference sequences A and B...")
    struct_A, _ = Fold(SEQUENCES["A"], mode="full")
    struct_B, _ = Fold(SEQUENCES["B"], mode="full")
    print()

    # 3) Run Genetic Algorithm
    print("=" * 60)
    print("Running Genetic Algorithm (GA)...")
    print("=" * 60)
    ga_engine = GAEngine(initial_pop, struct_A, struct_B)
    ga_final_pop, ga_history = ga_engine.run()
    ga_final_pop.to_csv("ga_results.csv", index=False)
    pd.DataFrame(ga_history).to_csv("ga_history.csv", index=False)
    print(f"GA complete. Results → ga_results.csv | History → ga_history.csv\n")

    # 4) Run Random Search
    print("=" * 60)
    print("Running Random Search (RS)...")
    print("=" * 60)
    rs_engine = RandomSearchEngine(initial_pop, struct_A, struct_B)
    rs_best_pool, rs_history = rs_engine.run()
    rs_best_pool.to_csv("rs_results.csv", index=False)
    pd.DataFrame(rs_history).to_csv("rs_history.csv", index=False)
    print(f"RS complete. Results → rs_results.csv | History → rs_history.csv\n")

    # 5) Generate comparison plot
    plot_comparison(ga_history, rs_history, output_path="comparison_plot.png")

    # 6) Summary
    ga_best = ga_history[-1]["best_score"] if ga_history else float("nan")
    rs_best = rs_history[-1]["best_score"] if rs_history else float("nan")
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  GA  final best score : {ga_best:.4f}")
    print(f"  RS  final best score : {rs_best:.4f}")
    if rs_best and rs_best != 0 and not (rs_best != rs_best):  # guard against 0 and NaN
        print(f"  GA advantage         : {(ga_best - rs_best):.4f} ({(ga_best / rs_best - 1) * 100:.1f}% better)")
    print("=" * 60)


if __name__ == "__main__":
    main()
