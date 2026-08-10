"""
Generate detailed charts comparing the top-10 individuals of GA vs Random Search.

Usage:
    python top10_charts.py

Inputs (expected in repository root):
    ga_results.csv
    rs_results.csv

Outputs:
    top10_metrics_by_rank.png   — Per-rank curves for key metrics
    top10_metric_means.png      — Bar chart of mean metrics across top-10
    top10_gc_vs_emean.png       — Scatter of GC vs Emean for top-10
"""

from pathlib import Path

import pandas as pd
import matplotlib

matplotlib.use("Agg")  # headless-safe backend
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
DATA_FILES = {"GA": ROOT / "ga_results.csv", "RS": ROOT / "rs_results.csv"}
METRIC_COLUMNS = [
    "Score",
    "Emean",
    "GC",
    "max_run",
    "Btop_dot",
    "Bbelow_dot",
    "Atop_dot",
    "Abelow_dot",
]
COLORS = {"GA": "steelblue", "RS": "tomato"}
MARKERS = {"GA": "o", "RS": "s"}


def load_top10(label: str, path: Path, n: int = 10) -> pd.DataFrame:
    """Load results CSV, sort by Score desc, and return top-n rows."""
    if not path.exists():
        raise FileNotFoundError(f"{path} not found. Please run compare.py first.")

    df = pd.read_csv(path)
    missing = [c for c in METRIC_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{path} missing required columns: {missing}")

    for col in METRIC_COLUMNS:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.sort_values("Score", ascending=False).head(n).reset_index(drop=True)
    df["rank"] = df.index + 1
    return df


def plot_metrics_by_rank(top_data: dict, output_path: Path) -> None:
    """Line plots by rank for key metrics."""
    metrics = [
        ("Score", "Score (↑ better)"),
        ("Emean", "Emean (↓ better)"),
        ("GC", "GC content"),
        ("Btop_dot", "Btop_dot"),
        ("Bbelow_dot", "Bbelow_dot"),
        ("Atop_dot", "Atop_dot"),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True)
    for (metric, title), ax in zip(metrics, axes.flat):
        for label, df in top_data.items():
            ax.plot(
                df["rank"],
                df[metric],
                marker=MARKERS[label],
                color=COLORS[label],
                label=label,
                linewidth=1.8,
            )
        ax.set_title(title)
        ax.set_xlabel("Rank (1 = best)")
        ax.set_ylabel(metric)
        ax.grid(True, linestyle="--", alpha=0.4)
    axes[0, 0].legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"Saved {output_path.name}")


def plot_metric_means(top_data: dict, output_path: Path) -> None:
    """Grouped bar chart of mean metrics across top-10."""
    metrics = [
        "Score",
        "Emean",
        "GC",
        "Btop_dot",
        "Bbelow_dot",
        "Atop_dot",
        "Abelow_dot",
        "max_run",
    ]
    means = pd.DataFrame({label: df[metrics].mean() for label, df in top_data.items()}).T
    means = means[metrics]

    fig, ax = plt.subplots(figsize=(12, 6))
    means.plot(kind="bar", ax=ax, color=["#4c78a8", "#f58518", "#54a24b", "#e45756",
                                         "#72b7b2", "#f2cf5b", "#b279a2", "#ff9da6"])
    ax.set_title("Top-10 Mean Metrics (GA vs RS)")
    ax.set_ylabel("Value")
    ax.set_xticklabels(means.index, rotation=0)
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"Saved {output_path.name}")


def plot_gc_vs_emean(top_data: dict, output_path: Path) -> None:
    """Scatter plot of GC vs Emean for the top-10 candidates."""
    fig, ax = plt.subplots(figsize=(7, 6))
    for label, df in top_data.items():
        ax.scatter(
            df["GC"],
            df["Emean"],
            s=60,
            color=COLORS[label],
            edgecolor="white",
            linewidth=0.8,
            alpha=0.9,
            label=label,
        )
    ax.set_title("Top-10 GC vs Emean")
    ax.set_xlabel("GC fraction")
    ax.set_ylabel("Emean (lower is better)")
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"Saved {output_path.name}")


def main():
    top_data = {
        label: load_top10(label, path)
        for label, path in DATA_FILES.items()
    }

    plot_metrics_by_rank(top_data, ROOT / "top10_metrics_by_rank.png")
    plot_metric_means(top_data, ROOT / "top10_metric_means.png")
    plot_gc_vs_emean(top_data, ROOT / "top10_gc_vs_emean.png")
    print("Charts generated. See PNG files in repository root.")


if __name__ == "__main__":
    main()
