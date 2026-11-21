"""Summarize paired gains, same-length scores, and a seed bootstrap interval."""

import csv
import json
import math
from pathlib import Path

import numpy as np


def load_summaries(root):
    found = {}
    for path in Path(root).glob("*/*/summary.json"):
        row = json.loads(path.read_text())
        found.setdefault(row["method"], {})[int(row["seed"])] = row
    return found


def load_scored(root):
    rows = []
    for path in Path(root).glob("*/*/events.csv"):
        with path.open(newline="") as handle:
            for record in csv.DictReader(handle):
                if record["status"] == "scored":
                    rows.append({
                        "method": record["method"],
                        "seed": int(record["seed"]),
                        "length": int(record["length"]),
                        "score": float(record["score"]),
                        "depth": int(record["edit_depth"]),
                    })
    return rows


def bootstrap_mean(values, draws=10000):
    sample = np.asarray(values, dtype=float)
    if len(sample) < 2:
        return None
    generator = np.random.default_rng(20261121)
    picked = generator.integers(0, len(sample), size=(draws, len(sample)))
    means = sample[picked].mean(axis=1)
    low, high = np.quantile(means, [0.025, 0.975])
    return {"mean": float(sample.mean()), "low": float(low), "high": float(high), "draws": draws}


def length_gap(rows, seeds):
    per_seed = []
    for seed in seeds:
        grouped = {}
        for method in ("ga", "crs"):
            selected = [row for row in rows if row["seed"] == seed and row["method"] == method]
            for row in selected:
                grouped.setdefault(method, {}).setdefault(row["length"], []).append(row["score"])
        common = set(grouped.get("ga", {})) & set(grouped.get("crs", {}))
        strata = []
        for length in sorted(common):
            ga_scores = grouped["ga"][length]
            crs_scores = grouped["crs"][length]
            weight = min(len(ga_scores), len(crs_scores))
            strata.append({
                "length": length,
                "ga_n": len(ga_scores),
                "crs_n": len(crs_scores),
                "ga_max": max(ga_scores),
                "crs_max": max(crs_scores),
                "difference": max(ga_scores) - max(crs_scores),
                "weight": weight,
            })
        weight = sum(item["weight"] for item in strata)
        weighted = math.fsum(item["weight"] * item["difference"] for item in strata) / weight if weight else None
        per_seed.append({"seed": seed, "strata": len(strata), "weighted_max_gap": weighted})
    values = [item["weighted_max_gap"] for item in per_seed if item["weighted_max_gap"] is not None]
    return {"per_seed": per_seed, "across_seeds": bootstrap_mean(values) if values else None}


def rank_correlation(rows, method):
    selected = [row for row in rows if row["method"] == method]
    if len(selected) < 3:
        return None
    length = np.array([row["length"] for row in selected], dtype=float)
    score = np.array([row["score"] for row in selected], dtype=float)
    if len(set(length)) < 2 or len(set(score)) < 2:
        return None
    return float(np.corrcoef(np.argsort(np.argsort(length)), np.argsort(np.argsort(score)))[0, 1])


def main(root):
    summaries = load_summaries(root)
    seeds = sorted(set(summaries.get("ga", {})) & set(summaries.get("crs", {})))
    if not seeds:
        raise SystemExit("no paired summaries found")
    gains = []
    details = []
    for seed in seeds:
        ga = summaries["ga"][seed]
        crs = summaries["crs"][seed]
        difference = ga["gain"] - crs["gain"]
        gains.append(difference)
        details.append({
            "seed": seed,
            "ga_gain": ga["gain"],
            "crs_gain": crs["gain"],
            "difference": difference,
            "ga_best": ga["best"],
            "crs_best": crs["best"],
            "ga_best_length": ga["best_length"],
            "crs_best_length": crs["best_length"],
            "ga_best_depth": ga["best_depth"],
            "crs_best_depth": crs["best_depth"],
            "ga_used": ga["used"],
            "crs_used": crs["used"],
            "ga_complete": ga["complete"],
            "crs_complete": crs["complete"],
            "ga_motif_contained": ga["motif_contained"],
            "ga_motif_checked": ga["motif_checked"],
            "crs_motif_contained": crs["motif_contained"],
            "crs_motif_checked": crs["motif_checked"],
            "ga_mean_depth": ga["mean_scored_depth"],
            "crs_mean_depth": crs["mean_scored_depth"],
        })
    scored = load_scored(root)
    output = {
        "seeds": seeds,
        "paired_gain_difference": bootstrap_mean(gains),
        "runs": details,
        "same_length_max_gap": length_gap(scored, seeds),
        "length_score_correlation": {"ga": rank_correlation(scored, "ga"), "crs": rank_correlation(scored, "crs")},
        "scored_rows": len(scored),
    }
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    import sys
    main(sys.argv[1])
