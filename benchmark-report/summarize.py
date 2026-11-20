"""Reproduce the numerical benchmark results from the adjacent summary tables."""

import json
from pathlib import Path

import numpy as np
import pandas as pd


def interval(values, bootstrap):
    values = np.asarray(values, dtype=float)
    indices = np.random.default_rng(bootstrap["seed"]).integers(
        0, len(values), size=(bootstrap["resamples"], len(values))
    )
    lo, hi = np.quantile(values[indices].mean(axis=1), [0.025, 0.975])
    return dict(mean=float(values.mean()), low=float(lo), high=float(hi))


def paired_endpoint(rows, field, bootstrap):
    ga, crs = rows[f"ga_{field}"], rows[f"crs_{field}"]
    difference = ga - crs
    return dict(ga_mean=float(ga.mean()), crs_mean=float(crs.mean()),
                ga_higher=int((difference > 0).sum()), **interval(difference, bootstrap))


def summarize(directory):
    protocol = json.loads((directory / "protocol.json").read_text())
    bootstrap = protocol["bootstrap"]
    runs = pd.read_csv(directory / "seed_outcomes.csv")
    strata = pd.read_csv(directory / "length_strata.csv")
    correlations = pd.read_csv(directory / "pooled_length_correlation.csv").set_index("method")
    components = pd.read_csv(directory / "winner_components.csv").set_index("metric")
    seeds = protocol["seeds"]
    if (runs.seed.tolist() != seeds or sorted(strata.seed.unique()) != sorted(seeds)
            or strata.duplicated(["seed", "length"]).any()):
        raise ValueError("Unexpected or incomplete seed/length coverage")
    if not np.array_equal(strata.comparison_count, strata[["ga_n", "crs_n"]].min(axis=1)):
        raise ValueError("Wrong weights for shared lengths")
    length_results = {}
    for field in ("max", "expected_max"):
        differences = strata[f"ga_{field}"] - strata[f"crs_{field}"]
        weighted = differences * strata.comparison_count
        per_seed = weighted.groupby(strata.seed).sum() / strata.groupby("seed").comparison_count.sum()
        length_results[field] = interval(per_seed.reindex(seeds), bootstrap)
    primary = paired_endpoint(runs, "gain", bootstrap)
    if not np.isclose(components.score_contribution.sum(), primary["mean"], atol=1e-8, rtol=0):
        raise ValueError("Component means do not account for the endpoint difference")
    scored = int(runs[["ga_scored", "crs_scored"]].sum().sum())
    requests = protocol["budget"] * len(runs) * len(protocol["methods"])
    return {
        "seed_pairs": len(runs), "requests": requests, "scored": scored, "failed": requests - scored,
        "primary": primary, "coordinate_check": paired_endpoint(runs, "contained_gain", bootstrap),
        "winners_passing_by_method": {
            m: int(runs[f"{m}_winner_contained"].sum()) for m in protocol["methods"]
        },
        "pairs_with_both_winners_passing": int(
            (runs.ga_winner_contained.astype(bool) & runs.crs_winner_contained.astype(bool)).sum()
        ),
        "length": {
            "shared_strata": len(strata), "unequal_counts": int((strata.ga_n != strata.crs_n).sum()),
            "pooled_spearman": correlations.spearman.to_dict(), **length_results,
        },
        "winner_weighted_score_contributions": components.score_contribution.dropna().to_dict(),
    }


if __name__ == "__main__":
    print(json.dumps(summarize(Path(__file__).resolve().parent), indent=2))
