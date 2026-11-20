"""Extract numerical summaries from a completed, paired search benchmark."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def correlation(lengths, scores):
    """Spearman's rho with average ranks for ties."""
    x = pd.Series(lengths).rank(method="average")
    y = pd.Series(scores).rank(method="average")
    return float(x.corr(y))


def expected_maximum(values, sample_size):
    """Expected maximum in a uniform sample drawn without replacement."""
    descending = np.sort(np.asarray(values, dtype=float))[::-1]
    n = len(descending)
    if not 1 <= sample_size <= n:
        raise ValueError("invalid sample size")
    # P(the first selected value is at rank j), traversing descending ranks.
    probability = sample_size / n
    expectation = 0.0
    for offset in range(n - sample_size + 1):
        expectation += probability * descending[offset]
        remaining = n - offset - 1
        if remaining >= sample_size:
            probability *= (remaining - sample_size + 1) / remaining
    return expectation


def read_run(run_dir, method, seed, budget, baseline):
    summary = json.loads((run_dir / "summary.json").read_text())
    if (summary["method"], summary["seed"], summary["complete"],
        summary["budget"], summary["used"]) != (method, seed, True, budget, budget):
        raise ValueError(f"Incomplete or unexpected run at {run_dir}")
    columns = ["anon_id", "status", "budget_used", "length", "score", "motif_contained"]
    all_events = pd.read_csv(run_dir / "events.csv", usecols=columns)
    spent = all_events[all_events.status.isin(["scored", "failed"])].copy()
    if not np.array_equal(spent.budget_used.to_numpy(), np.arange(1, budget + 1)):
        raise ValueError(f"Budget is discontinuous at {run_dir}")
    if spent.anon_id.duplicated().any():
        raise ValueError(f"Repeated charged record at {run_dir}")
    scored = spent[spent.status == "scored"].copy()
    flags = scored.motif_contained.astype(str).str.lower().map(
        {"0": False, "0.0": False, "false": False, "1": True, "1.0": True, "true": True}
    )
    if flags.isna().any() or scored.score.isna().any():
        raise ValueError(f"Missing numeric data at {run_dir}")
    scored["contained"] = flags.astype(bool)
    winner = scored.loc[scored.anon_id == summary["best_id"]]
    best = max(baseline, scored.score.max())
    if (len(winner) != 1 or abs(winner.score.iloc[0] - best) > 1e-8
            or abs(summary["best"] - best) > 1e-8
            or summary["scored"] != len(scored)):
        raise ValueError(f"Terminal record is inconsistent at {run_dir}")
    contained_best = max(baseline, scored.loc[scored.contained, "score"].max())
    outcomes = {
        "gain": best - baseline, "contained_gain": contained_best - baseline,
        "scored": len(scored), "scored_contained": int(scored.contained.sum()),
        "winner_contained": int(winner.contained.iloc[0]),
        "length_score_spearman": correlation(scored.length, scored.score),
    }
    gain_trace = np.maximum.accumulate(
        np.r_[baseline, spent.score.fillna(baseline).to_numpy()]
    ) - baseline
    return outcomes, gain_trace, scored[["length", "score"]]


def extract(archive, output, protocol):
    seeds, methods = protocol["seeds"], protocol["methods"]
    budget, baseline = protocol["budget"], protocol["starting_best"]
    stride = protocol["checkpoint_every"]
    for method in methods:
        found = {int(p.name[5:]) for p in (archive / method).glob("seed_*") if p.is_dir()}
        if found != set(seeds):
            raise ValueError(f"The {method} archive differs from the declared seed list")
    checks = np.arange(0, budget + 1, stride)
    if checks[-1] != budget:
        raise ValueError("Terminal checkpoint missing")
    per_seed, strata = [], []
    curves = {m: [] for m in methods}
    pooled = {m: [] for m in methods}
    for seed in seeds:
        row, by_length = {"seed": seed}, {}
        for method in methods:
            outcome, trace, scored = read_run(
                archive / method / f"seed_{seed}", method, seed, budget, baseline
            )
            row.update({f"{method}_{k}": v for k, v in outcome.items()})
            curves[method].append(trace[checks])
            pooled[method].append(scored)
            by_length[method] = {k: v.score.to_numpy() for k, v in scored.groupby("length")}
        per_seed.append(row)
        for length in sorted(by_length["ga"].keys() & by_length["crs"].keys()):
            groups = [by_length[m][length] for m in methods]
            size = min(len(g) for g in groups)
            entry = {"seed": seed, "length": int(length), "comparison_count": size}
            for method, values in zip(methods, groups):
                entry[f"{method}_n"] = len(values)
                entry[f"{method}_max"] = float(np.max(values))
                entry[f"{method}_expected_max"] = expected_maximum(values, size)
            strata.append(entry)
    checkpoints, correlations = {"evaluations": checks}, []
    for method in methods:
        for label, row in zip(("q25", "median", "q75"),
                              np.quantile(curves[method], [0.25, 0.50, 0.75], axis=0)):
            checkpoints[f"{method}_{label}"] = row
        records = pd.concat(pooled[method], ignore_index=True)
        correlations.append({"method": method, "scored": len(records),
                             "spearman": correlation(records.length, records.score)})
    output.mkdir(parents=True, exist_ok=True)
    for name, rows in [("seed_outcomes", per_seed), ("budget_checkpoints", checkpoints),
                       ("length_strata", strata), ("pooled_length_correlation", correlations)]:
        pd.DataFrame(rows).to_csv(output / f"{name}.csv", index=False, float_format="%.12g")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path, help="Formal event archive with ga/ and crs/ directories")
    parser.add_argument("output", type=Path, help="Output directory for numerical summaries")
    parser.add_argument("--protocol", type=Path, default=Path(__file__).with_name("protocol.json"))
    args = parser.parse_args()
    extract(args.archive, args.output, json.loads(args.protocol.read_text()))
