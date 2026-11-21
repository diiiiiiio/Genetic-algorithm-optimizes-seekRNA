"""Run GA and cumulative random search under one RNAfold budget.

The initial table was scored with RNAfold. References, the shared start, and
both methods use that same mode. Phy_mark still reads three reference slots, so
the single RNAfold structure is repeated three times. Scores and weights inside
Phy_mark are unchanged.
"""

import argparse
import csv
import json
import os
import time
from multiprocessing import Pool
from pathlib import Path

from config_loader import CONFIG, GA_PARAMS, SEQUENCES
from core.fair_search import COORD_KEYS, EVENT_FIELDS, copy_individual, run_method, score_record
from core.scoring import Phy_mark
from utils.folding import Fold

ROOT = Path(__file__).resolve().parents[1]


def ensure_tool_path():
    bindir = str(ROOT / ".venv" / "bin")
    os.environ["PATH"] = bindir + os.pathsep + os.environ.get("PATH", "")
    for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ.setdefault(name, "1")


def pad_reference(structures):
    if not structures:
        raise RuntimeError("reference folding returned no structure")
    if len(structures) >= 3:
        return list(structures[:3])
    return [structures[0], structures[0], structures[0]]


def load_stored_initials():
    import pandas as pd

    frame = pd.read_csv(ROOT / CONFIG.get("initial_population_csv", "initial_candidates.csv"))
    rows = []
    for record in frame.to_dict("records"):
        rows.append({
            "seq": record["seq"],
            "coords": {key: int(record[key]) for key in COORD_KEYS},
            "i": int(record["i"]),
            "stored_score": float(record["Score"]),
            "depth": 0,
        })
    if len({row["seq"] for row in rows}) != len(rows):
        raise RuntimeError("initial table contains a repeated sequence")
    return rows


def rescore_initials(stored, ref_a, ref_b, fold_mode):
    rescored = []
    gaps = []
    for row in stored:
        result = Phy_mark(row["seq"], ref_a, ref_b, row["i"], row["coords"], fold_mode=fold_mode)
        if result is None:
            raise RuntimeError("a shared start candidate failed the common scoring configuration")
        score = float(result["Score"])
        gaps.append(abs(score - row["stored_score"]))
        item = copy_individual({**row, "score": score})
        rescored.append(item)
    return rescored, max(gaps)


def write_events(path, rows_iter):
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = path.open("w", newline="")
    writer = csv.DictWriter(handle, fieldnames=EVENT_FIELDS)
    writer.writeheader()
    handle.flush()
    return handle, writer


def execute(output, seeds, budget, workers, pool_size, top_percent, target_pop_size, p_crossover, n_mut, fold_mode):
    output = Path(output)
    if output.exists():
        raise SystemExit(f"{output} already exists")
    output.mkdir(parents=True)
    ensure_tool_path()
    started = time.time()
    ref_a, _ = Fold(SEQUENCES["A"], mode=fold_mode)
    ref_b, _ = Fold(SEQUENCES["B"], mode=fold_mode)
    ref_a, ref_b = pad_reference(ref_a), pad_reference(ref_b)
    stored = load_stored_initials()
    initials, max_gap = rescore_initials(stored, ref_a, ref_b, fold_mode)
    if fold_mode == "cheap" and max_gap > 1e-6:
        raise RuntimeError(f"rescored starts differ from the initial table by {max_gap}")
    start_best = max(item["score"] for item in initials)
    manifest = {
        "fold_mode": fold_mode,
        "reference_policy": "same mode as candidates; RNAfold structure repeated to fill three scorer slots" if fold_mode == "cheap" else "structures returned by the requested mode",
        "scoring": "Phy_mark weights unchanged",
        "budget_unit": "first full score of a sequence, including failures",
        "methods": ["ga", "crs"],
        "seeds": list(seeds),
        "budget": budget,
        "workers": workers,
        "pool_size": pool_size,
        "top_percent": top_percent,
        "target_pop_size": target_pop_size,
        "p_crossover": p_crossover,
        "n_mut": n_mut,
        "start_best": start_best,
        "start_rescore_max_abs_gap": max_gap,
        "primary_metric": "best score at budget minus shared start best",
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2))
    with (output / "starts.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["anon_length", "score", "insertion", "depth"])
        writer.writeheader()
        for item in initials:
            writer.writerow({
                "anon_length": len(item["seq"]),
                "score": item["score"],
                "insertion": item["i"],
                "depth": item["depth"],
            })

    def progress(method, seed, round_index, used, best):
        print(f"{method} seed {seed} round {round_index} used {used}/{budget} best {best:.6f}", flush=True)

    pool = Pool(workers) if workers > 1 else None
    try:
        def score_many(batch):
            args = [(item["seq"], ref_a, ref_b, int(item["i"]), dict(item["coords"]), fold_mode) for item in batch]
            if pool is None:
                return [score_record(arg) for arg in args]
            return pool.map(score_record, args, chunksize=1)

        for seed in seeds:
            for method in ("ga", "crs"):
                destination = output / method / f"seed_{seed}"
                handle, writer = write_events(destination / "events.csv", None)
                try:
                    def emit(row, writer=writer, handle=handle):
                        writer.writerow(row)
                        handle.flush()

                    begun = time.perf_counter()
                    summary = run_method(
                        method, initials, score_many, seed=seed, budget=budget,
                        pool_size=pool_size, top_percent=top_percent,
                        target_pop_size=target_pop_size, p_crossover=p_crossover,
                        n_mut=n_mut, fold_mode=fold_mode, event_writer=emit,
                        progress=progress,
                    )
                    summary["seconds"] = time.perf_counter() - begun
                    summary.pop("events", None)
                    (destination / "summary.json").write_text(json.dumps(summary, indent=2))
                    print(f"{method} seed {seed} gain {summary['gain']:.6f} stop {summary['stop']} seconds {summary['seconds']:.1f}", flush=True)
                finally:
                    handle.close()
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    (output / "elapsed-seconds.txt").write_text(str(time.time() - started))
    print("FAIR RUN COMPLETE", flush=True)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--seeds", required=True, help="comma-separated integers")
    parser.add_argument("--budget", type=int, required=True)
    parser.add_argument("--workers", type=int, default=max(1, min(3, os.cpu_count() - 1)))
    parser.add_argument("--pool-size", type=int, default=GA_PARAMS["pool_size"])
    parser.add_argument("--top-percent", type=float, default=GA_PARAMS["top_percent"])
    parser.add_argument("--target-pop-size", type=int, default=GA_PARAMS["target_pop_size"])
    parser.add_argument("--p-crossover", type=float, default=GA_PARAMS["p_crossover"])
    parser.add_argument("--n-mut", type=int, default=GA_PARAMS["n_mut_positions"])
    parser.add_argument("--fold-mode", default="cheap")
    return parser.parse_args()


def main():
    args = parse_args()
    seeds = [int(part) for part in args.seeds.split(",") if part]
    execute(
        args.output, seeds, args.budget, args.workers, args.pool_size, args.top_percent,
        args.target_pop_size, args.p_crossover, args.n_mut, args.fold_mode,
    )


if __name__ == "__main__":
    main()
