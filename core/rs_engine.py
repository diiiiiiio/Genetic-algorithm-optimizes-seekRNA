import random
import pandas as pd
from multiprocessing import Pool, cpu_count
from core.scoring import cheap_filter, cheap_score, Phy_mark
from core.operators import mutate_in_targets
from config_loader import GA_PARAMS


def _evaluate_candidate_worker(args):
    """Worker function for parallel evaluation of Phy_mark"""
    seq, struct_ref_A, struct_ref_B, i, coords, fold_mode = args
    return Phy_mark(seq, struct_ref_A, struct_ref_B, i, coords, fold_mode=fold_mode)


class RandomSearchEngine:
    """Random Search baseline for comparison against GAEngine.

    Key differences from GA:
    - Parents are always sampled *uniformly* from the fixed initial population
      (no score-based selection pressure).
    - Crossover is disabled; every child is produced by pure mutation only.
    - The same ``mutate_in_targets`` function is used to guarantee coordinate
      validity, keeping the search space identical to GA.
    - Discovered high-scoring individuals are tracked in ``self.best_pool`` for
      reporting, but they are *never* used as parents in subsequent rounds.
    - The same computational budget (pool_size, top_percent, n_workers) is used
      as GA so that comparisons are fair.
    """

    def __init__(self, initial_pop_df, struct_ref_A, struct_ref_B):
        self.initial_pop_df = initial_pop_df.copy()
        self.struct_ref_A = struct_ref_A
        self.struct_ref_B = struct_ref_B
        self.history = []

        # Best individuals found across all rounds (for reporting only)
        self.best_pool = initial_pop_df.copy()

        # Reuse the same GA parameters to keep the budget identical
        self.rounds = GA_PARAMS["rounds"]
        self.pool_size = GA_PARAMS["pool_size"]
        self.top_percent = GA_PARAMS["top_percent"]
        self.target_pop_size = GA_PARAMS["target_pop_size"]
        self.fold_mode = GA_PARAMS["fold_mode"]
        self.n_workers = GA_PARAMS["n_workers"] or max(1, cpu_count() - 1)
        self.n_mut_positions = GA_PARAMS["n_mut_positions"]

    def run_round(self, round_idx):
        print(f"=== RS Round {round_idx + 1}/{self.rounds} ===")

        # 1) Uniform sampling from the *fixed* initial population (no selection pressure)
        initial_records = self.initial_pop_df.to_dict("records")

        # 2) Tier 1: Generate candidate pool with cheap filtering
        candidate_pool = []
        attempts = 0
        max_attempts = self.pool_size * 10

        while len(candidate_pool) < self.pool_size and attempts < max_attempts:
            attempts += 1
            # Uniform random parent from initial population — no score weighting.
            # Uses the same mutate_in_targets as GA to guarantee coordinate validity.
            p = random.choice(initial_records)
            coords_p = {
                k: int(p[k])
                for k in ["a1", "a2", "A1", "A2", "b1", "b2", "B1", "B2",
                          "a1m", "a2m", "A1m", "A2m", "b1m", "b2m", "B1m", "B2m"]
            }
            child_seq, child_coords = mutate_in_targets(
                p["seq"], coords_p, n_events=self.n_mut_positions
            )
            i_child = int(p["i"])

            if cheap_filter(child_seq, child_coords):
                candidate_pool.append({
                    "seq": child_seq,
                    "coords": child_coords,
                    "i": i_child,
                    "cheap_score": cheap_score(child_seq, child_coords),
                })

        print(f"  Generated {len(candidate_pool)} candidates (cheap filter pass)")

        # 3) Tier 2: Expensive evaluation on the top cheap-scored candidates
        candidate_pool.sort(key=lambda x: x["cheap_score"], reverse=True)
        n_top = max(1, int(len(candidate_pool) * self.top_percent))
        top_candidates = candidate_pool[:n_top]

        print(f"  Evaluating top {n_top} candidates with Phy_mark (mode={self.fold_mode})...")

        eval_args = [
            (c["seq"], self.struct_ref_A, self.struct_ref_B, c["i"], c["coords"], self.fold_mode)
            for c in top_candidates
        ]

        with Pool(self.n_workers) as pool:
            results = pool.map(_evaluate_candidate_worker, eval_args)

        new_individuals = []
        for cand, res in zip(top_candidates, results):
            if res is not None:
                res["i"] = cand["i"]
                new_individuals.append(res)

        print(f"  {len(new_individuals)} candidates passed Phy_mark")

        # 4) Update best_pool for tracking — NOT used as parents next round
        if new_individuals:
            self.best_pool = pd.concat(
                [self.best_pool, pd.DataFrame(new_individuals)], ignore_index=True
            )
        self.best_pool = (
            self.best_pool.drop_duplicates(subset=["seq"])
            .sort_values("Score", ascending=False)
            .head(self.target_pop_size)
            .reset_index(drop=True)
        )

        # 5) Compute per-round metrics
        round_scores = [ind["Score"] for ind in new_individuals]
        if round_scores:
            round_best = max(round_scores)
            round_avg = sum(round_scores) / len(round_scores)
        elif len(self.best_pool) > 0:
            round_best = float(self.best_pool["Score"].iloc[0])
            round_avg = float(self.best_pool["Score"].mean())
        else:
            round_best = 0.0
            round_avg = 0.0

        # best_score reported as the cumulative best ever seen
        cumulative_best = float(self.best_pool["Score"].iloc[0]) if len(self.best_pool) > 0 else round_best
        self.history.append({
            "round": round_idx + 1,
            "best_score": cumulative_best,
            "average_score": round_avg,
        })

        print(
            f"  Best pool size: {len(self.best_pool)}, "
            f"Cumulative best: {cumulative_best:.4f}, "
            f"Round avg: {round_avg:.4f}\n"
        )

    def run(self):
        self.history = []
        for r in range(self.rounds):
            self.run_round(r)
        return self.best_pool, self.history
