import random
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
from core.scoring import cheap_filter, cheap_score, Phy_mark
from core.operators import mutate_in_targets, crossover
from core.sampling import weighted_sample_df
from config_loader import GA_PARAMS

def _evaluate_candidate_worker(args):
    """Worker function for parallel evaluation of Phy_mark"""
    seq, struct_ref_A, struct_ref_B, i, coords, fold_mode = args
    return Phy_mark(seq, struct_ref_A, struct_ref_B, i, coords, fold_mode=fold_mode)

class GAEngine:
    def __init__(self, initial_pop_df, struct_ref_A, struct_ref_B):
        self.pop_df = initial_pop_df.copy()
        self.struct_ref_A = struct_ref_A
        self.struct_ref_B = struct_ref_B
        
        # Load params
        self.rounds = GA_PARAMS["rounds"]
        self.pool_size = GA_PARAMS["pool_size"]
        self.top_percent = GA_PARAMS["top_percent"]
        self.elite_size = GA_PARAMS["elite_size"]
        self.target_pop_size = GA_PARAMS["target_pop_size"]
        self.fold_mode = GA_PARAMS["fold_mode"]
        self.n_workers = GA_PARAMS["n_workers"] or max(1, cpu_count() - 1)
        self.p_crossover = GA_PARAMS["p_crossover"]
        self.n_mut_positions = GA_PARAMS["n_mut_positions"]

    def run_round(self, round_idx):
        print(f"=== Round {round_idx + 1}/{self.rounds} ===")
        
        # 1) Sort and prepare parents
        self.pop_df = self.pop_df.sort_values("Score", ascending=False).reset_index(drop=True)
        parents_df = self.pop_df.head(self.target_pop_size).copy()
        
        # Calculate cheap weights
        cheap_w = []
        for _, r in parents_df.iterrows():
            coords_r = {k: int(r[k]) for k in ["a1","a2","A1","A2","b1","b2","B1","B2","a1m","a2m","A1m","A2m","b1m","b2m","B1m","B2m"]}
            try:
                w = cheap_score(r["seq"], coords_r)
            except:
                w = 0.0
            cheap_w.append(max(0.0, w) + 0.05)
        parents_df["cheap_w"] = cheap_w

        # 2) Tier 1: Generate candidate pool with cheap filtering
        candidate_pool = []
        attempts = 0
        max_attempts = self.pool_size * 10
        
        while len(candidate_pool) < self.pool_size and attempts < max_attempts:
            attempts += 1
            if len(self.pop_df) >= 2 and random.random() < self.p_crossover:
                parents = weighted_sample_df(parents_df, 'cheap_w', k=2)
                res = crossover(parents[0], parents[1])
                if res is None: continue
                child_seq, child_coords, i_child = res
            else:
                p = weighted_sample_df(parents_df, 'cheap_w', k=1)[0]
                coords_p = {k: int(p[k]) for k in ["a1","a2","A1","A2","b1","b2","B1","B2","a1m","a2m","A1m","A2m","b1m","b2m","B1m","B2m"]}
                child_seq, child_coords = mutate_in_targets(p["seq"], coords_p, n_events=self.n_mut_positions)
                i_child = int(p["i"])

            if cheap_filter(child_seq, child_coords):
                candidate_pool.append({
                    "seq": child_seq, "coords": child_coords,
                    "i": i_child, "cheap_score": cheap_score(child_seq, child_coords)
                })

        print(f"  Generated {len(candidate_pool)} candidates (cheap filter pass)")

        # 3) Tier 2: Expensive evaluation
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

        # 4) Population update
        if new_individuals:
            self.pop_df = pd.concat([self.pop_df, pd.DataFrame(new_individuals)], ignore_index=True)
        
        self.pop_df = self.pop_df.drop_duplicates(subset=["seq"]).sort_values("Score", ascending=False)
        self.pop_df = self.pop_df.head(self.target_pop_size).reset_index(drop=True)
        
        print(f"  Population size: {len(self.pop_df)}, Best score: {self.pop_df['Score'].iloc[0]:.4f}\n")

    def run(self):
        for r in range(self.rounds):
            self.run_round(r)
        return self.pop_df

