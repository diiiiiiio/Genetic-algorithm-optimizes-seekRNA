"""Budget, replay, and edit-depth checks with a stand-in scorer."""

import unittest
from pathlib import Path

import pandas as pd

from core.fair_search import run_method

ROOT = Path(__file__).resolve().parents[1]
KEYS = ["a1", "a2", "A1", "A2", "b1", "b2", "B1", "B2", "a1m", "a2m", "A1m", "A2m", "b1m", "b2m", "B1m", "B2m"]


def initials():
    frame = pd.read_csv(ROOT / "initial_candidates.csv")
    rows = []
    for record in frame.to_dict("records"):
        rows.append({
            "seq": record["seq"],
            "coords": {key: int(record[key]) for key in KEYS},
            "i": int(record["i"]),
            "score": float(record["Score"]),
            "depth": 0,
        })
    return rows


def accepting(batch):
    return [({"Score": 10.0 + item["depth"]}, 0.0) for item in batch]


def rejecting(batch):
    return [(None, 0.0) for _ in batch]


def run(method, scorer, budget=40, seed=7):
    calls = []

    def score_many(batch):
        calls.extend(batch)
        return scorer(batch)

    summary = run_method(
        method, initials(), score_many, seed=seed, budget=budget, pool_size=20,
        top_percent=0.5, target_pop_size=8, p_crossover=0.2, n_mut=1,
        fold_mode="stand-in", keep_events=True,
    )
    return summary, calls


class FairSearchTests(unittest.TestCase):
    def test_budget_counts_failures_and_successes(self):
        failed, failed_calls = run("crs", rejecting, budget=12)
        passed, passed_calls = run("ga", accepting, budget=12)
        self.assertEqual(failed["used"], 12)
        self.assertEqual(failed["failed"], 12)
        self.assertEqual(failed["gain"], 0.0)
        self.assertEqual(len(failed_calls), 12)
        self.assertEqual(passed["used"], passed["scored"])
        self.assertEqual(len(passed_calls), 12)
        self.assertGreater(passed["gain"], 0.0)

    def test_cumulative_search_stacks_edits_and_replays(self):
        first, _ = run("crs", accepting, budget=40, seed=11)
        second, _ = run("crs", accepting, budget=40, seed=11)
        self.assertGreaterEqual(first["best_depth"], 2)
        self.assertEqual(first["events"], second["events"])
        self.assertTrue(first["complete"])

    def test_methods_replay_independently(self):
        ga_a, _ = run("ga", accepting, budget=24, seed=19)
        ga_b, _ = run("ga", accepting, budget=24, seed=19)
        crs, _ = run("crs", accepting, budget=24, seed=19)
        self.assertEqual(ga_a["events"], ga_b["events"])
        self.assertNotEqual(ga_a["events"], crs["events"])
