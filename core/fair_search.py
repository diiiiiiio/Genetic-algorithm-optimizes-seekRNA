"""Budgeted GA and cumulative random search on one scoring configuration.

The shared budget counts the first full scoring call of each sequence, including
calls that return no score. Repeats are cache hits. Sequences that fail the
cheap sequence filter do not use the budget. A cumulative random search draws
parents uniformly from every accepted candidate, so later edits can stack.
It does not use score-weighted selection or crossover.
"""

import hashlib
import random
from copy import deepcopy

import numpy as np

from core.operators import crossover, mutate_in_targets
from core.scoring import cheap_filter, cheap_score

COORD_KEYS = [
    "a1", "a2", "A1", "A2", "b1", "b2", "B1", "B2",
    "a1m", "a2m", "A1m", "A2m", "b1m", "b2m", "B1m", "B2m",
]
MOTIF_PAIRS = (
    ("a1", "a2", "a1m", "a2m"),
    ("A1", "A2", "A1m", "A2m"),
    ("b1", "b2", "b1m", "b2m"),
    ("B1", "B2", "B1m", "B2m"),
)
EVENT_FIELDS = [
    "method", "seed", "request", "round", "anon_id", "edit_depth", "length",
    "fold_mode", "status", "budget_used", "score", "seconds", "motif_contained",
]


def anon_id(seq):
    return hashlib.sha256(seq.encode()).hexdigest()[:16]


def motifs_inside(coords):
    for start, end, motif_start, motif_end in MOTIF_PAIRS:
        if not (coords[start] <= coords[motif_start] < coords[motif_end] <= coords[end]):
            return False
    return True


def signature(coords, insertion):
    return tuple(int(coords[key]) for key in COORD_KEYS) + (int(insertion),)


def copy_individual(ind):
    return {
        "seq": ind["seq"],
        "coords": {key: int(ind["coords"][key]) for key in COORD_KEYS},
        "i": int(ind["i"]),
        "score": float(ind["score"]),
        "depth": int(ind["depth"]),
    }


def make_streams(seed, method):
    """Independent Python and NumPy streams for one method and seed."""
    tag = 1 if method == "ga" else 2
    np_rng = np.random.default_rng(np.random.SeedSequence([int(seed), tag]))
    py_rng = random.Random(f"{int(seed)}:{method}")
    return np_rng, py_rng


def _weights(population):
    weights = []
    for ind in population:
        weights.append(max(0.0, float(cheap_score(ind["seq"], ind["coords"]))) + 0.05)
    return weights


def _weighted_pick(population, weights, rng):
    total = sum(weights)
    if total <= 0.0:
        return rng.choice(population)
    point = rng.random() * total
    running = 0.0
    for ind, weight in zip(population, weights):
        running += weight
        if point <= running:
            return ind
    return population[-1]


def _propose(method, population, weights, np_rng, py_rng, p_crossover, n_mut):
    if method == "ga" and len(population) >= 2 and py_rng.random() < p_crossover:
        left = _weighted_pick(population, weights, py_rng)
        right = _weighted_pick(population, weights, py_rng)
        result = crossover({**left["coords"], "seq": left["seq"], "i": left["i"]},
                           {**right["coords"], "seq": right["seq"], "i": right["i"]}, rng=py_rng)
        if result is None:
            return None
        seq, coords, insertion = result
        depth = max(int(left["depth"]), int(right["depth"])) + 1
        return seq, {key: int(coords[key]) for key in COORD_KEYS}, int(insertion), depth
    parent = population[py_rng.randrange(len(population))] if method != "ga" else _weighted_pick(population, weights, py_rng)
    seq, coords = mutate_in_targets(parent["seq"], dict(parent["coords"]), n_events=n_mut, rng=np_rng)
    return seq, {key: int(coords[key]) for key in COORD_KEYS}, int(parent["i"]), int(parent["depth"]) + int(n_mut)


def _retain(method, population, accepted, target_pop_size):
    if method != "ga":
        known = {ind["seq"] for ind in population}
        extended = list(population)
        for ind in accepted:
            if ind["seq"] not in known:
                extended.append(ind)
                known.add(ind["seq"])
        return extended
    best = {ind["seq"]: ind for ind in population}
    for ind in accepted:
        current = best.get(ind["seq"])
        if current is None or ind["score"] > current["score"]:
            best[ind["seq"]] = ind
    ranked = sorted(best.values(), key=lambda ind: ind["score"], reverse=True)
    return ranked[:target_pop_size]


def run_method(method, initials, score_many, *, seed, budget, pool_size, top_percent,
               target_pop_size, p_crossover, n_mut, fold_mode, event_writer=None,
               keep_events=False, progress=None):
    """Run one method until ``budget`` distinct full-score calls are used."""
    if method not in ("ga", "crs"):
        raise ValueError("method must be 'ga' or 'crs'")
    if budget < 1:
        raise ValueError("budget must be positive")
    population = [copy_individual(ind) for ind in initials]
    if not population:
        raise ValueError("initial population is empty")
    np_rng, py_rng = make_streams(seed, method)
    start = max(population, key=lambda ind: ind["score"])
    best_score = float(start["score"])
    best_depth = int(start["depth"])
    best_length = len(start["seq"])
    best_id = anon_id(start["seq"])
    best_seq = start["seq"]
    best_coords = dict(start["coords"])
    seen = {ind["seq"]: signature(ind["coords"], ind["i"]) for ind in population}
    quota_full = max(1, int(pool_size * top_percent))
    used = scored = failed = cached = conflicts = prefiltered = proposal_count = 0
    motif_ok = motif_scored = 0
    depth_sum = 0
    request = 0
    rounds = 0
    empty_rounds = 0
    events = [] if keep_events else None
    stop = "budget"

    def emit(row):
        if event_writer is not None:
            event_writer(row)
        if events is not None:
            events.append(row)

    while used < budget:
        rounds += 1
        quota = min(budget - used, quota_full)
        weights = _weights(population) if method == "ga" else None
        proposals = []
        attempts = 0
        limit = pool_size * 10
        while len(proposals) < pool_size and attempts < limit:
            attempts += 1
            proposed = _propose(method, population, weights, np_rng, py_rng, p_crossover, n_mut)
            if proposed is None:
                continue
            seq, coords, insertion, depth = proposed
            proposal_count += 1
            if not cheap_filter(seq, coords):
                prefiltered += 1
                request += 1
                emit(_event(method, seed, request, rounds, seq, depth, fold_mode, "prefiltered", used, "", 0.0, ""))
                continue
            proposals.append({
                "seq": seq, "coords": coords, "i": insertion, "depth": depth,
                "cheap": float(cheap_score(seq, coords)),
            })
        proposals.sort(key=lambda item: item["cheap"], reverse=True)
        batch = []
        for item in proposals:
            if len(batch) >= quota:
                break
            previous = seen.get(item["seq"])
            current = signature(item["coords"], item["i"])
            if previous is None:
                seen[item["seq"]] = current
                batch.append(item)
                continue
            request += 1
            status = "cached" if previous == current else "conflict"
            if previous == current:
                cached += 1
            else:
                conflicts += 1
            emit(_event(method, seed, request, rounds, item["seq"], item["depth"], fold_mode, status, used, "", 0.0, ""))
        if not batch:
            empty_rounds += 1
            if empty_rounds >= 5:
                stop = "exhausted"
                break
            continue
        empty_rounds = 0
        results = score_many(batch)
        if len(results) != len(batch):
            raise RuntimeError("scorer returned a different number of results")
        accepted = []
        for item, returned in zip(batch, results):
            result, seconds = returned
            used += 1
            request += 1
            contained = motifs_inside(item["coords"])
            motif_scored += 1
            motif_ok += int(contained)
            depth_sum += int(item["depth"])
            if result is None or "Score" not in result:
                failed += 1
                emit(_event(method, seed, request, rounds, item["seq"], item["depth"], fold_mode, "failed", used, "", seconds, contained))
                continue
            score = float(result["Score"])
            scored += 1
            if score > best_score:
                best_score = score
                best_depth = int(item["depth"])
                best_length = len(item["seq"])
                best_id = anon_id(item["seq"])
                best_seq = item["seq"]
                best_coords = dict(item["coords"])
            accepted.append({
                "seq": item["seq"], "coords": item["coords"], "i": item["i"],
                "score": score, "depth": int(item["depth"]),
            })
            emit(_event(method, seed, request, rounds, item["seq"], item["depth"], fold_mode, "scored", used, score, seconds, contained))
        population = _retain(method, population, accepted, target_pop_size)
        if progress is not None:
            progress(method, seed, rounds, used, best_score)

    return {
        "method": method,
        "seed": int(seed),
        "budget": int(budget),
        "used": used,
        "scored": scored,
        "failed": failed,
        "cached": cached,
        "conflicts": conflicts,
        "prefiltered": prefiltered,
        "proposals": proposal_count,
        "rounds": rounds,
        "stop": stop,
        "complete": used == budget,
        "start_best": float(start["score"]),
        "best": best_score,
        "gain": best_score - float(start["score"]),
        "best_depth": best_depth,
        "best_length": best_length,
        "best_id": best_id,
        "best_seq": best_seq,
        "best_coords": best_coords,
        "motif_contained": motif_ok,
        "motif_checked": motif_scored,
        "mean_scored_depth": (depth_sum / motif_scored) if motif_scored else 0.0,
        "events": events,
    }


def score_record(args):
    """Score one candidate. Used by worker processes."""
    import time
    from core.scoring import Phy_mark

    seq, ref_a, ref_b, insertion, coords, fold_mode = args
    started = time.perf_counter()
    result = Phy_mark(seq, ref_a, ref_b, insertion, coords, fold_mode=fold_mode)
    return result, time.perf_counter() - started


def _event(method, seed, request, round_index, seq, depth, fold_mode, status, used, score, seconds, contained):
    return {
        "method": method,
        "seed": int(seed),
        "request": int(request),
        "round": int(round_index),
        "anon_id": anon_id(seq),
        "edit_depth": int(depth),
        "length": len(seq),
        "fold_mode": fold_mode,
        "status": status,
        "budget_used": int(used),
        "score": score,
        "seconds": seconds,
        "motif_contained": contained,
    }
