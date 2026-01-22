import random
import bisect
import pandas as pd

def roulette_sampler(items, weights):
    """Return a sampler() that draws one item with probability ∝ weight.
    Fallback to uniform if all weights ≤ 0 or list empty.
    """
    items = list(items)
    weights = [max(0.0, float(w)) for w in weights]
    totals = []
    running = 0.0
    for w in weights:
        running += w
        totals.append(running)
    if running <= 0.0 or not items:
        def pick():
            return random.choice(items)
        return pick
    def pick():
        r = random.random() * running
        idx = bisect.bisect(totals, r)
        if idx >= len(items):
            idx = len(items) - 1
        return items[idx]
    return pick


def weighted_sample_df(df, weight_col: str, k: int = 1):
    """Sample k rows from DataFrame proportionally to df[weight_col].
    Returns a list of pandas Series (rows). Fallback to uniform if needed.
    """
    if df is None or len(df) == 0:
        return []
    items = [row for _, row in df.iterrows()]
    weights = [float(row.get(weight_col, 0.0)) for _, row in df.iterrows()]
    pick = roulette_sampler(items, weights)
    return [pick() for _ in range(k)]

