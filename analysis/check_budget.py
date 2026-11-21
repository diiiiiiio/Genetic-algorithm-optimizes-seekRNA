"""Exercise accounting and seed replay with numbered tickets and artificial scores."""

import hashlib
import json
import random


SEEDS = (104729, 130363, 155921)
BUDGET = 32
START_ID, START_SCORE = -1, 25


def tickets(seed, label):
    items = [START_ID, *range(1, 61), *range(1, 61)]
    random.Random(f'{seed}/{label}').shuffle(items)
    # An initial cache hit, a prefilter rejection, a failure and a cached failure.
    return [START_ID, 0, 7, 7, *items]


def count_events(items, budget):
    seen = {START_ID}
    best = START_SCORE
    calls = scored = failed = cached = filtered = 0
    trace = []
    for item in items:
        if calls == budget:
            break
        if item % 13 == 0:
            state = 'prefiltered'
            filtered += 1
        elif item in seen:
            state = 'cached'
            cached += 1
        else:
            seen.add(item)
            calls += 1
            if item % 7 == 0:
                state = 'failed'
                failed += 1
            else:
                state = 'scored'
                best = max(best, (11 * item) % 101)
                scored += 1
        trace.append((item, state, calls, best))
    assert calls == scored + failed
    assert len(trace) == calls + cached + filtered
    assert all(a[3] <= b[3] for a, b in zip(trace, trace[1:]))
    digest = hashlib.sha256(json.dumps(trace, separators=(',', ':')).encode()).hexdigest()
    return dict(requests=len(trace), used=calls, scored=scored, failed=failed, cached=cached,
                filtered=filtered, complete=calls == budget, gain=best - START_SCORE, trace_sha256=digest)


def main():
    example = count_events([START_ID, 0, 7, 7, 5, 8], 2)
    assert (example['requests'], example['used'], example['scored'], example['failed'],
            example['cached'], example['filtered'], example['gain']) == (5, 2, 1, 1, 2, 1, 30)
    assert not count_events([START_ID, 0, 7, 7], 2)['complete']
    runs = []
    for seed in SEEDS:
        streams = {key: count_events(tickets(seed, key), BUDGET) for key in ('A', 'B')}
        reversed_order = {key: count_events(tickets(seed, key), BUDGET) for key in ('B', 'A')}
        assert streams == reversed_order
        assert all(row['complete'] for row in streams.values())
        runs.append(dict(seed=seed, streams=streams))
    assert len({row['trace_sha256'] for r in runs for row in r['streams'].values()}) == 2 * len(SEEDS)
    print(json.dumps(dict(dataset='integer tickets only', budget=BUDGET, runs=runs), indent=2))


if __name__ == '__main__':
    main()
