"""Summarize length and stored scores of the three archive tables. No row-level output."""

import hashlib
import io
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
FILES = {'initial': 'initial_candidates.csv', 'GA': 'ga_results.csv', 'RS': 'rs_results.csv'}


def read_archive(filename):
    data = (ROOT / filename).read_bytes()
    table = pd.read_csv(io.BytesIO(data), usecols=['seq', 'Score', 'GC'])
    if table.empty or table['seq'].isna().any():
        raise ValueError(f'{filename} has missing sequences or no rows')
    table['length'] = table['seq'].str.len()
    if (table['length'] == 0).any() or not np.isfinite(table[['Score', 'GC']]).all().all():
        raise ValueError(f'{filename} has an empty record or a non-finite field')
    # This stored field identifies the archival inconsistencies noted in the prior audit.
    gc_from_text = table['seq'].str.count('[GC]') / table['length']
    table['consistent_gc'] = np.isclose(gc_from_text, table['GC'], rtol=0, atol=1e-8)
    return table, hashlib.sha256(data).hexdigest()


def summarize(table):
    length = table['length']
    score = table['Score']
    group_rows = []
    for nt, group in table.groupby('length', sort=True):
        points = group['Score']
        group_rows.append(dict(length=int(nt), n=len(group), mean=float(points.mean()),
                               q1=float(points.quantile(.25)), median=float(points.median()),
                               q3=float(points.quantile(.75))))
    rho = None
    if length.nunique() > 1 and score.nunique() > 1:
        rho = float(length.rank(method='average').corr(score.rank(method='average')))
    return dict(n=len(table), unique=int(table['seq'].nunique()), mean_score=float(score.mean()),
                mean_length=float(length.mean()), min_length=int(length.min()),
                max_length=int(length.max()), rho=rho, lengths=group_rows)


def compare(ga, rs):
    ga_by_length = {r['length']: r for r in ga['lengths']}
    rs_by_length = {r['length']: r for r in rs['lengths']}
    common = []
    for nt in sorted(ga_by_length.keys() & rs_by_length.keys()):
        a, b = ga_by_length[nt], rs_by_length[nt]
        common.append(dict(length=nt, ga_n=a['n'], rs_n=b['n'], ga_mean=a['mean'],
                           rs_mean=b['mean'], difference=a['mean'] - b['mean'],
                           weight=min(a['n'], b['n'])))

    def weighted_difference(rows):
        weight = sum(r['weight'] for r in rows)
        value = math.fsum(r['weight'] * r['difference'] for r in rows) / weight if weight else None
        return dict(strata=len(rows), total_weight=weight, difference=value)

    return dict(unadjusted=ga['mean_score'] - rs['mean_score'], common=common,
                weighted=weighted_difference(common),
                min_five=weighted_difference([r for r in common if r['weight'] >= 5]))


def main():
    provenance, full, sensitivity = {}, {}, {}
    for label, filename in FILES.items():
        table, digest = read_archive(filename)
        full[label] = summarize(table)
        valid = table.loc[table['consistent_gc']]
        duplicate = valid['seq'].duplicated(keep=False)
        provenance[filename] = dict(sha256=digest, rows=len(table),
                                    inconsistent_metadata_rows=list(map(int, table.index[~table['consistent_gc']])),
                                    repeated_rows_after_metadata_filter=int(duplicate.sum()))
        if label != 'initial':
            sensitivity[label] = summarize(valid.loc[~duplicate])
    output = dict(inputs=provenance, tables=full, comparison=compare(full['GA'], full['RS']),
                  sensitivity=dict(tables=sensitivity, comparison=compare(sensitivity['GA'], sensitivity['RS'])))
    print(json.dumps(output, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
