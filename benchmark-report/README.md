# Numerical results of the paired seekRNA search benchmark

This directory archives the statistical analysis of 20 paired runs of a genetic algorithm (GA) and cumulative random search (CRS). Each method had 7,500 first scoring requests per seed. The 11 shared starting records were scored separately. A request that returned no score consumed one unit of budget.

Run `python summarize.py` in an environment with NumPy and pandas to regenerate the published numerical summary from the CSV files in this directory. The manuscript uses four-decimal rounding. The full-precision outputs are in `reported_results.json`.

| File | Contents |
| --- | --- |
| `protocol.json` | Formal seed list, budget, start score, endpoint, bootstrap settings and software versions. |
| `seed_outcomes.csv` | One row per seed pair. Gains are terminal maxima minus the common start. Contained gains use the shared starts and records passing the annotation check in the original search history. Winner flags refer to the single stored terminal winner. |
| `budget_checkpoints.csv` | Quartiles and medians of cumulative best gain over the 20 seeds, at intervals of 250 requests. |
| `winner_components.csv` | Means and sample standard deviations across the 20 stored winners per method, paired mean differences, and additive contributions to the score difference. The existing scoring implementation was applied to the archived winners using their stored annotations. All 40 totals matched the recorded scores. The GC fraction is descriptive. |
| `length_strata.csv` | Shared seed-by-length groups, successful counts, observed maxima, and expected maxima under equally sized samples drawn without replacement. Within a seed, the analysis weights each shared length by the smaller count. Seeds receive equal weight. |
| `pooled_length_correlation.csv` | Length–score Spearman correlations over the successful evaluations of each method. All tied observations receive their average rank. |
| `audit_events.py` | Statistical extraction code for completed local event logs. It generated the seed, checkpoint, length-stratum and pooled-correlation tables. |
| `source_fingerprints.json` | SHA-256 fingerprints of inputs and implementation files from the source workspace used for the benchmark, together with the public prototype commit. |

The public source at the repository root is the original prototype. The formal search implementation and full event archive were held in the authors' analysis workspace. The code in this directory operates on numerical outputs. Bootstrap intervals resample complete seed pairs and quantify variation across runs for this input pair and the chosen score.
