# SeekRNA GA Project

A genetic algorithm framework for optimizing RNA sequences. It evolves candidates through mutation and crossover, evaluating fitness via two-tier scoring: fast sequence filtering and comprehensive structural analysis using RNA folding tools to maximize loop accessibility and stem conservation.

---

## Table of Contents

- [Project Structure](#project-structure)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Configuration](#configuration)
- [How to Run](#how-to-run)
- [Algorithm Details](#algorithm-details)
  - [1. Mutation Logic](#1-mutation-logic)
  - [2. Crossover Logic](#2-crossover-logic)
  - [3. Coordinate Logic](#3-coordinate-logic)
  - [4. Scoring Logic](#4-scoring-logic)

---

## Project Structure

```
SeekRNA_GA_Project/
├── main.py                 # Entry point for the GA
├── config.json             # Configuration for sequences, target sites, and GA parameters
├── config_loader.py        # Utility to load configurations
├── requirements.txt        # Python dependencies
├── initial_candidates_ISEc11_ISPa11.csv  # Initial population data
├── core/
│   ├── ga_engine.py        # Main GA loop and parallel processing
│   ├── operators.py        # Mutation and crossover logic
│   ├── scoring.py          # Tier 1 and Tier 2 scoring functions
│   └── sampling.py         # Weighted selection methods
└── utils/
    ├── sequence.py         # Basic RNA sequence utilities
    ├── coords.py           # Coordinate mapping and update logic
    └── folding.py          # External folding tool wrappers
```

---

## Prerequisites

### Python Dependencies

- **Python 3.7+**
- **pandas** - Data manipulation and analysis
- **numpy** - Numerical computing

### External RNA Folding Tools

The following tools must be installed and available in your system `PATH`:

| Tool | Description | Installation |
|------|-------------|--------------|
| **RNAfold** | ViennaRNA Package - thermodynamic structure prediction | [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) |
| **mxfold2** | Deep learning-based RNA folding | `pip install mxfold2` or [GitHub](https://github.com/mxfold/mxfold2) |
| **centroid_fold** | CentroidFold - ensemble-based prediction | [CentroidFold](https://github.com/satoken/centroid-rna-package) |

To verify installation, run:
```bash
RNAfold --version
mxfold2 --help
centroid_fold --help
```

---

## Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/diiiiiiio/Genetic-algorithm-optimizes-seekRNA.git
   cd Genetic-algorithm-optimizes-seekRNA
   ```

2. **Install Python dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

3. **Verify external tools are installed** (see [Prerequisites](#prerequisites))

4. **Place your initial candidates CSV file** in the project root:
   - The file should be named `initial_candidates_ISEc11_ISPa11.csv`
   - It must contain columns: `seq`, `i`, `a1`, `a2`, `A1`, `A2`, `b1`, `b2`, `B1`, `B2`, and optionally motif coordinates

---

## Configuration

Edit `config.json` to customize the GA parameters:

```json
{
    "sequences": {
        "A": "ACAGGUAAGACCGUGACUCACUAAACCUGAACAGUAUUUUGGGCUUGAAGUCCGCCGUGAAAAUAAGGGGUGAGUCGGCGAA",
        "B": "AACCUGGUUUAUACGUGGGCUCCCUGCUGCAGUGAAGCAAAGCCGUUAGGGCGAUCAGGU"
    },
    "target_sites": {
        "j1": 30, "j2": 40,
        "J1": 56, "J2": 65,
        "k1": 8,  "k2": 15,
        "K1": 46, "K2": 54
    },
    "motif_lengths": {
        "A_MOTIF_LEN": 6,
        "B_MOTIF_LEN": 3
    },
    "ga_parameters": {
        "n_mut_positions": 1,
        "p_crossover": 0.2,
        "rounds": 30,
        "pool_size": 1000,
        "top_percent": 0.25,
        "elite_size": 100,
        "target_pop_size": 200,
        "fold_mode": "full",
        "n_workers": 7
    }
}
```

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `rounds` | 30 | Number of GA generations |
| `pool_size` | 1000 | Candidates generated per round |
| `top_percent` | 0.25 | Fraction evaluated with expensive scoring |
| `target_pop_size` | 200 | Population size after selection |
| `p_crossover` | 0.2 | Probability of crossover vs mutation |
| `fold_mode` | "full" | Folding mode: "cheap", "medium", or "full" |
| `n_workers` | 7 | Parallel workers for Phy_mark evaluation |

---

## How to Run

```bash
python main.py
```

### Output

- Results are saved to `ga_results.csv`
- The top 5 candidates are printed with motif verification
- Each round displays progress: candidates generated, Phy_mark evaluations, best score

### Example Output

```
Folding reference sequences A and B...
=== Round 1/30 ===
  Generated 1000 candidates (cheap filter pass)
  Evaluating top 250 candidates with Phy_mark (mode=full)...
  187 candidates passed Phy_mark
  Population size: 200, Best score: 2.8451

...

GA complete. Results saved to ga_results.csv

==================================================
Motif Verification for Top 5 Candidates:
==================================================
```

---

## Algorithm Details

This section provides a detailed explanation of the **mutation**, **crossover**, **coordinate updates**, and **scoring** logic.

---

### 1. Mutation Logic

All mutations happen on a chimeric RNA **C = A[:i] + B + A[i:]**, and are confined to four recognition regions:

- **A top:** `[a1, a2)`
- **A bottom:** `[A1, A2)`
- **B top:** `[b1, b2)`
- **B bottom:** `[B1, B2)`

With shorter **motif sub-regions** inside each:

- A top motif: `[a1m, a2m)` (length 6)
- A bottom motif: `[A1m, A2m)` (length 6)
- B top motif: `[b1m, b2m)` (length 3)
- B bottom motif: `[B1m, B2m)` (length 3)

The main mutation driver is:

```python
mutate_in_targets(seq, coords, n_events=1, rng=None, p_pair=0.7,
                  single_probs=(0.6,0.2,0.2), pair_probs=(0.6,0.2,0.2))
```

#### 1.1 High-level mutation flow: `mutate_in_targets`

For each mutation event (there are `n_events`, usually 1):

- **Check coordinate sanity**  
  If `sanity_check_coords(seq_cur, coords_cur)` fails (any of the four recognition segments goes out of bounds), it stops mutating further.

- **Decide mutation type (paired vs single-strand)**
  - With probability `p_pair` (default **0.7**): do a **paired mutation** in the motif of A or B.
  - With probability `1 - p_pair` (default **0.3**): do a **single-strand mutation** somewhere in one recognition region but **outside** its motif.

- **If paired mutation:** use `mutate_paired_segments`
  - First randomly select **A_pair** or **B_pair** with equal probability.
  - For A_pair it calls:
    ```python
    mutate_paired_segments(
        seq_cur, coords_cur,
        seg_top_pair=("a1","a2"), seg_bot_pair=("A1","A2"),
        motif_top_pair=("a1m","a2m"), motif_bot_pair=("A1m","A2m"),
        rng=rng, p_sub=p_sub_pair, p_ins=p_ins_pair, p_del=p_del_pair
    )
    ```
  - For B_pair similarly but with `b1,b2,B1,B2` and motif keys `b1m,b2m,B1m,B2m`.

- **If single-strand mutation:** use `mutate_single_segment`
  - First choose **one of four** regions uniformly: "Atop", "Abelow", "Btop", "Bbelow".
  - Map that choice to a segment and its motif:
    - Atop → seg `("a1","a2")`, motif `("a1m","a2m")`
    - Abelow → seg `("A1","A2")`, motif `("A1m","A2m")`
    - Btop → seg `("b1","b2")`, motif `("b1m","b2m")`
    - Bbelow → seg `("B1","B2")`, motif `("B1m","B2m")`
  - Then call:
    ```python
    mutate_single_segment(seq_cur, coords_cur, seg, motif,
                          rng=rng, p_sub=p_sub_single,
                          p_ins=p_ins_single, p_del=p_del_single)
    ```

- Each event returns updated `(seq_cur, coords_cur)` and those are fed into the next event.

#### 1.2 Single-segment mutation: `mutate_single_segment`

**Goal:** Mutate **one recognition segment** (top OR bottom of A or B) but never directly touch the motif itself.

**Inputs:**
- `seg_key_pair`: e.g. `("a1", "a2")`
- `motif_key_pair`: e.g. `("a1m", "a2m")`
- `p_sub, p_ins, p_del`: probabilities for substitution, insertion, deletion.

**Steps:**

1. **Segment and motif ranges**
   ```python
   s_key, e_key = seg_key_pair
   s, e = coords[s_key], coords[e_key]
   length = e - s
   if length <= 1: return seq, coords
   
   m_start = coords[motif_key_pair[0]]
   m_end = coords[motif_key_pair[1]]
   ```

2. **Candidate positions (exclude motif)**  
   Only positions in `[s, e)` **outside** motif `[m_start, m_end)` are allowed:
   ```python
   candidates = [p for p in range(s, e) if not (m_start <= p < m_end)]
   if not candidates:
       return seq, coords
   pos = int(rng.choice(candidates))
   ```

3. **Choose edit type**  
   Normalize `p_sub, p_ins, p_del` and draw a random number `r`:
   ```python
   total_p = p_sub + p_ins + p_del
   p_sub_n = p_sub / total_p
   p_ins_n = p_ins / total_p
   # deletion prob = remaining
   
   if r < p_sub_n: edit_type = "sub"
   elif r < p_sub_n + p_ins_n: edit_type = "ins"
   else: edit_type = "del"
   ```

4. **Apply edit via `apply_edit`**
   - **Substitution** (`edit_type == "sub"`)
     - Pick a new base from "ACGU" different from the old base at pos.
     - `new_seq, cut_s, cut_e, delta = apply_edit(seq, pos, "sub", nts=new_nt)`
     - `# delta = 0, cut_s = pos, cut_e = pos+1`
   - **Insertion** (`edit_type == "ins"`)
     - Choose insertion length 1-2:
     - `ins_len = rng.integers(1, 3)`
     - `nts = ''.join(rng.choice(list("ACGU"), size=ins_len))`
     - `new_seq, cut_s, cut_e, delta = apply_edit(seq, pos, "ins", nts=nts)`
     - `# delta = +ins_len`
   - **Deletion** (`edit_type == "del"`)
     - Fixed `del_len = 1` (single base deletion), clamped to stay inside the segment.
     - `new_seq, cut_s, cut_e, delta = apply_edit(seq, pos, "del", del_len=1)`
     - `# delta = -1`

5. **Update all coordinates**  
   After any edit, recalculate **all** segment & motif coordinates:
   ```python
   new_coords = update_coords_after_edit(coords, cut_s, cut_e, delta, len(new_seq))
   ```

6. **Return** `new_seq, new_coords`

So `mutate_single_segment` is "local, single-strand edit in one recognition region, outside the motif, with coordinates globally updated afterwards."

#### 1.3 Paired motif mutation: `mutate_paired_segments`

**Goal:** Mutate **both top and bottom motifs together** so that they remain reverse-complementary.

**Inputs:**
- `seg_top_pair`: e.g. `("a1","a2")`
- `seg_bot_pair`: e.g. `("A1","A2")`
- `motif_top_pair`: e.g. `("a1m","a2m")`
- `motif_bot_pair`: e.g. `("A1m","A2m")`
- `p_sub, p_ins, p_del` for the **paired** edits.

**Steps:**

1. **Motif boundaries** (top & bottom):
   ```python
   t_m1, t_m2 = coords[motif_top_pair[0]], coords[motif_top_pair[1]]
   b_m1, b_m2 = coords[motif_bot_pair[0]], coords[motif_bot_pair[1]]
   L_motif = t_m2 - t_m1
   if L_motif <= 0: return seq, coords
   ```

2. **Pick edit type (sub / ins / del)**  
   Same normalization approach as above, but used for the **pair**.

3. **Pick offset within motif**
   - If insertion and motif length > 2, choose offset in `[1, L_motif-1)` (avoid ends).
   - Otherwise offset in `[0, L_motif)`.
   
   Top-strand position to edit: `t_pos = t_m1 + offset`

4. **Case 1: Paired substitution**
   - Choose `new_nt ≠ old_nt` at `t_pos`.
   - Apply top substitution:
     ```python
     new_seq, cut_s, cut_e, delta = apply_edit(seq, t_pos, "sub", nts=new_nt)
     new_coords = update_coords_after_edit(coords, cut_s, cut_e, delta, len(new_seq))
     # delta = 0 for substitution
     ```
   - Recompute bottom motif end in new coordinates:
     ```python
     b_m2_n = new_coords[motif_bot_pair[1]]
     b_pos = b_m2_n - 1 - offset
     ```
     This ensures the bottom position mirrors the same offset from the **end** of the bottom motif.
   - Complement base for bottom:
     ```python
     comp_nt = rev_comp(new_nt)  # single base rev-comp
     new_seq2, cut_s2, cut_e2, delta2 = apply_edit(new_seq, b_pos, "sub", nts=comp_nt)
     new_coords2 = update_coords_after_edit(new_coords, cut_s2, cut_e2, delta2, len(new_seq2))
     ```
   - Return `(new_seq2, new_coords2)`.

5. **Case 2: Paired insertion**
   - Insert 1-2 nucleotides on top motif at `t_pos`:
     ```python
     ins_len = rng.integers(1, 3)
     nts = ''.join(rng.choice(list("ACGU"), size=ins_len))
     new_seq, cut_s, cut_e, delta = apply_edit(seq, t_pos, "ins", nts=nts)
     new_coords = update_coords_after_edit(coords, cut_s, cut_e, delta, len(new_seq))
     ```
   - Recompute bottom motif end after insertion, then the bottom insertion point:
     ```python
     b_m2_n = new_coords[motif_bot_pair[1]]
     b_pos = b_m2_n - offset
     ```
   - Insert reverse-complement on bottom:
     ```python
     nts_bot = rev_comp(nts)
     new_seq2, cut_s2, cut_e2, delta2 = apply_edit(new_seq, b_pos, "ins", nts=nts_bot)
     new_coords2 = update_coords_after_edit(new_coords, cut_s2, cut_e2, delta2, len(new_seq2))
     ```

6. **Case 3: Paired deletion**
   - Always delete **1 base** on top at `t_pos`:
     ```python
     del_len = 1
     new_seq, cut_s, cut_e, delta = apply_edit(seq, t_pos, "del", del_len=del_len)
     new_coords = update_coords_after_edit(coords, cut_s, cut_e, delta, len(new_seq))
     ```
   - After top deletion, recompute bottom motif end and mirrored offset:
     ```python
     b_m2_n = new_coords[motif_bot_pair[1]]
     b_pos = b_m2_n - 1 - offset
     ```
   - Delete 1 base at `b_pos`:
     ```python
     new_seq2, cut_s2, cut_e2, delta2 = apply_edit(new_seq, b_pos, "del", del_len=del_len)
     new_coords2 = update_coords_after_edit(new_coords, cut_s2, cut_e2, delta2, len(new_seq2))
     ```

So paired mutations always keep top and bottom motif **length-matched** and **reverse-complementary**, with coordinates updated after each edit.

---

### 2. Crossover Logic

Crossover is very conservative: it only happens when two parents have **identical coordinates and insertion index**.

```python
def crossover(p1, p2):
    """Simple one-point crossover (only when parent coordinates are completely identical)"""
```

**Steps:**

1. **Check basic compatibility**
   ```python
   s1, s2 = p1["seq"], p2["seq"]
   if len(s1) != len(s2): return None
   i1, i2 = int(p1["i"]), int(p2["i"])
   if i1 != i2: return None
   ```

2. **Check all coordinate keys are identical**  
   The key list includes:
   ```python
   keys = ["a1","a2","A1","A2","b1","b2","B1","B2",
           "a1m","a2m","A1m","A2m","b1m","b2m","B1m","B2m"]
   ```
   If any coordinate differs, crossover is aborted (returns `None`).

3. **One-point crossover**
   - If `n <= 1`, return `None` (sequence too short for crossover).
   - Choose a cut point between 1 and n-1:
     ```python
     cut = random.randint(1, n-1)
     child_seq = s1[:cut] + s2[cut:]
     ```
   - **Coordinates are copied unchanged from p1** (since identical in both parents):
     ```python
     child_coords = {k: int(p1[k]) for k in keys}
     ```
   - Return `(child_seq, child_coords, i1)`.

**In the GA loop:**

- With probability `p_crossover` (~0.2), it attempts crossover using two parents drawn by **roulette sampling** on `cheap_w` (cheap_score).
- Otherwise it mutates a single parent using `mutate_in_targets`.
- After generating each child, it runs `cheap_filter` and `sanity_check_coords` before keeping it.

---

### 3. Coordinate Logic

There are three related pieces:

- **Mapping A/B intervals into the chimeric sequence C**
- **Tracking motif coordinates**
- **Updating all coordinates after every edit**

#### 3.1 Mapping reference intervals into C: `map_interval`

C is constructed as: **C = A[:i] + B + A[i:]**

`map_interval(start, end, i_insert, L_insert)` maps an interval on A (`[start,end)`) into C:

- If insertion index `i_insert <= start`: shift right by `L_insert`.
- If `i_insert >= end`: coordinates unchanged.
- If `start < i_insert < end`: the insertion splits the region → return `(None, None)`.

Used in `Phy_mark` when no explicit coords are provided:
```python
a1, a2 = map_interval(j1, j2, i, L)  # A top
A1, A2 = map_interval(J1, J2, i, L)  # A bottom
b1, b2 = i + k1, i + k2              # B top
B1, B2 = i + K1, i + K2              # B bottom
```

Motif coordinates when `coords` is None:

- A motifs (length 6):
  ```python
  a1m, a2m = map_interval(j2 - A_MOTIF_LEN, j2, i, L)
  A1m, A2m = map_interval(J2 - A_MOTIF_LEN, J2, i, L)
  ```
- B motifs (length 3):
  ```python
  b1m, b2m = i + (k2 - B_MOTIF_LEN), i + k2
  B1m, B2m = i + (K2 - B_MOTIF_LEN), i + K2
  ```

If any of these resolve to `None` or go out of bounds, the candidate is rejected.

#### 3.2 Motif column creation in the DataFrame

When initial candidates are loaded from CSV, they might lack motif columns. `add_motif_columns_if_missing` fills them:

```python
df["a1m"] = df["a2"] - A_MOTIF_LEN
df["a2m"] = df["a2"]
df["A1m"] = df["A2"] - A_MOTIF_LEN
df["A2m"] = df["A2"]
df["b1m"] = df["b2"] - B_MOTIF_LEN
df["b2m"] = df["b2"]
df["B1m"] = df["B2"] - B_MOTIF_LEN
df["B2m"] = df["B2"]
```

So each motif is a **fixed-length block at the 3′ end** of its recognition segment.

#### 3.3 Coordinate updates after edits

Every sub/ins/del returns:
```python
new_seq, cut_start, cut_end, delta
```

Where:
- `[cut_start, cut_end)` is the edited window on the _old_ sequence.
- `delta = len(new_seq) - len(old_seq)` (positive for insertions, negative for deletions, 0 for substitutions).

Then all intervals are updated via:

**(a) `update_interval_after_edit(start, end, cut_start, cut_end, delta)`**

This updates a single `[start,end)` interval:

- If the interval is **entirely to the right** of the edit (`start > cut_end`): shift both start and end by delta.
- If **entirely to the left** (`end <= cut_start`): unchanged.
- If overlapping:
  - If `delta > 0` (insertion): the code extends end by delta (inserted bases inside the interval).
  - If `delta < 0` (deletion): compute the number of bases removed from the intersection and shrink the interval accordingly, possibly shifting it left if deletion started before start.

It then returns new `(start, end)`.

**(b) `update_coords_after_edit(coords, cut_start, cut_end, delta, seqlen_after)`**

This applies `update_interval_after_edit` to **all** segments:

```python
keys = ["a1","a2","A1","A2","b1","b2","B1","B2",
        "a1m","a2m","A1m","A2m","b1m","b2m","B1m","B2m"]

new_coords = dict(coords)
for s_key, e_key in zip(keys[0::2], keys[1::2]):
    s = new_coords[s_key]
    e = new_coords[e_key]
    s2, e2 = update_interval_after_edit(s, e, cut_start, cut_end, delta)
    # Clamp to valid range and ensure e2 > s2
    s2 = max(0, min(s2, seqlen_after))
    e2 = max(s2+1, min(e2, seqlen_after))
    new_coords[s_key], new_coords[e_key] = s2, e2
```

So after every single edit, all eight intervals (4 segments + 4 motifs) are shifted/shortened/extended appropriately.

**(c) Sanity checks**

- `sanity_check_coords(seq, coords)` enforces for each of the four main recognition segments:
  - `0 <= start < end <= len(seq)`
- `cheap_filter` (when coords provided) additionally enforces validity of the four main segments.

These guards prevent mutation operations from drifting coordinates into nonsense ranges.

---

### 4. Scoring Logic

There are two layers:

- A **cheap sequence-only score** used for filtering and parent weighting.
- A **full structural score** `Phy_mark` used to compute the final Score for the GA.

#### 4.1 Cheap filter & cheap score

**`cheap_filter(C, coords=None)`**

Fast "hard reject" test:

- **Homopolymer run**  
  Reject if the longest run of identical bases (`max_same_length(C)`) is > 7.
- **GC range**  
  Compute `gc = gc_count(C)` (fraction of G or C).  
  Reject if `gc < 0.25` or `gc > 0.75`.
- **Coordinate sanity (if coords provided)**
  - Check all required keys exist: `a1,a2,A1,A2,b1,b2,B1,B2`.
  - For each of the four segments make sure `0 <= s < e <= len(C)`.

If any rule fails → `False` (reject), else `True`.

**`cheap_score(C, coords)`**

Fast approximate fitness, _ignoring structures_, just sequence:

```python
gc = gc_count(C)
max_run = max_same_length(C)
score = 1.0

# GC penalty
if gc < 0.30 or gc > 0.70:
    score -= 0.3

# Homopolymer penalty
if max_run > 4:
    score -= 0.2 * (max_run - 4)

# Bonus for sweet-spot GC
if 0.40 <= gc <= 0.60:
    score += 0.1
```

This `cheap_score` is used for:
- Sorting `candidate_pool` during cheap stage
- Computing `cheap_w` (weights) for **roulette sampling** of parents.

#### 4.2 Full structural scoring: `Phy_mark`

```python
Phy_mark(C, struct_ref_A, struct_ref_B, i, coords=None, fold_mode="full")
```

Returns a dict with `Score` and many sub-scores, or `None` if the candidate fails.

**Step 1 - Determine coordinates**

- If `coords` is `None`, it reconstructs `a1,a2,A1,A2,b1,b2,B1,B2` and motif coordinates using `map_interval` and the insertion index `i` (as in section 3.1).
- If coords are provided, it just reads all eight intervals from there and checks all are within `[0, len(C)]`.

If any interval is invalid, it returns `None`.

**Step 2 - Hard filter**

Calls `cheap_filter(C, coords_dict)` with full coordinate dict (including motifs). If that fails → `None`.

**Step 3 - Folding**

It calls `Fold(C, mode=fold_mode)`, which returns:

- `stru`: a list/ensemble of dot-bracket strings. The number depends on `fold_mode`:
  - `"cheap"`: 1 structure (RNAfold only)
  - `"medium"`: 2 structures (RNAfold + mxfold2)
  - `"full"`: 3 structures (RNAfold + centroid_fold + mxfold2)
- `E`: list of normalized energy scores (one per structure)

If folding fails or stru is empty: return `None`.

If fewer than 3 structures are returned, the code pads by duplicating the first structure until there are 3.

**Step 4 - Loop accessibility**

For each of the four recognition segments:

```python
Btop = dot_ratio_ensemble(stru, b1, b2)   # B top
Bbot = dot_ratio_ensemble(stru, B1, B2)   # B bottom
Atop = dot_ratio_ensemble(stru, a1, a2)   # A top
Abot = dot_ratio_ensemble(stru, A1, A2)   # A bottom
```

`dot_ratio_ensemble` computes the **fraction of dots '.'** (unpaired) in `[a,b)` averaged over all three structures. This measures accessibility of those sites.

**Step 5 - Energy term**

```python
Emean = mean(E) if E else 0.0
energy_score = (Emean + 1.0) / 2.0
```

So more favorable energies shift the score. It rescales into a roughly 0-1 range.

**Step 6 - Sequence penalty (again)**

Similar but not identical to `cheap_score`:

```python
gc = gc_count(C)
max_run = max_same_length(C)
penalty = 0.0

if gc < 0.30 or gc > 0.70:
    penalty += 0.2

if max_run > 4:
    penalty += 0.5
```

This penalty is subtracted from the final score.

**Step 7 - Stem similarity vs references**

At the top of the GA, reference foldings are precomputed:

```python
struct_A, e_A = Fold(A)
struct_A1, struct_A2, struct_A3 = struct_A[0], struct_A[1], struct_A[2]

struct_B, e_B = Fold(B)
struct_B1, struct_B2, struct_B3 = struct_B[0], struct_B[1], struct_B[2]
```

For a candidate C:

- Focus on the **stem regions between top and bottom recognition sites**:
  - A stem region on C: from `a2` to `A1`.
    ```python
    A_key_start = max(0, a2)
    A_key_end = min(n, A1)
    A_key_len = max(0, A_key_end - A_key_start)
    ```
    Reference stem on A: between `j2` and `J1`.
    ```python
    A_ref_len = min(A_key_len, J1 - j2)
    ```
  - B stem region on C: from `b2` to `B1`.
    ```python
    B_key_start = max(0, b2)
    B_key_end = min(n, B1)
    B_key_len = max(0, B_key_end - B_key_start)
    B_ref_len = min(B_key_len, K1 - k2)
    ```

- For each of the three structures (1,2,3) they take the **candidate substring** from the C structure and match it to the **corresponding substring** from the A or B reference structure:
  - A region:
    ```python
    A_cand_struct1 = stru_C1[A_key_start : A_key_start + A_ref_len]
    A_ref_struct1 = struct_A1[j2 : j2 + A_ref_len]
    A_stem_score1 = pair_distance_score(A_ref_struct1, A_cand_struct1)
    # similarly A_stem_score2, A_stem_score3 from struct_A2, struct_A3
    ```
  - B region:
    ```python
    B_cand_struct1 = stru_C1[B_key_start : B_key_start + B_ref_len]
    B_ref_struct1 = struct_B1[k2 : k2 + B_ref_len]
    B_stem_score1 = pair_distance_score(B_ref_struct1, B_cand_struct1)
    # similarly B_stem_score2, B_stem_score3
    ```

`pair_distance_score` = fraction of reference base pairs preserved in the candidate:
- Parse dot-bracket into base-pair sets.
- `Score = (number of reference pairs that also appear in candidate) / (number of reference pairs)`

So it measures **how well the stems in C mimic those in A and B**.

**Step 7.5 - One-vote veto**

```python
STEM_THRESHOLD = 0.70

if max(A_stem_score1, A_stem_score2, A_stem_score3) < STEM_THRESHOLD \
   or max(B_stem_score1, B_stem_score2, B_stem_score3) < STEM_THRESHOLD:
    return None
```

So **both** A and B stem regions must have **at least one structure** with stem score ≥ 0.7, otherwise the candidate is discarded.

**Step 8 - Aggregate scores**

```python
avg_A_stem = mean(A_stem_scores)  # average of 3 A stem scores
avg_B_stem = mean(B_stem_scores)  # average of 3 B stem scores

stem_score = (avg_A_stem + avg_B_stem) / 2.0
loop_score = (Btop + Bbot + Atop + Abot) / 4.0

A_stem_score = avg_A_stem  # stored in output
B_stem_score = avg_B_stem  # stored in output

total = 1.5 * loop_score + energy_score - penalty + 1.5 * stem_score
```

So:
- **Loop accessibility** and **stem similarity** are each weighted by **1.5**.
- **Energy** contributes positively.
- **Sequence penalties** for GC and homopolymers subtract from the score.

The function returns:

```python
{
    "seq": C, "GC": gc, "max_run": max_run,
    "Btop_dot": Btop, "Bbelow_dot": Bbot,
    "Atop_dot": Atop, "Abelow_dot": Abot,
    "Emean": Emean, "penalty": penalty,
    "A_stem_score": A_stem_score, "B_stem_score": B_stem_score,
    "Score": total,
    # plus all coordinates a1,a2,...,B2m
}
```

#### 4.3 How GA uses these scores

Per GA round:

1. Sort population by `Score` descending.
2. Select `parents_df` = top `target_pop_size` as parent pool.
3. Compute `cheap_w = cheap_score(seq) + 0.05` for parents and use those weights for roulette sampling in mutations and crossovers.
4. Generate a `candidate_pool` (size `pool_size`) by:
   - If population size >= 2 and `random < p_crossover`: attempt crossover.
   - Otherwise: mutation via `mutate_in_targets`.
   - Each child filtered by `cheap_filter`, then assigned a `cheap_score`.
5. Sort `candidate_pool` by `cheap_score`, keep top `top_percent` fraction; evaluate them with **full `Phy_mark`** (in parallel).
6. Collect those with non-None result (they have a valid `Score`), add to population, drop duplicate sequences, sort by `Score`, truncate to `target_pop_size`.

This is the closed loop connecting **mutation / crossover**, **coordinate tracking**, and **scoring**.

---

## License

This project is provided as-is for research and educational purposes.

## Contact

For questions or contributions, please open an issue on GitHub.
