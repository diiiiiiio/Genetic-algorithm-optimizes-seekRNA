from statistics import mean
from utils.sequence import gc_count, max_same_length
from utils.coords import map_interval
from utils.folding import Fold
from config_loader import TARGET_SITES, MOTIF_LENGTHS

# Load constants from config
j1, j2 = TARGET_SITES["j1"], TARGET_SITES["j2"]
J1, J2 = TARGET_SITES["J1"], TARGET_SITES["J2"]
k1, k2 = TARGET_SITES["k1"], TARGET_SITES["k2"]
K1, K2 = TARGET_SITES["K1"], TARGET_SITES["K2"]
A_MOTIF_LEN = MOTIF_LENGTHS["A_MOTIF_LEN"]
B_MOTIF_LEN = MOTIF_LENGTHS["B_MOTIF_LEN"]

def dot_ratio(S: str, a: int, b: int) -> float:
    """Proportion of '.' in [a, b) interval of a single structure string"""
    x = S[a:b]
    return x.count('.') / len(x) if x else 0.0


def dot_ratio_ensemble(stru_C, a: int, b: int) -> float:
    """Average loop accessibility across multiple structures"""
    if not stru_C:
        return 0.0
    return mean(dot_ratio(S, a, b) for S in stru_C)


def parse_pairs(struct_str: str):
    """Convert dot-bracket structure to a list of base pairs (i,j)"""
    stack = []
    pairs = []
    for i, ch in enumerate(struct_str):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if stack:
                j = stack.pop()
                pairs.append((j, i))
    pairs.sort()
    return pairs


def pair_distance_score(ref_struct: str, cand_struct: str) -> float:
    """Compare similarity of pairing patterns between two structure strings"""
    if not ref_struct or not cand_struct:
        return 0.0

    L = min(len(ref_struct), len(cand_struct))
    ref = ref_struct[:L]
    cand = cand_struct[:L]

    pairs_ref = parse_pairs(ref)
    if not pairs_ref:
        return 0.0

    pairs_cand = set(parse_pairs(cand))

    keep = 0
    for p, q in pairs_ref:
        if (p, q) in pairs_cand:
            keep += 1

    return keep / len(pairs_ref)


def cheap_filter(C, coords=None):
    """Tier 1: Fast sequence-only evaluation"""
    if max_same_length(C) > 7:
        return False
    gc = gc_count(C)
    if gc < 0.25 or gc > 0.75:
        return False
    
    if coords:
        n = len(C)
        required_keys = ["a1", "a2", "A1", "A2", "b1", "b2", "B1", "B2"]
        if not all(k in coords for k in required_keys):
            return False
        
        segs = [(coords[k1], coords[k2]) for k1, k2 in zip(required_keys[0::2], required_keys[1::2])]
        for s, e in segs:
            if not (0 <= s < e <= n):
                return False
    return True


def cheap_score(C, coords):
    """Tier 1: Fast sequence-only scoring"""
    gc = gc_count(C)
    max_run = max_same_length(C)
    
    score = 1.0
    if gc < 0.30 or gc > 0.70:
        score -= 0.3
    if max_run > 4:
        score -= 0.2 * (max_run - 4)
    if 0.40 <= gc <= 0.60:
        score += 0.1
    return score


def Phy_mark(C, struct_ref_A, struct_ref_B, i, coords=None, fold_mode="full"):
    """Tier 2: Comprehensive evaluation combining structure, energy, and sequence features"""
    L_B_ref = len(struct_ref_B[0]) if struct_ref_B else 0
    
    # 1) Coordinate determination
    if coords is None:
        a1, a2 = map_interval(j1, j2, i, L_B_ref)
        A1, A2 = map_interval(J1, J2, i, L_B_ref)
        if j2 <= i <= J1 or a1 is None or A1 is None:
            return None

        b1, b2 = i + k1, i + k2
        B1, B2 = i + K1, i + K2

        a1m, a2m = map_interval(j2 - A_MOTIF_LEN, j2, i, L_B_ref)
        A1m, A2m = map_interval(J2 - A_MOTIF_LEN, J2, i, L_B_ref)
        if any(x is None for x in (a1m, a2m, A1m, A2m)):
            return None

        b1m, b2m = i + (k2 - B_MOTIF_LEN), i + k2
        B1m, B2m = i + (K2 - B_MOTIF_LEN), i + K2
    else:
        a1, a2, A1, A2 = coords["a1"], coords["a2"], coords["A1"], coords["A2"]
        b1, b2, B1, B2 = coords["b1"], coords["b2"], coords["B1"], coords["B2"]
        a1m, a2m, A1m, A2m = coords["a1m"], coords["a2m"], coords["A1m"], coords["A2m"]
        b1m, b2m, B1m, B2m = coords["b1m"], coords["b2m"], coords["B1m"], coords["B2m"]
        
        n = len(C)
        segs = [(a1, a2), (A1, A2), (b1, b2), (B1, B2), (a1m, a2m), (A1m, A2m), (b1m, b2m), (B1m, B2m)]
        if any(s is None or e is None or not (0 <= s < e <= n) for s, e in segs):
            return None

    # 2) Hard filter
    if not cheap_filter(C, coords):
        return None

    # 3) Folding
    stru_C, E = Fold(C, mode=fold_mode)
    if len(stru_C) < 1:
        return None
    
    # Pad stru_C if needed for indices
    while len(stru_C) < 3:
        stru_C.append(stru_C[0])

    # 4) Loop accessibility
    Btop = dot_ratio_ensemble(stru_C, b1, b2)
    Bbot = dot_ratio_ensemble(stru_C, B1, B2)
    Atop = dot_ratio_ensemble(stru_C, a1, a2)
    Abot = dot_ratio_ensemble(stru_C, A1, A2)

    # 5) Energy
    Emean = mean(E) if E else 0.0
    energy_score = (Emean + 1.0) / 2.0

    # 6) Sequence penalty
    gc = gc_count(C)
    max_run = max_same_length(C)
    penalty = 0.0
    if gc < 0.30 or gc > 0.70: penalty += 0.2
    if max_run > 4: penalty += 0.5

    # 7) Stem similarity
    n = len(C)
    # A region stem
    A_key_start, A_key_end = max(0, a2), min(n, A1)
    A_ref_len = min(max(0, A_key_end - A_key_start), J1 - j2)
    
    A_stem_scores = []
    for i_stru in range(3):
        cand_s = stru_C[i_stru][A_key_start : A_key_start + A_ref_len]
        ref_s = struct_ref_A[i_stru][j2 : j2 + A_ref_len]
        A_stem_scores.append(pair_distance_score(ref_s, cand_s))
    
    # B region stem
    B_key_start, B_key_end = max(0, b2), min(n, B1)
    B_ref_len = min(max(0, B_key_end - B_key_start), K1 - k2)
    
    B_stem_scores = []
    for i_stru in range(3):
        cand_s = stru_C[i_stru][B_key_start : B_key_start + B_ref_len]
        ref_s = struct_ref_B[i_stru][k2 : k2 + B_ref_len]
        B_stem_scores.append(pair_distance_score(ref_s, cand_s))

    # One-vote veto
    STEM_THRESHOLD = 0.70
    if max(A_stem_scores) < STEM_THRESHOLD or max(B_stem_scores) < STEM_THRESHOLD:
        return None

    # 8) Final Score
    avg_A_stem = mean(A_stem_scores)
    avg_B_stem = mean(B_stem_scores)
    stem_score = (avg_A_stem + avg_B_stem) / 2.0
    loop_score = (Btop + Bbot + Atop + Abot) / 4.0

    total = 1.5 * loop_score + energy_score - penalty + 1.5 * stem_score

    return {
        "seq": C, "GC": gc, "max_run": max_run,
        "Btop_dot": Btop, "Bbelow_dot": Bbot, "Atop_dot": Atop, "Abelow_dot": Abot,
        "Emean": Emean, "penalty": penalty,
        "A_stem_score": avg_A_stem, "B_stem_score": avg_B_stem, "Score": total,
        "a1": a1, "a2": a2, "A1": A1, "A2": A2,
        "b1": b1, "b2": b2, "B1": B1, "B2": B2,
        "a1m": a1m, "a2m": a2m, "A1m": A1m, "A2m": A2m,
        "b1m": b1m, "b2m": b2m, "B1m": B1m, "B2m": B2m,
    }

