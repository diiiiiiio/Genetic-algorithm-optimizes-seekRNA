import random
import numpy as np
from utils.sequence import rev_comp
from utils.coords import update_coords_after_edit, sanity_check_coords
from config_loader import NUCS

def apply_edit(seq: str, pos: int, edit_type: str, nts: str = "", del_len: int = 1):
    """Perform one edit on seq (sub/ins/del)"""
    n = len(seq)
    pos = max(0, min(pos, n))

    if edit_type == "sub":
        if not nts:
            nts = random.choice(NUCS)
        old = seq[pos]
        if nts == old:
            nts = random.choice([b for b in NUCS if b != old])
        new_seq = seq[:pos] + nts + seq[pos+1:]
        return new_seq, pos, pos+1, 0

    elif edit_type == "ins":
        if not nts:
            k = random.choice([1, 2])
            nts = "".join(random.choices(NUCS, k=k))
        new_seq = seq[:pos] + nts + seq[pos:]
        return new_seq, pos, pos, len(nts)

    elif edit_type == "del":
        k = max(1, del_len)
        if pos + k > n:
            k = max(1, n - pos)
        new_seq = seq[:pos] + seq[pos+k:]
        return new_seq, pos, pos+k, -k

    return seq, pos, pos, 0


def mutate_single_segment(seq, coords, seg_key_pair, motif_key_pair, rng, p_sub=0.6, p_ins=0.2, p_del=0.2):
    """Mutate within a single recognition region, avoiding the motif"""
    s, e = coords[seg_key_pair[0]], coords[seg_key_pair[1]]
    if e - s <= 1: return seq, coords

    m_start, m_end = coords[motif_key_pair[0]], coords[motif_key_pair[1]]
    candidates = [p for p in range(s, e) if not (m_start <= p < m_end)]
    if not candidates: return seq, coords

    pos = int(rng.choice(candidates))
    r = rng.random()
    total_p = p_sub + p_ins + p_del
    
    if r < p_sub / total_p:
        old_nt = seq[pos]
        nts = rng.choice([n for n in NUCS if n != old_nt])
        new_seq, cut_s, cut_e, delta = apply_edit(seq, pos, "sub", nts=nts)
    elif r < (p_sub + p_ins) / total_p:
        ins_len = int(rng.integers(1, 3))
        nts = "".join(rng.choice(list(NUCS), size=ins_len))
        new_seq, cut_s, cut_e, delta = apply_edit(seq, pos, "ins", nts=nts)
    else:
        new_seq, cut_s, cut_e, delta = apply_edit(seq, pos, "del", del_len=1)

    new_coords = update_coords_after_edit(coords, cut_s, cut_e, delta, len(new_seq))
    return new_seq, new_coords


def mutate_paired_segments(seq, coords, seg_top_pair, seg_bot_pair, motif_top_pair, motif_bot_pair, rng, p_sub=0.6, p_ins=0.2, p_del=0.2):
    """Perform paired mutations within the motif region of complementary strands"""
    t_m1, t_m2 = coords[motif_top_pair[0]], coords[motif_top_pair[1]]
    b_m1, b_m2 = coords[motif_bot_pair[0]], coords[motif_bot_pair[1]]
    L_motif = t_m2 - t_m1
    if L_motif <= 0: return seq, coords

    r = rng.random()
    total_p = p_sub + p_ins + p_del
    
    if r < p_sub / total_p:
        edit_type = "sub"
        offset = int(rng.integers(0, L_motif))
    elif r < (p_sub + p_ins) / total_p:
        edit_type = "ins"
        offset = int(rng.integers(1, L_motif - 1)) if L_motif > 2 else int(rng.integers(0, L_motif))
    else:
        edit_type = "del"
        offset = int(rng.integers(0, L_motif))

    t_pos = t_m1 + offset

    if edit_type == "sub":
        new_nt = rng.choice([n for n in NUCS if n != seq[t_pos]])
        new_seq, cut_s, cut_e, delta = apply_edit(seq, t_pos, "sub", nts=new_nt)
        new_coords = update_coords_after_edit(coords, cut_s, cut_e, delta, len(new_seq))
        
        b_pos = new_coords[motif_bot_pair[1]] - 1 - offset
        new_seq2, cut_s2, cut_e2, delta2 = apply_edit(new_seq, b_pos, "sub", nts=rev_comp(new_nt))
        return new_seq2, update_coords_after_edit(new_coords, cut_s2, cut_e2, delta2, len(new_seq2))

    elif edit_type == "ins":
        ins_len = int(rng.integers(1, 3))
        nts = "".join(rng.choice(list(NUCS), size=ins_len))
        new_seq, cut_s, cut_e, delta = apply_edit(seq, t_pos, "ins", nts=nts)
        new_coords = update_coords_after_edit(coords, cut_s, cut_e, delta, len(new_seq))
        
        b_pos = new_coords[motif_bot_pair[1]] - offset
        new_seq2, cut_s2, cut_e2, delta2 = apply_edit(new_seq, b_pos, "ins", nts=rev_comp(nts))
        return new_seq2, update_coords_after_edit(new_coords, cut_s2, cut_e2, delta2, len(new_seq2))

    else: # del
        new_seq, cut_s, cut_e, delta = apply_edit(seq, t_pos, "del", del_len=1)
        new_coords = update_coords_after_edit(coords, cut_s, cut_e, delta, len(new_seq))
        
        b_pos = new_coords[motif_bot_pair[1]] - 1 - offset
        new_seq2, cut_s2, cut_e2, delta2 = apply_edit(new_seq, b_pos, "del", del_len=1)
        return new_seq2, update_coords_after_edit(new_coords, cut_s2, cut_e2, delta2, len(new_seq2))


def mutate_in_targets(seq, coords, n_events=1, rng=None, p_pair=0.7, single_probs=(0.6, 0.2, 0.2), pair_probs=(0.6, 0.2, 0.2)):
    """Main mutation entry point: selects between single-strand and paired mutations"""
    if rng is None: rng = np.random.default_rng()
    seq_cur, coords_cur = seq, dict(coords)

    for _ in range(n_events):
        if not sanity_check_coords(seq_cur, coords_cur): break

        if rng.random() < p_pair:
            choice = rng.choice(["A_pair", "B_pair"])
            if choice == "A_pair":
                seq_cur, coords_cur = mutate_paired_segments(seq_cur, coords_cur, ("a1", "a2"), ("A1", "A2"), ("a1m", "a2m"), ("A1m", "A2m"), rng, *pair_probs)
            else:
                seq_cur, coords_cur = mutate_paired_segments(seq_cur, coords_cur, ("b1", "b2"), ("B1", "B2"), ("b1m", "b2m"), ("B1m", "B2m"), rng, *pair_probs)
        else:
            choice = rng.choice(["Atop", "Abelow", "Btop", "Bbelow"])
            mappings = {
                "Atop": (("a1", "a2"), ("a1m", "a2m")),
                "Abelow": (("A1", "A2"), ("A1m", "A2m")),
                "Btop": (("b1", "b2"), ("b1m", "b2m")),
                "Bbelow": (("B1", "B2"), ("B1m", "B2m"))
            }
            seg, motif = mappings[choice]
            seq_cur, coords_cur = mutate_single_segment(seq_cur, coords_cur, seg, motif, rng, *single_probs)

    return seq_cur, coords_cur


def crossover(p1, p2):
    """One-point crossover, only if parents have identical coordinates"""
    s1, s2 = p1["seq"], p2["seq"]
    if len(s1) != len(s2) or p1["i"] != p2["i"]: return None

    keys = ["a1","a2","A1","A2","b1","b2","B1","B2","a1m","a2m","A1m","A2m","b1m","b2m","B1m","B2m"]
    if any(int(p1[k]) != int(p2[k]) for k in keys): return None

    n = len(s1)
    if n <= 1: return None
    cut = random.randint(1, n - 1)
    child_seq = s1[:cut] + s2[cut:]
    return child_seq, {k: int(p1[k]) for k in keys}, int(p1["i"])

