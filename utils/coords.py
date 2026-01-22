def map_interval(start: int, end: int, i_insert: int, L_insert: int):
    """
    Map the [start, end) interval on A to coordinates on C = A[:i] + B + A[i:]:
      - If insertion point is to the left of or at the left end (i <= start): shift entire segment right by L_insert
      - If insertion point is to the right of or at the right end (i >= end): coordinates unchanged
      - If insertion point falls inside the interval (start < i < end): segment is broken by B, return (None, None)
    """
    if i_insert <= start:
        return start + L_insert, end + L_insert
    elif i_insert >= end:
        return start, end
    else:
        return None, None


def update_interval_after_edit(start: int, end: int, cut_start: int, cut_end: int, delta: int):
    """Update interval coordinates after an edit (substitution, insertion, deletion)"""
    if start > cut_end:
        return start + delta, end + delta
    if end <= cut_start:
        return start, end
    if delta > 0:
        return start, end + delta
    removed = max(0, min(end, cut_end) - max(start, cut_start))
    new_start = start
    new_end = end - removed
    if cut_start < start:
        shift = min(-delta, start - cut_start)
        new_start -= shift
        new_end -= shift
    return new_start, new_end


def update_coords_after_edit(coords, cut_start, cut_end, delta, seqlen_after):
    """Update all recognition and motif coordinates after a sequence edit"""
    keys = ["a1", "a2", "A1", "A2", "b1", "b2", "B1", "B2",
            "a1m", "a2m", "A1m", "A2m", "b1m", "b2m", "B1m", "B2m"]
    new_coords = dict(coords)
    for s_key, e_key in zip(keys[0::2], keys[1::2]):
        s = new_coords[s_key]
        e = new_coords[e_key]
        s2, e2 = update_interval_after_edit(s, e, cut_start, cut_end, delta)
        s2 = max(0, min(s2, seqlen_after))
        e2 = max(s2+1, min(e2, seqlen_after))
        new_coords[s_key], new_coords[e_key] = s2, e2
    return new_coords


def sanity_check_coords(seq, coords):
    """Verify that recognition regions are within the current sequence bounds"""
    n = len(seq)
    for s, e in [(coords["a1"], coords["a2"]), (coords["A1"], coords["A2"]),
                 (coords["b1"], coords["b2"]), (coords["B1"], coords["B2"])]:
        if not (0 <= s < e <= n):
            return False
    return True

