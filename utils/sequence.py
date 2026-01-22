def gc_count(x: str) -> float:
    """Calculate GC content"""
    x = x.upper()
    if not x:
        return 0.0
    gc = x.count("G") + x.count("C")
    return gc / len(x)


def max_same_length(x: str) -> int:
    """Maximum consecutive length of the same base"""
    if not x:
        return 0
    max_run = 1
    run = 1
    last = x[0]
    for c in x[1:]:
        if c == last:
            run += 1
            max_run = max(max_run, run)
        else:
            last = c
            run = 1
    return max_run


def rev_comp(seq: str) -> str:
    """Reverse complement sequence (handles RNA U)"""
    tbl = str.maketrans("ACGUacgu", "UGCAugca")
    return seq.translate(tbl)[::-1]

