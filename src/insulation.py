"""
insulation.py

Stage 1: insulation-score breakpoint detection (Crane et al. 2015),
the standard, existing method this project uses as its fast candidate-
boundary filter. Unchanged from current field practice -- the actual
contribution of this project is Stage 2 (persistence.py), which
refines these candidates. This file's job is just to produce the
same kind of boundary calls any standard tool (cooltools, FAN-C)
would, so Stage 2 and the benchmark have something real to compare
against.

Validated against synthetic data with known planted boundaries before
being pointed at real data -- same discipline as synthetic_validation.py
for Stage 2. See JOURNAL.md for the validation result.
"""

import numpy as np

from io_utils import load_hic, get_matrix


def insulation_score(matrix: np.ndarray, window: int) -> np.ndarray:
    """
    Compute the insulation score at each bin: the mean contact value
    within a window x window square that straddles the diagonal at
    that position (window bins upstream, window bins downstream).
    This is low at a real boundary (few contacts cross it) and high
    within a domain (many contacts on both sides freely mix).

    Edge bins that don't have a full window on both sides are left as
    NaN rather than computed from a truncated, biased window.
    """
    n = matrix.shape[0]
    scores = np.full(n, np.nan)
    for i in range(window, n - window):
        block = matrix[i - window:i, i:i + window]
        scores[i] = block.mean()
    return scores


def call_boundaries(scores: np.ndarray, min_prominence: float = 0.1):
    """
    Call candidate boundaries as local minima of the log2 insulation
    score, following the delta-vector approach in Crane et al.: a
    boundary is a bin where the score dips below both neighbors, kept
    only if the dip is prominent enough (min_prominence) to not just
    be noise. Prominence here is the insulation-score equivalent of
    persistence in Stage 2 -- both are a magnitude attached to a call,
    though only persistence (Stage 2) comes with a formal robustness
    guarantee (see the paper, Section 5).

    Returns a list of (bin_index, prominence) tuples, unsorted.
    """
    log_scores = np.log2(scores + 1e-9)
    valid = ~np.isnan(log_scores)
    idx = np.where(valid)[0]
    vals = log_scores[idx]

    boundaries = []
    for k in range(1, len(vals) - 1):
        if vals[k] < vals[k - 1] and vals[k] < vals[k + 1]:
            prominence = min(vals[k - 1], vals[k + 1]) - vals[k]
            if prominence >= min_prominence:
                boundaries.append((idx[k], prominence))
    return boundaries


def call_boundaries_multiscale(matrix: np.ndarray, windows: list, min_prominence: float = 0.1):
    """
    Run Stage 1 at multiple window sizes, matching the paper's
    description (Section 4.1) of candidate boundaries at multiple
    scales, which Stage 2 then examines per-domain. Returns a dict
    of {window_size: [(bin_index, prominence), ...]}.
    """
    results = {}
    for w in windows:
        scores = insulation_score(matrix, w)
        results[w] = call_boundaries(scores, min_prominence)
    return results


if __name__ == "__main__":
    hic = load_hic()

    chrom = "21"
    resolution = 10000
    start, end = 20_000_000, 22_000_000

    print(f"Fetching chr{chrom}:{start}-{end} at {resolution} bp...")
    try:
        matrix = get_matrix(hic, chrom, start, end, resolution, normalization="KR")
    except Exception:
        matrix = get_matrix(hic, chrom, start, end, resolution, normalization="NONE")

    window_bins = 8  # 80kb window at 10kb resolution, a reasonable single-scale start
    scores = insulation_score(matrix, window_bins)
    boundaries = call_boundaries(scores, min_prominence=0.15)
    boundaries.sort(key=lambda x: -x[1])

    print(f"\nFound {len(boundaries)} candidate boundaries (window={window_bins} bins)")
    print("Top 10 by prominence:")
    for bin_idx, prom in boundaries[:10]:
        genomic_pos = start + bin_idx * resolution
        print(f"  bin={bin_idx}  position={genomic_pos:,} bp  prominence={prom:.3f}")
