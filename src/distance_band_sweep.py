"""
distance_band_sweep.py

Tests whether the ~300kb distance-band cutoff (used in
real_data_persistence.py to fix the long-range noise artifact) is a
stable, defensible choice -- or whether it was just what happened to
look reasonable on one region. This directly addresses paper
Limitation 6, which states the cutoff was chosen by direct
observation on a single region and needs a sensitivity check before
being treated as validated.

Same region as before (chr21:20-22Mb, GM12878, 10kb) is swept across
several candidate cutoffs. If the resulting scree curves look similar
in shape (elbow in a similar place, similar overall structure) across
a reasonable range of cutoffs, that's real evidence 300kb wasn't an
arbitrary lucky pick. If the curves change shape substantially with
the cutoff, that's the opposite finding -- also worth knowing and
reporting honestly.
"""

import numpy as np
import gudhi
import matplotlib.pyplot as plt

from io_utils import load_hic, get_matrix


def mask_diagonal(matrix: np.ndarray) -> np.ndarray:
    masked = matrix.copy()
    off_diagonal_values = matrix[~np.eye(matrix.shape[0], dtype=bool)]
    background_level = np.median(off_diagonal_values[off_diagonal_values > 0])
    np.fill_diagonal(masked, background_level)
    return masked


def mask_long_range(matrix: np.ndarray, max_distance_bins: int) -> np.ndarray:
    n = matrix.shape[0]
    idx = np.arange(n)
    dist = np.abs(idx[:, None] - idx[None, :])
    masked = matrix.copy()
    background = np.median(matrix[dist <= max_distance_bins])
    masked[dist > max_distance_bins] = background
    return masked


def run_cubical_persistence(matrix: np.ndarray):
    filtration_values = -matrix
    cc = gudhi.CubicalComplex(top_dimensional_cells=filtration_values)
    cc.compute_persistence()
    return cc.persistence_intervals_in_dimension(0)


def sorted_persistences(h0):
    finite = h0[np.isfinite(h0[:, 1])]
    return np.sort(finite[:, 1] - finite[:, 0])[::-1]


if __name__ == "__main__":
    hic = load_hic()

    chrom = "21"
    resolution = 10000
    start, end = 20_000_000, 22_000_000

    print(f"Fetching chr{chrom}:{start}-{end} at {resolution} bp (O/E)...")
    try:
        matrix = get_matrix(hic, chrom, start, end, resolution, normalization="KR", data_type="oe")
    except Exception:
        matrix = get_matrix(hic, chrom, start, end, resolution, normalization="NONE", data_type="oe")

    cutoffs_kb = [100, 200, 300, 500, 1000]
    cutoffs_bins = [int(kb * 1000 / resolution) for kb in cutoffs_kb]

    plt.figure(figsize=(9, 6))
    summary = []

    for kb, bins in zip(cutoffs_kb, cutoffs_bins):
        masked = mask_diagonal(matrix)
        masked = mask_long_range(masked, bins)
        h0 = run_cubical_persistence(masked)
        pers = sorted_persistences(h0)

        plt.plot(pers[:100], marker='o', markersize=2, label=f"{kb}kb cutoff")

        top10 = pers[:10]
        rest = pers[10:60] if len(pers) > 60 else pers[10:]
        gap = (top10.min() - rest.max()) if len(rest) > 0 else float("nan")
        summary.append((kb, pers[0], gap))
        print(f"cutoff={kb}kb  top_persistence={pers[0]:.3f}  "
              f"gap(top10 vs next50)={gap:.3f}  n_features={len(pers)}")

    plt.xlabel("Rank (most to least persistent)")
    plt.ylabel("Persistence")
    plt.title("Distance-band cutoff sensitivity -- chr21:20-22Mb, GM12878, 10kb")
    plt.legend()
    plt.tight_layout()
    plt.savefig("figures/distance_band_sweep.png", dpi=150)
    plt.close()
    print("\nSaved figures/distance_band_sweep.png")

    print("\nSummary:")
    print(f"{'cutoff (kb)':>12} {'top persistence':>16} {'gap':>10}")
    for kb, top, gap in summary:
        print(f"{kb:>12} {top:>16.3f} {gap:>10.3f}")
    print("\nIf top persistence and gap are roughly stable across cutoffs from")
    print("~200-500kb, that's real support for 300kb as a reasonable choice,")
    print("not an arbitrary one. Large swings would mean the opposite --")
    print("also a genuine, reportable finding, not a failure.")
