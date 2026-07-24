"""
real_data_persistence.py

Applies the corrected Stage 2 persistent homology method (validated on
synthetic data in synthetic_validation.py -- diagonal masked to avoid
the corner-connection artifact) to a real Hi-C contact matrix.

This is exploratory: does real data show the same kind of clean
persistence gap (real domains vs. noise) that the synthetic test did,
or is real biological structure messier? No claim is made here about
which boundaries are "correct" -- there's no ground truth for real
data, that's the whole reason we validated on synthetic data first.
"""

import numpy as np
import gudhi
import matplotlib.pyplot as plt

from io_utils import load_hic, get_matrix


def mask_diagonal(matrix: np.ndarray) -> np.ndarray:
    """
    Mask the exact diagonal (self-contacts) to the matrix's background
    level, same as the synthetic validation -- this is standard Hi-C
    practice (self-contacts are an assay artifact) and avoids the
    trivial corner-connection issue in the cubical filtration.
    """
    masked = matrix.copy()
    off_diagonal_values = matrix[~np.eye(matrix.shape[0], dtype=bool)]
    background_level = np.median(off_diagonal_values[off_diagonal_values > 0])
    np.fill_diagonal(masked, background_level)
    return masked


def mask_long_range(matrix: np.ndarray, max_distance_bins: int) -> np.ndarray:
    """
    Mask out any bin-pair farther apart than max_distance_bins, setting
    it to the matrix median. This is the actual fix for the corner-noise
    artifact: the far corners of ANY square crop are the longest-range
    pairs WITHIN that crop, regardless of crop size, and O/E ratios at
    long range are dominated by sampling noise (expected count shrinks
    toward zero, so a single stray read inflates the ratio hugely).
    Real TAD structure lives well within 1Mb; capping distance directly,
    rather than relying on crop size, is what actually removes this.
    """
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
    h0 = cc.persistence_intervals_in_dimension(0)
    pairs = cc.cofaces_of_persistence_pairs()
    finite_pairs = pairs[0][0] if len(pairs[0]) > 0 else np.empty((0, 2), dtype=int)
    return h0, finite_pairs


def summarize_and_locate(matrix, h0, finite_pairs, bin_start, resolution, top_n=10):
    finite_mask = np.isfinite(h0[:, 1])
    finite_intervals = h0[finite_mask]
    persistences = finite_intervals[:, 1] - finite_intervals[:, 0]
    order = np.argsort(-persistences)

    print(f"Total H0 features: {len(h0)}  (finite: {len(finite_intervals)})")
    print(f"\nTop {top_n} most persistent features:")
    print("(birth cell location is a representative point INSIDE a persistent")
    print(" domain, not a precise boundary -- domains aren't perfectly uniform)")

    results = []
    for i in order[:top_n]:
        birth_idx, death_idx = finite_pairs[i]
        row, col = np.unravel_index(birth_idx, matrix.shape)
        genomic_row = bin_start + row * resolution
        genomic_col = bin_start + col * resolution
        pers = persistences[i]
        print(f"  persistence={pers:.2f}  representative point near "
              f"{genomic_row:,}-{genomic_col:,} bp")
        results.append((pers, row, col))

    if len(persistences) > top_n:
        top = np.sort(persistences)[-top_n:]
        rest = np.sort(persistences)[:-top_n]
        gap = top.min() - rest.max()
        print(f"\nGap between top {top_n} and the rest: {gap:.2f}")
        print("(On synthetic data this gap was ~8. A much smaller or negative")
        print(" gap here would mean real data doesn't separate as cleanly --")
        print(" itself a real, useful finding, not a failure.)")

    return results


def plot_with_detections(matrix, results, title, out_path):
    plt.figure(figsize=(7, 7))
    plt.imshow(np.log1p(matrix), cmap="Reds", origin="upper")
    for pers, row, col in results:
        plt.scatter(col, row, s=80, facecolors="none", edgecolors="blue", linewidths=1.5)
    plt.title(title)
    plt.colorbar(label="log(1+contact)")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"\nSaved {out_path}")


def plot_persistence_scree(h0, out_path):
    """
    Plot all H0 persistence values sorted descending, looking for a
    natural elbow -- a point where values drop steeply then flatten
    into noise. This is a fairer diagnostic on real data than a fixed
    top-N cutoff, since we don't know in advance how many real domains
    exist in the window (unlike the synthetic test, where we planted
    an exact, known number).
    """
    import matplotlib.pyplot as plt
    finite = h0[np.isfinite(h0[:, 1])]
    persistences = np.sort(finite[:, 1] - finite[:, 0])[::-1]

    plt.figure(figsize=(8, 5))
    plt.plot(persistences[:200], marker='o', markersize=3)
    plt.xlabel("Rank (most to least persistent)")
    plt.ylabel("Persistence")
    plt.title("Sorted persistence values -- look for a natural elbow")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved {out_path}")


if __name__ == "__main__":
    hic = load_hic()

    chrom = "21"
    resolution = 10000
    start, end = 20_000_000, 22_000_000

    print(f"Fetching chr{chrom}:{start}-{end} at {resolution} bp...")
    try:
        matrix = get_matrix(hic, chrom, start, end, resolution, normalization="KR", data_type="oe")
    except Exception:
        matrix = get_matrix(hic, chrom, start, end, resolution, normalization="NONE", data_type="oe")

    print(f"Matrix shape: {matrix.shape}\n")

    masked = mask_diagonal(matrix)
    masked = mask_long_range(masked, max_distance_bins=30)  # cap at 30 bins = 300kb, closer to single-TAD scale  # cap at 100 bins = 1Mb at 10kb resolution
    h0, finite_pairs = run_cubical_persistence(masked)
    results = summarize_and_locate(masked, h0, finite_pairs, start, resolution)
    plot_persistence_scree(h0, "figures/persistence_scree.png")

    plot_with_detections(
        matrix, results,
        f"chr{chrom}:{start}-{end}, GM12878, {resolution}bp -- top persistence points",
        "figures/real_data_persistence.png"
    )


