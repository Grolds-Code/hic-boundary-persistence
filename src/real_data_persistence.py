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


if __name__ == "__main__":
    hic = load_hic()

    chrom = "21"
    resolution = 10000
    start, end = 20_000_000, 25_000_000

    print(f"Fetching chr{chrom}:{start}-{end} at {resolution} bp...")
    try:
        matrix = get_matrix(hic, chrom, start, end, resolution, normalization="KR")
    except Exception:
        matrix = get_matrix(hic, chrom, start, end, resolution, normalization="NONE")

    print(f"Matrix shape: {matrix.shape}\n")

    masked = mask_diagonal(matrix)
    h0, finite_pairs = run_cubical_persistence(masked)
    results = summarize_and_locate(masked, h0, finite_pairs, start, resolution)

    plot_with_detections(
        matrix, results,
        f"chr{chrom}:{start}-{end}, GM12878, {resolution}bp -- top persistence points",
        "figures/real_data_persistence.png"
    )
