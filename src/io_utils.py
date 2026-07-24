"""
io_utils.py

Thin wrapper around hicstraw for loading .hic files and pulling
contact matrices for a genomic region at a given resolution.

This is the first, purely exploratory script for the project:
load one real .hic file, list what's in it, and pull one small
region as a dense matrix so we can actually look at real data
before writing any insulation-score or persistent-homology code. Something of the sort..
"""

import numpy as np
import hicstraw


# Public GM12878 in-situ Hi-C file (Rao et al. 2014), combined replicates,
# MAPQ >= 30. Hosted by the Aiden lab; straw can stream directly from
# this URL without downloading the whole file first.
DEFAULT_HIC_URL = "https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined_30.hic"


def load_hic(path_or_url: str = DEFAULT_HIC_URL) -> hicstraw.HiCFile:
    """
    Open a .hic file, local path or URL.
    """
    hic = hicstraw.HiCFile(path_or_url)
    return hic


def describe(hic: hicstraw.HiCFile) -> None:
    """
    Print basic info about a loaded .hic file: genome build,
    chromosomes, and available resolutions.
    """
    print(f"Genome ID: {hic.getGenomeID()}")

    chroms = hic.getChromosomes()
    print(f"Chromosomes ({len(chroms)}):")
    for c in chroms:
        print(f"  {c.name}\tlength={c.length}")

    resolutions = hic.getResolutions()
    print(f"Available resolutions (bp): {resolutions}")


def get_matrix(
    hic: hicstraw.HiCFile,
    chrom: str,
    start: int,
    end: int,
    resolution: int,
    normalization: str = "KR",
    data_type: str = "observed",
) -> np.ndarray:
    """
    Fetch a dense contact matrix for a single genomic region on one
    chromosome, at a given resolution.

    chrom: chromosome name as it appears in hic.getChromosomes(), e.g. "21"
    start, end: region bounds in base pairs
    resolution: bin size in bp; must be one of hic.getResolutions()
    normalization: 'NONE', 'VC', 'VC_SQRT', 'KR', or 'SCALE'
        (must actually exist in the .hic file for this resolution --
        straw only reads what's already stored, it does not compute
        normalization itself)
    data_type: 'observed' or 'oe' (observed/expected)
    """
    mzd = hic.getMatrixZoomData(chrom, chrom, data_type, normalization, "BP", resolution)
    matrix = mzd.getRecordsAsMatrix(start, end, start, end)
    return matrix


if __name__ == "__main__":
    # Quick smoke test: load the file, list what's in it, then pull
    # a small region on chr21 to confirm everything actually works
    # end to end against real data.
    hic = load_hic()
    describe(hic)

    chrom = "21"
    resolution = 25000  # 25kb: coarse enough to load fast for a first look
    start, end = 20_000_000, 25_000_000  # a 5 Mb window

    print(f"\nFetching chr{chrom}:{start}-{end} at {resolution} bp resolution...")
    try:
        matrix = get_matrix(hic, chrom, start, end, resolution, normalization="KR")
    except Exception as e:
        print(f"KR normalization failed ({e}), retrying with NONE...")
        matrix = get_matrix(hic, chrom, start, end, resolution, normalization="NONE")

    print(f"Matrix shape: {matrix.shape}")
    print(f"Non-zero entries: {np.count_nonzero(matrix)} / {matrix.size}")
    print(f"Max value: {matrix.max():.2f}")
    print(f"Mean value (non-zero only): {matrix[matrix > 0].mean():.2f}")

def plot_matrix(matrix: np.ndarray, title: str, out_path: str, log_scale: bool = True) -> None:
    """
    Render a contact matrix as a heatmap and save it to disk.

    log_scale: Hi-C matrices span several orders of magnitude between
    the diagonal and off-diagonal, so log1p makes structure visible
    across the whole matrix instead of just a bright diagonal line.
    """
    import matplotlib.pyplot as plt

    data = np.log1p(matrix) if log_scale else matrix

    plt.figure(figsize=(7, 7))
    plt.imshow(data, cmap="Reds", origin="upper")
    plt.title(title)
    plt.colorbar(label="log(1 + contacts)" if log_scale else "contacts")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved plot to {out_path}")


def plot_resolution_comparison(hic, chrom: str, start: int, end: int, resolutions: list, out_path: str) -> None:
    """
    Fetch the same genomic region at several resolutions and plot them
    side by side, so structure (or its breakdown) across scales is
    directly visible in one image.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    n = len(resolutions)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5))
    if n == 1:
        axes = [axes]

    for ax, res in zip(axes, resolutions):
        try:
            matrix = get_matrix(hic, chrom, start, end, res, normalization="KR")
        except Exception:
            matrix = get_matrix(hic, chrom, start, end, res, normalization="NONE")
        data = np.log1p(matrix)
        ax.imshow(data, cmap="Reds", origin="upper")
        ax.set_title(f"{res} bp\n{matrix.shape[0]}x{matrix.shape[1]} bins")
        ax.axis("off")

    plt.suptitle(f"chr{chrom}:{start}-{end}, GM12878, across resolutions")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved comparison to {out_path}")
