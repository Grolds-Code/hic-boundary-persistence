import modal

app = modal.App("hic-boundary-persistence")

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("libcurl4-openssl-dev", "build-essential", "zlib1g-dev")
    .pip_install("hic-straw", "gudhi", "numpy", "scipy", "matplotlib")
)


@app.function(image=image, timeout=1800)
def run_sweep():
    import io

    import gudhi
    import hicstraw
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    def load_hic(url="https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined_30.hic"):
        return hicstraw.HiCFile(url)

    def get_matrix(hic, chrom, start, end, resolution, normalization="KR", data_type="observed"):
        mzd = hic.getMatrixZoomData(chrom, chrom, data_type, normalization, "BP", resolution)
        return mzd.getRecordsAsMatrix(start, end, start, end)

    def mask_diagonal(matrix):
        masked = matrix.copy()
        off = matrix[~np.eye(matrix.shape[0], dtype=bool)]
        bg = np.median(off[off > 0])
        np.fill_diagonal(masked, bg)
        return masked

    def mask_long_range(matrix, max_distance_bins):
        n = matrix.shape[0]
        idx = np.arange(n)
        dist = np.abs(idx[:, None] - idx[None, :])
        masked = matrix.copy()
        bg = np.median(matrix[dist <= max_distance_bins])
        masked[dist > max_distance_bins] = bg
        return masked

    def run_cubical_persistence(matrix):
        filt = -matrix
        cc = gudhi.CubicalComplex(top_dimensional_cells=filt)
        cc.compute_persistence()
        return cc.persistence_intervals_in_dimension(0)

    def sorted_persistences(h0):
        finite = h0[np.isfinite(h0[:, 1])]
        return np.sort(finite[:, 1] - finite[:, 0])[::-1]

    hic = load_hic()
    chrom, resolution = "21", 10000
    start, end = 20_000_000, 22_000_000
    matrix = get_matrix(hic, chrom, start, end, resolution, normalization="KR", data_type="oe")

    cutoffs_kb = [20, 50, 100, 200, 300, 500, 1000]
    cutoffs_bins = [int(kb * 1000 / resolution) for kb in cutoffs_kb]

    plt.figure(figsize=(9, 6))
    summary_lines = []
    for kb, bins in zip(cutoffs_kb, cutoffs_bins):
        m = mask_diagonal(matrix)
        m = mask_long_range(m, bins)
        h0 = run_cubical_persistence(m)
        pers = sorted_persistences(h0)

        plt.plot(pers[:100], marker="o", markersize=2, label=f"{kb}kb")

        top10 = pers[:10]
        rest = pers[10:60] if len(pers) > 60 else pers[10:]
        gap = (top10.min() - rest.max()) if len(rest) > 0 else float("nan")
        summary_lines.append(
            f"cutoff={kb}kb  top_persistence={pers[0]:.3f}  gap={gap:.3f}  n_features={len(pers)}"
        )

    plt.xlabel("Rank (most to least persistent)")
    plt.ylabel("Persistence")
    plt.title("Distance-band sweep (extended), chr21:20-22Mb, GM12878, 10kb")
    plt.legend()
    plt.tight_layout()

    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=150)
    plt.close()
    buf.seek(0)

    return buf.getvalue(), "\n".join(summary_lines)


@app.function(image=image, timeout=1800)
def run_zscore_comparison():
    import io

    import gudhi
    import hicstraw
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    def load_hic(url="https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined_30.hic"):
        return hicstraw.HiCFile(url)

    def get_matrix(hic, chrom, start, end, resolution, normalization="KR", data_type="observed"):
        mzd = hic.getMatrixZoomData(chrom, chrom, data_type, normalization, "BP", resolution)
        return mzd.getRecordsAsMatrix(start, end, start, end)

    def mask_diagonal(matrix):
        masked = matrix.copy()
        off = matrix[~np.eye(matrix.shape[0], dtype=bool)]
        bg = np.median(off[off > 0])
        np.fill_diagonal(masked, bg)
        return masked

    def distance_stratified_zscore(matrix):
        n = matrix.shape[0]
        idx = np.arange(n)
        dist = np.abs(idx[:, None] - idx[None, :])
        z = np.zeros_like(matrix)
        for d in range(n):
            mask = dist == d
            vals = matrix[mask]
            if len(vals) < 2:
                continue
            mean_d, std_d = vals.mean(), vals.std()
            if std_d > 0:
                z[mask] = (matrix[mask] - mean_d) / std_d
            else:
                z[mask] = 0
        return z

    def run_cubical_persistence(matrix):
        filt = -matrix
        cc = gudhi.CubicalComplex(top_dimensional_cells=filt)
        cc.compute_persistence()
        return cc.persistence_intervals_in_dimension(0)

    def sorted_persistences(h0):
        finite = h0[np.isfinite(h0[:, 1])]
        return np.sort(finite[:, 1] - finite[:, 0])[::-1]

    hic = load_hic()
    chrom, resolution = "21", 10000

    windows_mb = [1, 2, 3, 5]
    summary_lines = []
    plt.figure(figsize=(9, 6))

    for mb in windows_mb:
        start, end = 20_000_000, 20_000_000 + mb * 1_000_000
        matrix = get_matrix(hic, chrom, start, end, resolution, normalization="KR", data_type="oe")
        m = mask_diagonal(matrix)
        z = distance_stratified_zscore(m)
        h0 = run_cubical_persistence(z)
        pers = sorted_persistences(h0)

        plt.plot(pers[:100], marker="o", markersize=2, label=f"{mb}Mb window")
        summary_lines.append(
            f"window={mb}Mb  top_persistence={pers[0]:.3f}  n_features={len(pers)}"
        )

    plt.xlabel("Rank (most to least persistent)")
    plt.ylabel("Persistence (z-scored)")
    plt.title("Distance-stratified z-score: does persistence stay flat as window grows?")
    plt.legend()
    plt.tight_layout()

    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=150)
    plt.close()
    buf.seek(0)

    return buf.getvalue(), "\n".join(summary_lines)


@app.local_entrypoint()
def main():
    image_bytes, summary = run_sweep.remote()
    with open("figures/distance_band_sweep_extended.png", "wb") as f:
        f.write(image_bytes)
    print("=== Cutoff sweep ===")
    print(summary)
    print("Saved figures/distance_band_sweep_extended.png\n")

    zimage_bytes, zsummary = run_zscore_comparison.remote()
    with open("figures/zscore_window_test.png", "wb") as f:
        f.write(zimage_bytes)
    print("=== Z-score window-size test ===")
    print(zsummary)
    print("Saved figures/zscore_window_test.png")
