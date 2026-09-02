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


@app.function(image=image, timeout=1800)
def run_insulation():
    import io

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

    def insulation_score(matrix, window):
        n = matrix.shape[0]
        scores = np.full(n, np.nan)
        for i in range(window, n - window):
            block = matrix[i - window:i, i:i + window]
            scores[i] = block.mean()
        return scores

    def call_boundaries(scores, percentile=50):
        # percentile-based, not a fixed absolute number -- real
        # prominence scale depends heavily on the data and window size
        # (biological insulation dips are far subtler than the
        # deliberately dramatic synthetic test), so a fixed threshold
        # tuned on one dataset does not transfer to another
        log_scores = np.log2(scores + 1e-9)
        valid = ~np.isnan(log_scores)
        idx = np.where(valid)[0]
        vals = log_scores[idx]
        all_boundaries = []
        for k in range(1, len(vals) - 1):
            if vals[k] < vals[k - 1] and vals[k] < vals[k + 1]:
                prominence = min(vals[k - 1], vals[k + 1]) - vals[k]
                all_boundaries.append((idx[k], prominence))
        if not all_boundaries:
            return []
        cutoff = np.percentile([p for _, p in all_boundaries], percentile)
        return [(b, p) for b, p in all_boundaries if p >= cutoff]

    hic = load_hic()
    chrom, resolution = "21", 10000
    start, end = 20_000_000, 22_000_000
    matrix = get_matrix(hic, chrom, start, end, resolution, normalization="KR")

    window_bins = 8
    scores = insulation_score(matrix, window_bins)
    boundaries = call_boundaries(scores, percentile=50)  # keep top half by prominence
    boundaries.sort(key=lambda x: -x[1])

    summary_lines = [f"Found {len(boundaries)} candidate boundaries (window={window_bins} bins)"]
    summary_lines.append("Top 10 by prominence:")
    for bin_idx, prom in boundaries[:10]:
        genomic_pos = start + bin_idx * resolution
        summary_lines.append(f"  bin={bin_idx}  position={genomic_pos:,} bp  prominence={prom:.3f}")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    ax1.imshow(np.log1p(matrix), cmap="Reds", origin="upper")
    for bin_idx, prom in boundaries[:15]:
        ax1.axvline(bin_idx, color="blue", linewidth=0.8, alpha=0.6)
        ax1.axhline(bin_idx, color="blue", linewidth=0.8, alpha=0.6)
    ax1.set_title(f"chr{chrom}:{start}-{end}, GM12878, {resolution}bp -- top insulation boundaries")

    ax2.plot(scores)
    for bin_idx, prom in boundaries[:15]:
        ax2.axvline(bin_idx, color="blue", linewidth=0.8, alpha=0.6)
    ax2.set_ylabel("Insulation score")
    ax2.set_xlabel("Bin")

    plt.tight_layout()
    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=150)
    plt.close()
    buf.seek(0)

    return buf.getvalue(), "\n".join(summary_lines)


@app.function(image=image, timeout=600)
def diagnose_insulation():
    import hicstraw
    import numpy as np

    def load_hic(url="https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined_30.hic"):
        return hicstraw.HiCFile(url)

    def get_matrix(hic, chrom, start, end, resolution, normalization="KR", data_type="observed"):
        mzd = hic.getMatrixZoomData(chrom, chrom, data_type, normalization, "BP", resolution)
        return mzd.getRecordsAsMatrix(start, end, start, end)

    def insulation_score(matrix, window):
        n = matrix.shape[0]
        scores = np.full(n, np.nan)
        for i in range(window, n - window):
            block = matrix[i - window:i, i:i + window]
            scores[i] = block.mean()
        return scores

    hic = load_hic()
    chrom, resolution = "21", 10000
    start, end = 20_000_000, 22_000_000
    matrix = get_matrix(hic, chrom, start, end, resolution, normalization="KR")

    lines = [f"matrix shape: {matrix.shape}"]
    lines.append(f"matrix has NaN: {np.isnan(matrix).any()}  count: {np.isnan(matrix).sum()}")
    lines.append(f"matrix min/max/mean: {np.nanmin(matrix):.3f} / {np.nanmax(matrix):.3f} / {np.nanmean(matrix):.3f}")

    for window in [4, 8, 15]:
        scores = insulation_score(matrix, window)
        valid = scores[~np.isnan(scores)]
        log_scores = np.log2(valid + 1e-9)

        # count raw local minima with NO threshold at all
        raw_minima = 0
        for k in range(1, len(log_scores) - 1):
            if log_scores[k] < log_scores[k-1] and log_scores[k] < log_scores[k+1]:
                raw_minima += 1

        # actual prominence values for those raw minima
        prominences = []
        for k in range(1, len(log_scores) - 1):
            if log_scores[k] < log_scores[k-1] and log_scores[k] < log_scores[k+1]:
                prominences.append(min(log_scores[k-1], log_scores[k+1]) - log_scores[k])

        lines.append(f"\nwindow={window}: valid_bins={len(valid)}  raw_local_minima={raw_minima}")
        if prominences:
            lines.append(f"  prominence distribution: min={min(prominences):.4f} max={max(prominences):.4f} "
                          f"median={np.median(prominences):.4f} mean={np.mean(prominences):.4f}")
        lines.append(f"  log2 score range: min={log_scores.min():.3f} max={log_scores.max():.3f} std={log_scores.std():.3f}")

    return "\n".join(lines)


@app.function(image=image, timeout=1800)
def run_integration():
    import io

    import gudhi
    import hicstraw
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.stats import spearmanr

    def load_hic(url="https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined_30.hic"):
        return hicstraw.HiCFile(url)

    def get_matrix(hic, chrom, start, end, resolution, normalization="KR", data_type="observed"):
        mzd = hic.getMatrixZoomData(chrom, chrom, data_type, normalization, "BP", resolution)
        return mzd.getRecordsAsMatrix(start, end, start, end)

    def insulation_score(matrix, window):
        n = matrix.shape[0]
        scores = np.full(n, np.nan)
        for i in range(window, n - window):
            block = matrix[i - window:i, i:i + window]
            scores[i] = block.mean()
        return scores

    def call_boundaries(scores, percentile=50):
        log_scores = np.log2(scores + 1e-9)
        valid = ~np.isnan(log_scores)
        idx = np.where(valid)[0]
        vals = log_scores[idx]
        all_b = []
        for k in range(1, len(vals) - 1):
            if vals[k] < vals[k - 1] and vals[k] < vals[k + 1]:
                prom = min(vals[k - 1], vals[k + 1]) - vals[k]
                all_b.append((idx[k], prom))
        if not all_b:
            return []
        cutoff = np.percentile([p for _, p in all_b], percentile)
        return [(b, p) for b, p in all_b if p >= cutoff]

    def mask_diagonal(matrix):
        masked = matrix.copy()
        off = matrix[~np.eye(matrix.shape[0], dtype=bool)]
        bg = np.median(off[off > 0]) if np.any(off > 0) else 0
        np.fill_diagonal(masked, bg)
        return masked

    def domain_persistence(matrix, b_start, b_end, margin=3):
        n = matrix.shape[0]
        lo = max(0, b_start - margin)
        hi = min(n, b_end + margin)
        sub = matrix[lo:hi, lo:hi]
        sub = mask_diagonal(sub)
        filt = -sub
        cc = gudhi.CubicalComplex(top_dimensional_cells=filt)
        cc.compute_persistence()
        h0 = cc.persistence_intervals_in_dimension(0)
        finite = h0[np.isfinite(h0[:, 1])]
        if len(finite) == 0:
            return 0.0
        pers = finite[:, 1] - finite[:, 0]
        return pers.max()

    hic = load_hic()
    chrom, resolution = "21", 10000
    start, end = 20_000_000, 22_000_000

    # Stage 1 uses raw KR-normalized observed counts, matching the
    # original design (insulation.py) -- unchanged from current field practice
    kr_matrix = get_matrix(hic, chrom, start, end, resolution, normalization="KR")
    # Stage 2 uses O/E, per Findings 1 and 2 -- raw counts are dominated
    # by distance decay, which O/E corrects for
    oe_matrix = get_matrix(hic, chrom, start, end, resolution, normalization="KR", data_type="oe")

    window_bins = 8
    scores = insulation_score(kr_matrix, window_bins)
    boundaries = call_boundaries(scores, percentile=50)
    boundaries.sort(key=lambda x: x[0])

    positions = sorted(set([0] + [b for b, p in boundaries] + [kr_matrix.shape[0]]))

    results = []
    for i in range(len(positions) - 1):
        b_start, b_end = positions[i], positions[i + 1]
        if b_end - b_start < 2:
            continue
        pers = domain_persistence(oe_matrix, b_start, b_end)
        prom_candidates = [p for b, p in boundaries if b in (b_start, b_end)]
        weak_prom = min(prom_candidates) if prom_candidates else float("nan")
        results.append((b_start, b_end, weak_prom, pers))

    # domains under 5 bins barely have room for a meaningful persistence
    # computation and are likely just noise -- exclude them rather than
    # let them dilute the comparison between the two methods
    min_domain_size = 5
    n_before = len(results)
    results = [r for r in results if (r[1] - r[0]) >= min_domain_size]
    n_dropped = n_before - len(results)

    summary_lines = [f"Stage 1 found {len(boundaries)} boundaries -> {n_before} candidate domains "
                      f"({n_dropped} dropped for being under {min_domain_size} bins, {len(results)} kept)\n"]
    summary_lines.append(f"{'domain':>15} {'stage1_prominence':>18} {'stage2_persistence':>18}")
    for b_start, b_end, prom, pers in results:
        prom_str = f"{prom:.4f}" if not np.isnan(prom) else "n/a"
        summary_lines.append(f"[{b_start:>4},{b_end:>4}] {prom_str:>18} {pers:>18.3f}")

    proms = [r[2] for r in results if not np.isnan(r[2])]
    perss = [r[3] for r in results if not np.isnan(r[2])]
    if len(proms) > 2:
        rho, pval = spearmanr(proms, perss)
        summary_lines.append(f"\nSpearman correlation (Stage1 prominence vs Stage2 persistence): "
                              f"rho={rho:.3f}  p={pval:.3f}")
        summary_lines.append("(This is Claim 1's actual question, on one region: do the two methods")
        summary_lines.append("rank boundaries the same way? Low/insignificant correlation would mean")
        summary_lines.append("persistence is capturing something insulation score does not -- which")
        summary_lines.append("could support persistence as a genuinely different signal, not just a")
        summary_lines.append("noisier version of the same one. This is one region, not the full")
        summary_lines.append("benchmark -- a real result needs many regions and the CTCF/cohesin")
        summary_lines.append("proxy metrics from the paper's Section 4.5.")

    fig, ax = plt.subplots(figsize=(7, 6))
    if len(proms) > 0:
        ax.scatter(proms, perss)
        ax.set_xlabel("Stage 1: insulation-score prominence")
        ax.set_ylabel("Stage 2: persistence")
        ax.set_title("Do the two methods rank domains the same way?")
    plt.tight_layout()
    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=150)
    plt.close()
    buf.seek(0)

    return buf.getvalue(), "\n".join(summary_lines)


@app.local_entrypoint()
def main(task: str = "insulation"):
    if task == "integrate":
        image_bytes, summary = run_integration.remote()
        with open("figures/stage1_stage2_integration.png", "wb") as f:
            f.write(image_bytes)
        print(summary)
        print("Saved figures/stage1_stage2_integration.png")

    elif task == "diagnose":
        print(diagnose_insulation.remote())

    elif task == "insulation":
        image_bytes, summary = run_insulation.remote()
        with open("figures/insulation_boundaries.png", "wb") as f:
            f.write(image_bytes)
        print(summary)
        print("Saved figures/insulation_boundaries.png")

    elif task == "sweep":
        image_bytes, summary = run_sweep.remote()
        with open("figures/distance_band_sweep_extended.png", "wb") as f:
            f.write(image_bytes)
        print(summary)
        print("Saved figures/distance_band_sweep_extended.png")

    elif task == "zscore":
        image_bytes, summary = run_zscore_comparison.remote()
        with open("figures/zscore_window_test.png", "wb") as f:
            f.write(image_bytes)
        print(summary)
        print("Saved figures/zscore_window_test.png")

    elif task == "all":
        for t in ["insulation", "sweep", "zscore"]:
            main(task=t)

    else:
        print(f"Unknown task '{task}'. Choose from: insulation, sweep, zscore, all")
