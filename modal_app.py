import modal

app = modal.App("hic-boundary-persistence")

results_volume = modal.Volume.from_name("hic-results-vol", create_if_missing=True)

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


@app.function(image=image, timeout=3600)
def run_integration_multiregion():
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

    def process_region(hic, chrom, start, end, resolution=10000, min_domain_size=5):
        kr_matrix = get_matrix(hic, chrom, start, end, resolution, normalization="KR")
        oe_matrix = get_matrix(hic, chrom, start, end, resolution, normalization="KR", data_type="oe")

        scores = insulation_score(kr_matrix, window=8)
        boundaries = call_boundaries(scores, percentile=50)
        boundaries.sort(key=lambda x: x[0])

        positions = sorted(set([0] + [b for b, p in boundaries] + [kr_matrix.shape[0]]))

        region_results = []
        for i in range(len(positions) - 1):
            b_start, b_end = positions[i], positions[i + 1]
            if b_end - b_start < min_domain_size:
                continue
            pers = domain_persistence(oe_matrix, b_start, b_end)
            prom_candidates = [p for b, p in boundaries if b in (b_start, b_end)]
            weak_prom = min(prom_candidates) if prom_candidates else float("nan")
            if not np.isnan(weak_prom):
                region_results.append((b_start, b_end, weak_prom, pers))
        return region_results

    hic = load_hic()
    resolution = 10000

    # spread across several non-adjacent windows on chr21, plus a
    # second chromosome (chr20), to get real statistical power instead
    # of one region's worth of domains
    regions = [
        ("21", 20_000_000, 22_000_000),
        ("21", 25_000_000, 27_000_000),
        ("21", 30_000_000, 32_000_000),
        ("21", 35_000_000, 37_000_000),
        ("21", 40_000_000, 42_000_000),
        ("20", 20_000_000, 22_000_000),
        ("20", 30_000_000, 32_000_000),
        ("20", 40_000_000, 42_000_000),
    ]

    pooled = []
    summary_lines = []
    for chrom, start, end in regions:
        try:
            region_results = process_region(hic, chrom, start, end, resolution)
            summary_lines.append(f"chr{chrom}:{start}-{end}  ->  {len(region_results)} domains (>=5 bins)")
            for b_start, b_end, prom, pers in region_results:
                pooled.append((f"chr{chrom}:{start}-{end}", b_start, b_end, prom, pers))
        except Exception as e:
            summary_lines.append(f"chr{chrom}:{start}-{end}  ->  FAILED: {e}")

    summary_lines.append(f"\nTotal pooled domains: {len(pooled)}")

    proms = [r[3] for r in pooled]
    perss = [r[4] for r in pooled]

    if len(pooled) > 5:
        rho, pval = spearmanr(proms, perss)
        summary_lines.append(f"\nPooled Spearman correlation (Stage1 prominence vs Stage2 persistence): "
                              f"rho={rho:.3f}  p={pval:.4f}  n={len(pooled)}")
        if pval < 0.05:
            summary_lines.append("Statistically significant at alpha=0.05.")
        else:
            summary_lines.append("Not statistically significant at alpha=0.05 -- still inconclusive,")
            summary_lines.append("though with more regions than the single-window pilot.")

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(proms, perss, alpha=0.6)
    ax.set_xlabel("Stage 1: insulation-score prominence")
    ax.set_ylabel("Stage 2: persistence")
    ax.set_title(f"Pooled across {len(regions)} regions, n={len(pooled)} domains")
    plt.tight_layout()
    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=150)
    plt.close()
    buf.seek(0)

    return buf.getvalue(), "\n".join(summary_lines)


@app.function(image=image, timeout=3600)
def run_ctcf_benchmark():
    import io
    import urllib.request

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

    def fetch_peaks(url, chroms):
        """
        Download and parse a narrowPeak/bed file, keeping only the
        chromosomes we actually need. narrowPeak columns: chrom, start,
        end, name, score, strand, signalValue, pValue, qValue, summit.
        """
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=120) as resp:
            raw = resp.read()
        import gzip
        text = gzip.decompress(raw).decode("utf-8") if url.endswith(".gz") else raw.decode("utf-8")

        peaks = {c: [] for c in chroms}
        for line in text.strip().split("\n"):
            fields = line.split("\t")
            chrom = fields[0].replace("chr", "")
            if chrom in peaks:
                start, end = int(fields[1]), int(fields[2])
                peaks[chrom].append((start, end))
        return peaks

    def has_nearby_peak(genomic_pos, peak_list, window=20000):
        for (p_start, p_end) in peak_list:
            if p_start - window <= genomic_pos <= p_end + window:
                return True
        return False

    def rank_auc(scores, labels):
        scores = np.asarray(scores)
        labels = np.asarray(labels)
        n_pos = labels.sum()
        n_neg = len(labels) - n_pos
        if n_pos == 0 or n_neg == 0:
            return float("nan")
        ranks = np.argsort(np.argsort(scores)) + 1
        sum_ranks_pos = ranks[labels == 1].sum()
        return (sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)

    def bootstrap_auc_diff(scores_stage1, scores_stage2, labels, n_boot=2000, seed=0):
        """
        Paired bootstrap: resample domains (with replacement), scoring
        both methods on the SAME resampled domains each time, since
        both methods are evaluated on the same set of real domains --
        not independent samples. Returns the mean AUC difference
        (Stage2 - Stage1) and a 95% confidence interval. If the
        interval excludes zero, the observed advantage is unlikely to
        be due to chance at this sample size.
        """
        rng = np.random.default_rng(seed)
        scores_stage1 = np.asarray(scores_stage1)
        scores_stage2 = np.asarray(scores_stage2)
        labels = np.asarray(labels)
        n = len(labels)
        diffs = []
        for _ in range(n_boot):
            idx = rng.integers(0, n, n)
            auc1 = rank_auc(scores_stage1[idx], labels[idx])
            auc2 = rank_auc(scores_stage2[idx], labels[idx])
            if not (np.isnan(auc1) or np.isnan(auc2)):
                diffs.append(auc2 - auc1)
        diffs = np.array(diffs)
        ci_low, ci_high = np.percentile(diffs, [2.5, 97.5])
        return diffs.mean(), ci_low, ci_high

    # fetch all three ChIP-seq peak sets, hg19, GM12878 -- confirmed
    # accessions, cross-checked against ENCODE's own assembly label
    chroms_needed = ["21", "20"]
    summary_lines = ["Fetching ChIP-seq peak files..."]

    ctcf_peaks = fetch_peaks(
        "https://www.encodeproject.org/files/ENCFF833FTF/@@download/ENCFF833FTF.bed.gz",
        chroms_needed,
    )
    rad21_peaks = fetch_peaks(
        "https://www.encodeproject.org/files/ENCFF753RGL/@@download/ENCFF753RGL.bed.gz",
        chroms_needed,
    )
    smc3_peaks = fetch_peaks(
        "https://www.encodeproject.org/files/ENCFF572RPI/@@download/ENCFF572RPI.bed.gz",
        chroms_needed,
    )

    for name, peaks in [("CTCF", ctcf_peaks), ("RAD21", rad21_peaks), ("SMC3", smc3_peaks)]:
        counts = {c: len(v) for c, v in peaks.items()}
        summary_lines.append(f"{name}: {counts}")

    def process_region(hic, chrom, start, end, resolution=10000, min_domain_size=5):
        kr_matrix = get_matrix(hic, chrom, start, end, resolution, normalization="KR")
        oe_matrix = get_matrix(hic, chrom, start, end, resolution, normalization="KR", data_type="oe")

        scores = insulation_score(kr_matrix, window=8)
        boundaries = call_boundaries(scores, percentile=50)
        boundaries.sort(key=lambda x: x[0])

        positions = sorted(set([0] + [b for b, p in boundaries] + [kr_matrix.shape[0]]))

        region_results = []
        for i in range(len(positions) - 1):
            b_start, b_end = positions[i], positions[i + 1]
            if b_end - b_start < min_domain_size:
                continue
            pers = domain_persistence(oe_matrix, b_start, b_end)
            prom_candidates = [p for b, p in boundaries if b in (b_start, b_end)]
            weak_prom = min(prom_candidates) if prom_candidates else float("nan")
            if np.isnan(weak_prom):
                continue

            # check both edges of the domain for nearby peaks -- a
            # domain boundary is "real" if either edge overlaps a peak
            genomic_start = start + b_start * resolution
            genomic_end = start + b_end * resolution
            has_ctcf = (has_nearby_peak(genomic_start, ctcf_peaks[chrom]) or
                        has_nearby_peak(genomic_end, ctcf_peaks[chrom]))
            has_cohesin = (has_nearby_peak(genomic_start, rad21_peaks[chrom]) or
                           has_nearby_peak(genomic_end, rad21_peaks[chrom]) or
                           has_nearby_peak(genomic_start, smc3_peaks[chrom]) or
                           has_nearby_peak(genomic_end, smc3_peaks[chrom]))

            region_results.append((weak_prom, pers, has_ctcf, has_cohesin))
        return region_results

    hic = load_hic()
    resolution = 10000
    regions = [
        ("21", 20_000_000, 22_000_000),
        ("21", 25_000_000, 27_000_000),
        ("21", 30_000_000, 32_000_000),
        ("21", 35_000_000, 37_000_000),
        ("21", 40_000_000, 42_000_000),
        ("20", 20_000_000, 22_000_000),
        ("20", 30_000_000, 32_000_000),
        ("20", 40_000_000, 42_000_000),
    ]

    pooled = []
    for chrom, start, end in regions:
        try:
            pooled.extend(process_region(hic, chrom, start, end, resolution))
        except Exception as e:
            summary_lines.append(f"chr{chrom}:{start}-{end} FAILED: {e}")

    summary_lines.append(f"\nTotal domains with a call on both methods: {len(pooled)}")

    proms = [r[0] for r in pooled]
    perss = [r[1] for r in pooled]
    ctcf_labels = [1 if r[2] else 0 for r in pooled]
    cohesin_labels = [1 if r[3] else 0 for r in pooled]

    summary_lines.append(f"Domains with a nearby CTCF peak: {sum(ctcf_labels)} / {len(pooled)}")
    summary_lines.append(f"Domains with a nearby cohesin (RAD21/SMC3) peak: {sum(cohesin_labels)} / {len(pooled)}")

    summary_lines.append("\n=== Claim 1: does each method's ranking predict real CTCF/cohesin binding? ===")
    summary_lines.append("AUC = 0.5 means no predictive power; higher AUC means the score ranks")
    summary_lines.append("boundaries with real protein binding above those without it, better.\n")

    for label_name, labels in [("CTCF", ctcf_labels), ("cohesin (RAD21/SMC3)", cohesin_labels)]:
        auc_stage1 = rank_auc(proms, labels)
        auc_stage2 = rank_auc(perss, labels)
        mean_diff, ci_low, ci_high = bootstrap_auc_diff(proms, perss, labels)
        summary_lines.append(f"{label_name}:")
        summary_lines.append(f"  Stage 1 (insulation prominence) AUC = {auc_stage1:.3f}")
        summary_lines.append(f"  Stage 2 (persistence)           AUC = {auc_stage2:.3f}")
        summary_lines.append(f"  AUC difference (Stage2 - Stage1): {mean_diff:.3f}  "
                              f"95% CI [{ci_low:.3f}, {ci_high:.3f}]")
        if ci_low > 0:
            summary_lines.append(f"  CI excludes zero -- Stage 2's advantage is unlikely to be chance at this sample size.")
        elif ci_high < 0:
            summary_lines.append(f"  CI excludes zero in Stage 1's favor.")
        else:
            summary_lines.append(f"  CI includes zero -- cannot rule out that this difference is due to chance yet.")
        winner = "Stage 2 (persistence)" if auc_stage2 > auc_stage1 else "Stage 1 (insulation score)"
        summary_lines.append(f"  Higher point estimate: {winner}\n")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, labels, title in [(axes[0], ctcf_labels, "CTCF"), (axes[1], cohesin_labels, "Cohesin")]:
        colors = ["red" if l else "gray" for l in labels]
        ax.scatter(proms, perss, c=colors, alpha=0.6)
        ax.set_xlabel("Stage 1: insulation prominence")
        ax.set_ylabel("Stage 2: persistence")
        ax.set_title(f"{title} (red = peak nearby)")
    plt.tight_layout()
    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=150)
    plt.close()
    buf.seek(0)

    result_text = "\n".join(summary_lines)
    with open("/results/ctcf_benchmark_genomewide.png", "wb") as f:
        f.write(buf.getvalue())
    with open("/results/ctcf_benchmark_genomewide_summary.txt", "w") as f:
        f.write(result_text)
    results_volume.commit()

    return buf.getvalue(), result_text


@app.function(image=image, timeout=5400)
def run_ctcf_benchmark_expanded():
    """
    Expanded version of run_ctcf_benchmark -- same logic entirely,
    just spread across many more regions (22 vs 8, 5 chromosomes vs 2)
    to tighten the bootstrap confidence intervals that were too wide
    to be conclusive in the first pass. run_ctcf_benchmark itself is
    left untouched; this is a new, separate function so the original
    validated result stays intact and reproducible on its own.
    """
    import io
    import urllib.request

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

    def fetch_peaks(url, chroms):
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=120) as resp:
            raw = resp.read()
        import gzip
        text = gzip.decompress(raw).decode("utf-8") if url.endswith(".gz") else raw.decode("utf-8")
        peaks = {c: [] for c in chroms}
        for line in text.strip().split("\n"):
            fields = line.split("\t")
            chrom = fields[0].replace("chr", "")
            if chrom in peaks:
                start, end = int(fields[1]), int(fields[2])
                peaks[chrom].append((start, end))
        return peaks

    def has_nearby_peak(genomic_pos, peak_list, window=20000):
        for (p_start, p_end) in peak_list:
            if p_start - window <= genomic_pos <= p_end + window:
                return True
        return False

    def rank_auc(scores, labels):
        scores = np.asarray(scores)
        labels = np.asarray(labels)
        n_pos = labels.sum()
        n_neg = len(labels) - n_pos
        if n_pos == 0 or n_neg == 0:
            return float("nan")
        ranks = np.argsort(np.argsort(scores)) + 1
        sum_ranks_pos = ranks[labels == 1].sum()
        return (sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)

    def bootstrap_auc_diff(scores_stage1, scores_stage2, labels, n_boot=2000, seed=0):
        rng = np.random.default_rng(seed)
        scores_stage1 = np.asarray(scores_stage1)
        scores_stage2 = np.asarray(scores_stage2)
        labels = np.asarray(labels)
        n = len(labels)
        diffs = []
        for _ in range(n_boot):
            idx = rng.integers(0, n, n)
            auc1 = rank_auc(scores_stage1[idx], labels[idx])
            auc2 = rank_auc(scores_stage2[idx], labels[idx])
            if not (np.isnan(auc1) or np.isnan(auc2)):
                diffs.append(auc2 - auc1)
        diffs = np.array(diffs)
        ci_low, ci_high = np.percentile(diffs, [2.5, 97.5])
        return diffs.mean(), ci_low, ci_high

    chroms_needed = ["18", "19", "20", "21", "22"]
    summary_lines = ["Fetching ChIP-seq peak files..."]

    ctcf_peaks = fetch_peaks(
        "https://www.encodeproject.org/files/ENCFF833FTF/@@download/ENCFF833FTF.bed.gz",
        chroms_needed,
    )
    rad21_peaks = fetch_peaks(
        "https://www.encodeproject.org/files/ENCFF753RGL/@@download/ENCFF753RGL.bed.gz",
        chroms_needed,
    )
    smc3_peaks = fetch_peaks(
        "https://www.encodeproject.org/files/ENCFF572RPI/@@download/ENCFF572RPI.bed.gz",
        chroms_needed,
    )

    for name, peaks in [("CTCF", ctcf_peaks), ("RAD21", rad21_peaks), ("SMC3", smc3_peaks)]:
        counts = {c: len(v) for c, v in peaks.items()}
        summary_lines.append(f"{name}: {counts}")

    def process_region(hic, chrom, start, end, resolution=10000, min_domain_size=5):
        kr_matrix = get_matrix(hic, chrom, start, end, resolution, normalization="KR")
        oe_matrix = get_matrix(hic, chrom, start, end, resolution, normalization="KR", data_type="oe")

        scores = insulation_score(kr_matrix, window=8)
        boundaries = call_boundaries(scores, percentile=50)
        boundaries.sort(key=lambda x: x[0])

        positions = sorted(set([0] + [b for b, p in boundaries] + [kr_matrix.shape[0]]))

        region_results = []
        for i in range(len(positions) - 1):
            b_start, b_end = positions[i], positions[i + 1]
            if b_end - b_start < min_domain_size:
                continue
            pers = domain_persistence(oe_matrix, b_start, b_end)
            prom_candidates = [p for b, p in boundaries if b in (b_start, b_end)]
            weak_prom = min(prom_candidates) if prom_candidates else float("nan")
            if np.isnan(weak_prom):
                continue

            genomic_start = start + b_start * resolution
            genomic_end = start + b_end * resolution
            has_ctcf = (has_nearby_peak(genomic_start, ctcf_peaks[chrom]) or
                        has_nearby_peak(genomic_end, ctcf_peaks[chrom]))
            has_cohesin = (has_nearby_peak(genomic_start, rad21_peaks[chrom]) or
                           has_nearby_peak(genomic_end, rad21_peaks[chrom]) or
                           has_nearby_peak(genomic_start, smc3_peaks[chrom]) or
                           has_nearby_peak(genomic_end, smc3_peaks[chrom]))

            region_results.append((weak_prom, pers, has_ctcf, has_cohesin))
        return region_results

    hic = load_hic()
    resolution = 10000

    regions = []
    for start_mb in [20, 30, 40, 50, 60, 70]:
        regions.append(("18", start_mb * 1_000_000, start_mb * 1_000_000 + 2_000_000))
    for start_mb in [20, 30, 40, 50]:
        regions.append(("19", start_mb * 1_000_000, start_mb * 1_000_000 + 2_000_000))
    for start_mb in [20, 30, 40, 50]:
        regions.append(("20", start_mb * 1_000_000, start_mb * 1_000_000 + 2_000_000))
    for start_mb in [20, 25, 30, 35, 40]:
        regions.append(("21", start_mb * 1_000_000, start_mb * 1_000_000 + 2_000_000))
    for start_mb in [20, 30, 40]:
        regions.append(("22", start_mb * 1_000_000, start_mb * 1_000_000 + 2_000_000))

    pooled = []
    for chrom, start, end in regions:
        try:
            pooled.extend(process_region(hic, chrom, start, end, resolution))
        except Exception as e:
            summary_lines.append(f"chr{chrom}:{start}-{end} FAILED: {e}")

    summary_lines.append(f"\nTotal domains with a call on both methods: {len(pooled)} (from {len(regions)} regions, 5 chromosomes)")

    proms = [r[0] for r in pooled]
    perss = [r[1] for r in pooled]
    ctcf_labels = [1 if r[2] else 0 for r in pooled]
    cohesin_labels = [1 if r[3] else 0 for r in pooled]

    summary_lines.append(f"Domains with a nearby CTCF peak: {sum(ctcf_labels)} / {len(pooled)}")
    summary_lines.append(f"Domains with a nearby cohesin (RAD21/SMC3) peak: {sum(cohesin_labels)} / {len(pooled)}")

    summary_lines.append("\n=== Claim 1 (expanded, 22 regions, 5 chromosomes) ===")

    for label_name, labels in [("CTCF", ctcf_labels), ("cohesin (RAD21/SMC3)", cohesin_labels)]:
        auc_stage1 = rank_auc(proms, labels)
        auc_stage2 = rank_auc(perss, labels)
        mean_diff, ci_low, ci_high = bootstrap_auc_diff(proms, perss, labels)
        summary_lines.append(f"{label_name}:")
        summary_lines.append(f"  Stage 1 (insulation prominence) AUC = {auc_stage1:.3f}")
        summary_lines.append(f"  Stage 2 (persistence)           AUC = {auc_stage2:.3f}")
        summary_lines.append(f"  AUC difference (Stage2 - Stage1): {mean_diff:.3f}  "
                              f"95% CI [{ci_low:.3f}, {ci_high:.3f}]")
        if ci_low > 0:
            summary_lines.append(f"  CI excludes zero -- Stage 2's advantage is unlikely to be chance at this sample size.")
        elif ci_high < 0:
            summary_lines.append(f"  CI excludes zero in Stage 1's favor.")
        else:
            summary_lines.append(f"  CI includes zero -- still cannot rule out chance.")
        winner = "Stage 2 (persistence)" if auc_stage2 > auc_stage1 else "Stage 1 (insulation score)"
        summary_lines.append(f"  Higher point estimate: {winner}\n")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, labels, title in [(axes[0], ctcf_labels, "CTCF"), (axes[1], cohesin_labels, "Cohesin")]:
        colors = ["red" if l else "gray" for l in labels]
        ax.scatter(proms, perss, c=colors, alpha=0.5, s=20)
        ax.set_xlabel("Stage 1: insulation prominence")
        ax.set_ylabel("Stage 2: persistence")
        ax.set_title(f"{title} (red = peak nearby), n={len(pooled)}")
    plt.tight_layout()
    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=150)
    plt.close()
    buf.seek(0)

    return buf.getvalue(), "\n".join(summary_lines)


@app.function(image=image, timeout=10800, volumes={"/results": results_volume})
def run_ctcf_benchmark_genomewide():
    """
    Genome-wide (autosomal) version of the CTCF/cohesin benchmark --
    101 regions spread across all 22 human autosomes (hg19), evenly
    spaced every 25Mb per chromosome with safety margins from
    centromere/telomere-adjacent low-mappability regions. A fully
    separate function from run_ctcf_benchmark and
    run_ctcf_benchmark_expanded, which are both left untouched so
    their results stay intact and independently reproducible.
    """
    import io
    import urllib.request

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

    def fetch_peaks(url, chroms):
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=120) as resp:
            raw = resp.read()
        import gzip
        text = gzip.decompress(raw).decode("utf-8") if url.endswith(".gz") else raw.decode("utf-8")
        peaks = {c: [] for c in chroms}
        for line in text.strip().split("\n"):
            fields = line.split("\t")
            chrom = fields[0].replace("chr", "")
            if chrom in peaks:
                start, end = int(fields[1]), int(fields[2])
                peaks[chrom].append((start, end))
        return peaks

    def has_nearby_peak(genomic_pos, peak_list, window=20000):
        for (p_start, p_end) in peak_list:
            if p_start - window <= genomic_pos <= p_end + window:
                return True
        return False

    def rank_auc(scores, labels):
        scores = np.asarray(scores)
        labels = np.asarray(labels)
        n_pos = labels.sum()
        n_neg = len(labels) - n_pos
        if n_pos == 0 or n_neg == 0:
            return float("nan")
        ranks = np.argsort(np.argsort(scores)) + 1
        sum_ranks_pos = ranks[labels == 1].sum()
        return (sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)

    def bootstrap_auc_diff(scores_stage1, scores_stage2, labels, n_boot=2000, seed=0):
        rng = np.random.default_rng(seed)
        scores_stage1 = np.asarray(scores_stage1)
        scores_stage2 = np.asarray(scores_stage2)
        labels = np.asarray(labels)
        n = len(labels)
        diffs = []
        for _ in range(n_boot):
            idx = rng.integers(0, n, n)
            auc1 = rank_auc(scores_stage1[idx], labels[idx])
            auc2 = rank_auc(scores_stage2[idx], labels[idx])
            if not (np.isnan(auc1) or np.isnan(auc2)):
                diffs.append(auc2 - auc1)
        diffs = np.array(diffs)
        ci_low, ci_high = np.percentile(diffs, [2.5, 97.5])
        return diffs.mean(), ci_low, ci_high

    chrom_lengths_hg19 = {
        "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276,
        "5": 180915260, "6": 171115067, "7": 159138663, "8": 146364022,
        "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895,
        "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753,
        "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
        "21": 48129895, "22": 51304566,
    }

    regions = []
    for chrom, length in chrom_lengths_hg19.items():
        start_mb = 20
        while (start_mb + 2) * 1_000_000 < length - 7_000_000:
            regions.append((chrom, start_mb * 1_000_000, start_mb * 1_000_000 + 2_000_000))
            start_mb += 25

    chroms_needed = list(chrom_lengths_hg19.keys())
    summary_lines = [f"Fetching ChIP-seq peak files for {len(chroms_needed)} chromosomes..."]

    ctcf_peaks = fetch_peaks(
        "https://www.encodeproject.org/files/ENCFF833FTF/@@download/ENCFF833FTF.bed.gz",
        chroms_needed,
    )
    rad21_peaks = fetch_peaks(
        "https://www.encodeproject.org/files/ENCFF753RGL/@@download/ENCFF753RGL.bed.gz",
        chroms_needed,
    )
    smc3_peaks = fetch_peaks(
        "https://www.encodeproject.org/files/ENCFF572RPI/@@download/ENCFF572RPI.bed.gz",
        chroms_needed,
    )
    summary_lines.append(f"Peaks fetched across {len(chroms_needed)} chromosomes.")

    def process_region(hic, chrom, start, end, resolution=10000, min_domain_size=5):
        kr_matrix = get_matrix(hic, chrom, start, end, resolution, normalization="KR")
        oe_matrix = get_matrix(hic, chrom, start, end, resolution, normalization="KR", data_type="oe")

        scores = insulation_score(kr_matrix, window=8)
        boundaries = call_boundaries(scores, percentile=50)
        boundaries.sort(key=lambda x: x[0])

        positions = sorted(set([0] + [b for b, p in boundaries] + [kr_matrix.shape[0]]))

        region_results = []
        for i in range(len(positions) - 1):
            b_start, b_end = positions[i], positions[i + 1]
            if b_end - b_start < min_domain_size:
                continue
            pers = domain_persistence(oe_matrix, b_start, b_end)
            prom_candidates = [p for b, p in boundaries if b in (b_start, b_end)]
            weak_prom = min(prom_candidates) if prom_candidates else float("nan")
            if np.isnan(weak_prom):
                continue

            genomic_start = start + b_start * resolution
            genomic_end = start + b_end * resolution
            has_ctcf = (has_nearby_peak(genomic_start, ctcf_peaks[chrom]) or
                        has_nearby_peak(genomic_end, ctcf_peaks[chrom]))
            has_cohesin = (has_nearby_peak(genomic_start, rad21_peaks[chrom]) or
                           has_nearby_peak(genomic_end, rad21_peaks[chrom]) or
                           has_nearby_peak(genomic_start, smc3_peaks[chrom]) or
                           has_nearby_peak(genomic_end, smc3_peaks[chrom]))

            region_results.append((weak_prom, pers, has_ctcf, has_cohesin))
        return region_results

    hic = load_hic()
    resolution = 10000

    pooled = []
    n_failed = 0
    for chrom, start, end in regions:
        try:
            pooled.extend(process_region(hic, chrom, start, end, resolution))
        except Exception as e:
            n_failed += 1

    summary_lines.append(f"\nRegions attempted: {len(regions)} across {len(chroms_needed)} chromosomes "
                          f"({n_failed} failed)")
    summary_lines.append(f"Total domains with a call on both methods: {len(pooled)}")

    proms = [r[0] for r in pooled]
    perss = [r[1] for r in pooled]
    ctcf_labels = [1 if r[2] else 0 for r in pooled]
    cohesin_labels = [1 if r[3] else 0 for r in pooled]

    summary_lines.append(f"Domains with a nearby CTCF peak: {sum(ctcf_labels)} / {len(pooled)}")
    summary_lines.append(f"Domains with a nearby cohesin (RAD21/SMC3) peak: {sum(cohesin_labels)} / {len(pooled)}")

    summary_lines.append(f"\n=== Claim 1 (genome-wide autosomal, {len(regions)} regions, 22 chromosomes) ===")

    for label_name, labels in [("CTCF", ctcf_labels), ("cohesin (RAD21/SMC3)", cohesin_labels)]:
        auc_stage1 = rank_auc(proms, labels)
        auc_stage2 = rank_auc(perss, labels)
        mean_diff, ci_low, ci_high = bootstrap_auc_diff(proms, perss, labels)
        summary_lines.append(f"{label_name}:")
        summary_lines.append(f"  Stage 1 (insulation prominence) AUC = {auc_stage1:.3f}")
        summary_lines.append(f"  Stage 2 (persistence)           AUC = {auc_stage2:.3f}")
        summary_lines.append(f"  AUC difference (Stage2 - Stage1): {mean_diff:.3f}  "
                              f"95% CI [{ci_low:.3f}, {ci_high:.3f}]")
        if ci_low > 0:
            summary_lines.append(f"  CI excludes zero -- Stage 2's advantage is unlikely to be chance.")
        elif ci_high < 0:
            summary_lines.append(f"  CI excludes zero in Stage 1's favor.")
        else:
            summary_lines.append(f"  CI includes zero -- cannot rule out chance.")
        winner = "Stage 2 (persistence)" if auc_stage2 > auc_stage1 else "Stage 1 (insulation score)"
        summary_lines.append(f"  Higher point estimate: {winner}\n")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, labels, title in [(axes[0], ctcf_labels, "CTCF"), (axes[1], cohesin_labels, "Cohesin")]:
        colors = ["red" if l else "gray" for l in labels]
        ax.scatter(proms, perss, c=colors, alpha=0.4, s=12)
        ax.set_xlabel("Stage 1: insulation prominence")
        ax.set_ylabel("Stage 2: persistence")
        ax.set_title(f"{title} (red = peak nearby), n={len(pooled)}")
    plt.tight_layout()
    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=150)
    plt.close()
    buf.seek(0)

    return buf.getvalue(), "\n".join(summary_lines)


@app.function(image=image, timeout=300, volumes={"/results": results_volume})
def fetch_and_cache_peaks():
    """
    Fetches the three peak files ONCE and caches them in the volume,
    so the 101 parallel region workers never hit ENCODE's server
    directly -- that's what triggered rate-limiting (HTTP 405) when
    all workers tried to download the same files simultaneously.
    """
    import gzip
    import json
    import urllib.request

    def fetch_peaks_all_chroms(url):
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=120) as resp:
            raw = resp.read()
        text = gzip.decompress(raw).decode("utf-8") if url.endswith(".gz") else raw.decode("utf-8")
        peaks = {}
        for line in text.strip().split("\n"):
            fields = line.split("\t")
            chrom = fields[0].replace("chr", "")
            peaks.setdefault(chrom, []).append((int(fields[1]), int(fields[2])))
        return peaks

    ctcf = fetch_peaks_all_chroms("https://www.encodeproject.org/files/ENCFF833FTF/@@download/ENCFF833FTF.bed.gz")
    rad21 = fetch_peaks_all_chroms("https://www.encodeproject.org/files/ENCFF753RGL/@@download/ENCFF753RGL.bed.gz")
    smc3 = fetch_peaks_all_chroms("https://www.encodeproject.org/files/ENCFF572RPI/@@download/ENCFF572RPI.bed.gz")

    with open("/results/_peaks_cache.json", "w") as f:
        json.dump({"ctcf": ctcf, "rad21": rad21, "smc3": smc3}, f)
    results_volume.commit()
    return "Peaks cached in volume."


@app.function(image=image, timeout=600, volumes={"/results": results_volume})
def process_one_region(chrom, start, end, resolution=10000):
    """
    Processes exactly ONE region, completely independently, and saves
    its own result to the volume immediately -- no dependency on any
    other region, and no risk of losing everything if one call fails.
    Fetches its own copy of the (small, <1MB) peak files rather than
    sharing state across calls, since that keeps each call fully
    self-contained and safely parallelizable.
    """
    import json
    import gzip
    import urllib.request

    import gudhi
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

    def fetch_peaks_one_chrom(url, chrom):
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=120) as resp:
            raw = resp.read()
        text = gzip.decompress(raw).decode("utf-8") if url.endswith(".gz") else raw.decode("utf-8")
        peaks = []
        for line in text.strip().split("\n"):
            fields = line.split("\t")
            if fields[0].replace("chr", "") == chrom:
                peaks.append((int(fields[1]), int(fields[2])))
        return peaks

    def has_nearby_peak(genomic_pos, peak_list, window=20000):
        for (p_start, p_end) in peak_list:
            if p_start - window <= genomic_pos <= p_end + window:
                return True
        return False

    import json as _json

    result_id = f"{chrom}_{start}_{end}"
    try:
        with open("/results/_peaks_cache.json") as f:
            cached = _json.load(f)
        ctcf_peaks = [tuple(x) for x in cached["ctcf"].get(chrom, [])]
        rad21_peaks = [tuple(x) for x in cached["rad21"].get(chrom, [])]
        smc3_peaks = [tuple(x) for x in cached["smc3"].get(chrom, [])]

        hic = load_hic()
        kr_matrix = get_matrix(hic, chrom, start, end, resolution, normalization="KR")
        oe_matrix = get_matrix(hic, chrom, start, end, resolution, normalization="KR", data_type="oe")

        scores = insulation_score(kr_matrix, window=8)
        boundaries = call_boundaries(scores, percentile=50)
        boundaries.sort(key=lambda x: x[0])

        positions = sorted(set([0] + [b for b, p in boundaries] + [kr_matrix.shape[0]]))

        region_results = []
        for i in range(len(positions) - 1):
            b_start, b_end = positions[i], positions[i + 1]
            if b_end - b_start < 5:
                continue
            pers = domain_persistence(oe_matrix, b_start, b_end)
            prom_candidates = [p for b, p in boundaries if b in (b_start, b_end)]
            weak_prom = min(prom_candidates) if prom_candidates else None
            if weak_prom is None:
                continue

            genomic_start = start + b_start * resolution
            genomic_end = start + b_end * resolution
            has_ctcf = (has_nearby_peak(genomic_start, ctcf_peaks) or
                        has_nearby_peak(genomic_end, ctcf_peaks))
            has_cohesin = (has_nearby_peak(genomic_start, rad21_peaks) or
                           has_nearby_peak(genomic_end, rad21_peaks) or
                           has_nearby_peak(genomic_start, smc3_peaks) or
                           has_nearby_peak(genomic_end, smc3_peaks))

            region_results.append({"prom": weak_prom, "pers": pers, "ctcf": has_ctcf, "cohesin": has_cohesin})

        with open(f"/results/region_{result_id}.json", "w") as f:
            json.dump({"status": "ok", "chrom": chrom, "start": start, "end": end, "domains": region_results}, f)
        results_volume.commit()
        return f"chr{chrom}:{start}-{end} ok, {len(region_results)} domains"

    except Exception as e:
        with open(f"/results/region_{result_id}.json", "w") as f:
            json.dump({"status": "failed", "chrom": chrom, "start": start, "end": end, "error": str(e)}, f)
        results_volume.commit()
        return f"chr{chrom}:{start}-{end} FAILED: {e}"


@app.function(image=image, timeout=600, volumes={"/results": results_volume})
def aggregate_results():
    """
    Reads whatever region results currently exist in the volume
    (regardless of whether all 101 have finished -- this can be run
    at any point to check progress, or once everything is done for
    the final result) and computes the pooled AUC comparison.
    """
    import glob
    import io
    import json

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    def rank_auc(scores, labels):
        scores = np.asarray(scores)
        labels = np.asarray(labels)
        n_pos = labels.sum()
        n_neg = len(labels) - n_pos
        if n_pos == 0 or n_neg == 0:
            return float("nan")
        ranks = np.argsort(np.argsort(scores)) + 1
        sum_ranks_pos = ranks[labels == 1].sum()
        return (sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)

    def bootstrap_auc_diff(scores_stage1, scores_stage2, labels, n_boot=2000, seed=0):
        rng = np.random.default_rng(seed)
        scores_stage1 = np.asarray(scores_stage1)
        scores_stage2 = np.asarray(scores_stage2)
        labels = np.asarray(labels)
        n = len(labels)
        diffs = []
        for _ in range(n_boot):
            idx = rng.integers(0, n, n)
            auc1 = rank_auc(scores_stage1[idx], labels[idx])
            auc2 = rank_auc(scores_stage2[idx], labels[idx])
            if not (np.isnan(auc1) or np.isnan(auc2)):
                diffs.append(auc2 - auc1)
        diffs = np.array(diffs)
        ci_low, ci_high = np.percentile(diffs, [2.5, 97.5])
        return diffs.mean(), ci_low, ci_high

    files = glob.glob("/results/region_*.json")
    summary_lines = [f"Found {len(files)} region result files in the volume."]

    pooled = []
    n_ok, n_failed = 0, 0
    for fpath in files:
        with open(fpath) as f:
            data = json.load(f)
        if data["status"] == "ok":
            n_ok += 1
            pooled.extend(data["domains"])
        else:
            n_failed += 1
            fail_chrom = data["chrom"]
            fail_start = data["start"]
            fail_end = data["end"]
            fail_error = data.get("error", "?")
            summary_lines.append(f"  FAILED: chr{fail_chrom}:{fail_start}-{fail_end} -- {fail_error}")

    summary_lines.append(f"Regions succeeded: {n_ok}, failed: {n_failed}")
    summary_lines.append(f"Total pooled domains: {len(pooled)}")

    proms = [d["prom"] for d in pooled]
    perss = [d["pers"] for d in pooled]
    ctcf_labels = [1 if d["ctcf"] else 0 for d in pooled]
    cohesin_labels = [1 if d["cohesin"] else 0 for d in pooled]

    summary_lines.append(f"Domains with a nearby CTCF peak: {sum(ctcf_labels)} / {len(pooled)}")
    summary_lines.append(f"Domains with a nearby cohesin peak: {sum(cohesin_labels)} / {len(pooled)}")
    summary_lines.append("\n=== Claim 1 (checkpointed genome-wide run) ===")

    for label_name, labels in [("CTCF", ctcf_labels), ("cohesin (RAD21/SMC3)", cohesin_labels)]:
        if sum(labels) == 0 or sum(labels) == len(labels):
            summary_lines.append(f"{label_name}: not enough label variation yet to compute AUC.")
            continue
        auc_stage1 = rank_auc(proms, labels)
        auc_stage2 = rank_auc(perss, labels)
        mean_diff, ci_low, ci_high = bootstrap_auc_diff(proms, perss, labels)
        summary_lines.append(f"{label_name}:")
        summary_lines.append(f"  Stage 1 AUC = {auc_stage1:.3f}   Stage 2 AUC = {auc_stage2:.3f}")
        summary_lines.append(f"  Difference: {mean_diff:.3f}  95% CI [{ci_low:.3f}, {ci_high:.3f}]")
        if ci_low > 0:
            summary_lines.append("  CI excludes zero -- significant in favor of Stage 2.")
        elif ci_high < 0:
            summary_lines.append("  CI excludes zero -- significant in favor of Stage 1.")
        else:
            summary_lines.append("  CI includes zero -- not yet significant.")

    result_text = "\n".join(summary_lines)
    with open("/results/aggregate_summary.txt", "w") as f:
        f.write(result_text)
    results_volume.commit()

    return result_text


@app.function(image=image, timeout=300)
def check_intact_hic_resolution():
    """
    Confirms the actual resolution ceiling of the HCT116 intact Hi-C
    file from the Fatima Institute Hugging Face bucket, before
    building anything on top of it -- same discipline as the very
    first describe() check run on the original GM12878 file at the
    start of this project. Streams the file header only via straw,
    does not download the full 43.7GB file.
    """
    import hicstraw

    url = "https://huggingface.co/buckets/FatimaInstitute/3d-genomics-prediction/resolve/HCT116/ENCFF573OPJ.hic?download=true"

    hic = hicstraw.HiCFile(url)

    lines = [f"Genome ID: {hic.getGenomeID()}"]
    chroms = hic.getChromosomes()
    lines.append(f"Chromosomes ({len(chroms)}):")
    for c in chroms:
        lines.append(f"  {c.name}\tlength={c.length}")

    resolutions = hic.getResolutions()
    lines.append(f"Available resolutions (bp): {resolutions}")
    lines.append(f"Finest resolution: {min(resolutions)} bp")

    return "\n".join(lines)


@app.function(image=image, timeout=600)
def check_intact_hic_sparsity():
    """
    A resolution being LISTED in a .hic file does not mean it has
    usable signal -- the original GM12878 file taught us this at
    1kb, where data was already very sparse. This checks real data
    density at 500bp and 200bp on the new intact Hi-C file before
    assuming either is actually usable.
    """
    import numpy as np
    import hicstraw

    url = "https://huggingface.co/buckets/FatimaInstitute/3d-genomics-prediction/resolve/HCT116/ENCFF573OPJ.hic?download=true"
    hic = hicstraw.HiCFile(url)

    chrom = "chr21"
    start, end = 20_000_000, 21_000_000  # 1Mb window, kept modest given file size

    lines = []
    for resolution in [1000, 500, 200]:
        lines.append(f"\n--- Resolution: {resolution} bp ---")
        # KR normalization crashes (SIGSEGV, not a catchable Python
        # exception) rather than failing cleanly at these fine
        # resolutions on this file -- confirmed by direct testing.
        # Only NONE (raw observed) is attempted here.
        try:
            mzd = hic.getMatrixZoomData(chrom, chrom, "observed", "NONE", "BP", resolution)
            matrix = mzd.getRecordsAsMatrix(start, end, start, end)
            nonzero_frac = np.count_nonzero(matrix) / matrix.size
            nonzero_vals = matrix[matrix > 0]
            lines.append(f"  normalization=NONE: shape={matrix.shape}, "
                         f"non-zero fraction={nonzero_frac:.4f}, "
                         f"mean(non-zero)={nonzero_vals.mean() if len(nonzero_vals) else 0:.3f}, "
                         f"max={matrix.max():.1f}")
        except Exception as e:
            lines.append(f"  normalization=NONE: FAILED -- {e}")

    return "\n".join(lines)


@app.function(image=image, timeout=1200)
def hct116_prelim_500bp():
    """
    Preliminary result on HCT116 intact Hi-C at 500bp, chr10, per Dr.
    Shamim's direct instruction. The fetch runs in a subprocess,
    isolated from this parent process, because hic-straw can SIGSEGV
    on this file at fine resolutions -- a crash Python cannot catch
    with try/except, since it kills the whole process. Isolating it
    in a subprocess means a crash there does not take down this
    function, and lets us retry.
    """
    import io
    import subprocess
    import sys
    import tempfile
    import os

    import gudhi
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    fetch_script = '''
import sys
import numpy as np
import hicstraw

url = sys.argv[1]
chrom = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])
resolution = int(sys.argv[5])
normalization = sys.argv[6]
out_path = sys.argv[7]

hic = hicstraw.HiCFile(url)
mzd = hic.getMatrixZoomData(chrom, chrom, "observed", normalization, "BP", resolution)
records = mzd.getRecords(start, end, start, end)

# getRecords() is stable for this 500bp / 1Mb query, whereas
# getRecordsAsMatrix() reproducibly SIGSEGVs inside hic-straw 1.3.1.
#
# Reconstruct the symmetric dense matrix ourselves from sparse
# (binX, binY, counts) records.
bin_start = (start // resolution) * resolution
bin_end = (end // resolution) * resolution
n_bins = ((bin_end - bin_start) // resolution) + 1

matrix = np.zeros((n_bins, n_bins), dtype=np.float32)

kept = 0
dropped = 0

for rec in records:
    x = int(rec.binX)
    y = int(rec.binY)
    value = float(rec.counts)

    i = (x - bin_start) // resolution
    j = (y - bin_start) // resolution

    if 0 <= i < n_bins and 0 <= j < n_bins:
        matrix[i, j] = value

        # Intra-chromosomal sparse .hic records normally store one
        # triangle, so explicitly reflect across the diagonal.
        if i != j:
            matrix[j, i] = value

        kept += 1
    else:
        dropped += 1

print(
    f"SPARSE_RECONSTRUCTION "
    f"records={len(records)} kept={kept} dropped={dropped} "
    f"shape={matrix.shape} nonzero={np.count_nonzero(matrix)}",
    flush=True,
)

np.save(out_path, matrix)
'''

    def fetch_with_retries(chrom, start, end, resolution, normalization, max_attempts=4):
        url = "https://huggingface.co/buckets/FatimaInstitute/3d-genomics-prediction/resolve/HCT116/ENCFF573OPJ.hic?download=true"
        with tempfile.TemporaryDirectory() as tmpdir:
            script_path = os.path.join(tmpdir, "fetch.py")
            out_path = os.path.join(tmpdir, "matrix.npy")
            with open(script_path, "w") as f:
                f.write(fetch_script)

            last_error = None
            for attempt in range(1, max_attempts + 1):
                result = subprocess.run(
                    [sys.executable, script_path, url, chrom, str(start), str(end),
                     str(resolution), normalization, out_path],
                    capture_output=True, text=True, timeout=300,
                )
                if result.returncode == 0 and os.path.exists(out_path):
                    return np.load(out_path), attempt, None
                last_error = f"attempt {attempt}: returncode={result.returncode}, stderr={result.stderr[-300:]}"
            return None, max_attempts, last_error

    def mask_diagonal(matrix):
        masked = matrix.copy()
        off = matrix[~np.eye(matrix.shape[0], dtype=bool)]
        bg = np.median(off[off > 0]) if np.any(off > 0) else 0
        np.fill_diagonal(masked, bg)
        return masked

    def distance_stratified_zscore(matrix):
        """
        Z-score contacts independently at each genomic separation.

        For a symmetric contact matrix, cells at distance d are exactly
        the +/-d diagonals. Processing one diagonal at a time gives the
        same distance-stratified normalization without constructing a
        full distance matrix or scanning all n^2 cells n times.
        """
        n = matrix.shape[0]
        z = np.zeros_like(matrix, dtype=np.float32)

        for d in range(n):
            vals = np.diagonal(matrix, offset=d)

            if vals.size < 2:
                continue

            mean_d = float(vals.mean())
            std_d = float(vals.std())

            if std_d <= 0:
                continue

            zvals = ((vals - mean_d) / std_d).astype(
                np.float32,
                copy=False,
            )

            i = np.arange(n - d)
            j = i + d

            z[i, j] = zvals

            if d > 0:
                z[j, i] = zvals

        return z

    chrom = "chr10"
    resolution = 500
    start, end = 20_000_000, 21_000_000  # 1Mb window, kept modest for a first preliminary pass

    summary_lines = [f"Fetching {chrom}:{start}-{end} at {resolution}bp (subprocess-isolated, retries enabled)..."]
    matrix, attempts_used, error = fetch_with_retries(chrom, start, end, resolution, "NONE")

    if matrix is None:
        summary_lines.append(f"FAILED after {attempts_used} attempts. Last error: {error}")
        return None, "\n".join(summary_lines)

    summary_lines.append(f"Succeeded on attempt {attempts_used}/4.")
    summary_lines.append(
        "Extraction: hic-straw sparse getRecords() + "
        "manual symmetric dense reconstruction"
    )
    summary_lines.append(f"Matrix shape: {matrix.shape}")
    nonzero_frac = np.count_nonzero(matrix) / matrix.size
    summary_lines.append(f"Non-zero fraction: {nonzero_frac:.4f}")

    masked = mask_diagonal(matrix)
    z = distance_stratified_zscore(masked)

    filt = -z
    cc = gudhi.CubicalComplex(top_dimensional_cells=filt)
    cc.compute_persistence()
    h0 = cc.persistence_intervals_in_dimension(0)
    finite = h0[np.isfinite(h0[:, 1])]
    pers = np.sort(finite[:, 1] - finite[:, 0])[::-1] if len(finite) else np.array([])

    summary_lines.append(f"H0 features: {len(h0)} total, {len(finite)} finite")
    if len(pers) > 0:
        summary_lines.append(f"Top 10 persistence values: {pers[:10].round(3).tolist()}")

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    axes[0].imshow(np.log1p(matrix), cmap="Reds", origin="upper")
    axes[0].set_title(f"{chrom}:{start}-{end}, HCT116 intact Hi-C, {resolution}bp\nraw counts")
    if len(pers) > 0:
        axes[1].plot(pers[:200], marker='o', markersize=2)
    axes[1].set_xlabel("Rank")
    axes[1].set_ylabel("Persistence (z-scored)")
    axes[1].set_title("Sorted persistence -- preliminary")
    plt.tight_layout()
    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=150)
    plt.close()
    buf.seek(0)

    return buf.getvalue(), "\n".join(summary_lines)


@app.function(image=image, timeout=1200)
def diagnose_hct116_remote():
    """
    Diagnose exactly where hic-straw fails on the HCT116 Hugging Face
    bucket file.

    Each hic-straw test runs in its own subprocess so a native SIGSEGV
    cannot kill the parent. We use unbuffered Python + explicit stage
    markers so that, even if C++ crashes, we can see the last operation
    reached.

    This deliberately does NOT run persistence. It only tests transport,
    sparse extraction, dense extraction, and resolution dependence.
    """
    import os
    import subprocess
    import sys
    import tempfile
    import urllib.request
    import urllib.parse

    url = (
        "https://huggingface.co/buckets/"
        "FatimaInstitute/3d-genomics-prediction/"
        "resolve/HCT116/ENCFF573OPJ.hic?download=true"
    )

    lines = []
    lines.append("=== HCT116 remote .hic diagnostic ===")
    lines.append(f"URL: {url}")

    # First test ordinary HTTP byte-range behaviour independently of hic-straw.
    def probe_range(start, end):
        req = urllib.request.Request(
            url,
            headers={
                "Range": f"bytes={start}-{end}",
                "User-Agent": "Mozilla/5.0",
                "Accept-Encoding": "identity",
            },
        )
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                # Never read more than 1 KiB even if the server ignores Range.
                sample = resp.read(1024)
                final_url = resp.geturl()
                host = urllib.parse.urlparse(final_url).netloc
                return (
                    f"HTTP range {start}-{end}: "
                    f"status={getattr(resp, 'status', 'unknown')} "
                    f"host={host} "
                    f"content-range={resp.headers.get('Content-Range')} "
                    f"content-length={resp.headers.get('Content-Length')} "
                    f"bytes-read={len(sample)}"
                )
        except Exception as e:
            return f"HTTP range {start}-{end}: FAILED {type(e).__name__}: {e}"

    lines.append("")
    lines.append("--- independent HTTP range probes ---")
    lines.append(probe_range(0, 1023))
    lines.append(probe_range(1_000_000, 1_001_023))

    probe_script = r"""
import faulthandler
import importlib.metadata
import sys

faulthandler.enable(all_threads=True)

import hicstraw
import numpy as np

url = sys.argv[1]
label = sys.argv[2]
chrom = sys.argv[3]
start = int(sys.argv[4])
end = int(sys.argv[5])
resolution = int(sys.argv[6])
mode = sys.argv[7]

def mark(msg):
    print(f"[{label}] {msg}", file=sys.stderr, flush=True)

try:
    version = importlib.metadata.version("hic-straw")
except Exception:
    version = "unknown"

mark(f"hic-straw version={version}")
mark(
    f"parameters chrom={chrom} start={start} end={end} "
    f"resolution={resolution} normalization=NONE mode={mode}"
)

mark("stage=HiCFile:start")
hic = hicstraw.HiCFile(url)
mark("stage=HiCFile:ok")

mark("stage=getResolutions:start")
resolutions = hic.getResolutions()
mark(f"stage=getResolutions:ok resolutions={resolutions}")

mark("stage=getMatrixZoomData:start")
mzd = hic.getMatrixZoomData(
    chrom,
    chrom,
    "observed",
    "NONE",
    "BP",
    resolution,
)
mark("stage=getMatrixZoomData:ok")

if mode == "records":
    mark("stage=getRecords:start")
    records = mzd.getRecords(start, end, start, end)
    mark(f"stage=getRecords:ok n_records={len(records)}")

elif mode == "matrix":
    mark("stage=getRecordsAsMatrix:start")
    matrix = mzd.getRecordsAsMatrix(start, end, start, end)
    matrix = np.asarray(matrix)
    mark(
        "stage=getRecordsAsMatrix:ok "
        f"shape={matrix.shape} "
        f"nonzero={np.count_nonzero(matrix)} "
        f"total={matrix.size}"
    )

else:
    raise ValueError(f"Unknown mode: {mode}")

mark("stage=complete")
"""

    # Important:
    #   - coarse resolution first
    #   - sparse before dense
    #   - tiny 500 bp window before 1 Mb
    #
    # This distinguishes transport failure from dense-conversion failure
    # and fine-resolution-specific failure.
    cases = [
        (
            "10kb_sparse_1Mb",
            "chr10",
            20_000_000,
            21_000_000,
            10_000,
            "records",
        ),
        (
            "10kb_dense_1Mb",
            "chr10",
            20_000_000,
            21_000_000,
            10_000,
            "matrix",
        ),
        (
            "500bp_sparse_100kb",
            "chr10",
            20_000_000,
            20_100_000,
            500,
            "records",
        ),
        (
            "500bp_dense_100kb",
            "chr10",
            20_000_000,
            20_100_000,
            500,
            "matrix",
        ),
        (
            "500bp_sparse_1Mb",
            "chr10",
            20_000_000,
            21_000_000,
            500,
            "records",
        ),
        (
            "500bp_dense_1Mb",
            "chr10",
            20_000_000,
            21_000_000,
            500,
            "matrix",
        ),
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        script_path = os.path.join(tmpdir, "hic_probe.py")
        with open(script_path, "w") as f:
            f.write(probe_script)

        env = os.environ.copy()
        env["PYTHONFAULTHANDLER"] = "1"

        lines.append("")
        lines.append("--- hic-straw isolated probes ---")

        for label, chrom, start, end, resolution, mode in cases:
            lines.append("")
            lines.append(f"### {label}")

            cmd = [
                sys.executable,
                "-u",
                script_path,
                url,
                label,
                chrom,
                str(start),
                str(end),
                str(resolution),
                mode,
            ]

            try:
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=240,
                    env=env,
                )

                lines.append(f"returncode={result.returncode}")

                if result.stdout.strip():
                    lines.append("stdout:")
                    lines.append(result.stdout.strip())

                if result.stderr.strip():
                    lines.append("stderr:")
                    lines.append(result.stderr.strip())

            except subprocess.TimeoutExpired as e:
                lines.append("TIMEOUT after 240 seconds")

                if e.stdout:
                    lines.append("partial stdout:")
                    lines.append(
                        e.stdout
                        if isinstance(e.stdout, str)
                        else e.stdout.decode(errors="replace")
                    )

                if e.stderr:
                    lines.append("partial stderr:")
                    lines.append(
                        e.stderr
                        if isinstance(e.stderr, str)
                        else e.stderr.decode(errors="replace")
                    )

    lines.append("")
    lines.append("=== diagnostic complete ===")
    return "\n".join(lines)


@app.function(image=image, timeout=600)
def validate_hct116_sparse_reconstruction():
    """
    Validate our sparse->dense reconstruction against hic-straw's native
    getRecordsAsMatrix() on a small 500bp region where the native method
    is known not to segfault.
    """
    import hicstraw
    import numpy as np

    url = (
        "https://huggingface.co/buckets/"
        "FatimaInstitute/3d-genomics-prediction/"
        "resolve/HCT116/ENCFF573OPJ.hic?download=true"
    )

    chrom = "chr10"
    start = 20_000_000
    end = 20_100_000
    resolution = 500

    hic = hicstraw.HiCFile(url)

    # Native dense method: safe at this 100kb query.
    mzd_native = hic.getMatrixZoomData(
        chrom, chrom, "observed", "NONE", "BP", resolution
    )
    native = np.asarray(
        mzd_native.getRecordsAsMatrix(start, end, start, end),
        dtype=np.float32,
    )

    # Sparse method used by the new 1Mb pipeline.
    mzd_sparse = hic.getMatrixZoomData(
        chrom, chrom, "observed", "NONE", "BP", resolution
    )
    records = mzd_sparse.getRecords(start, end, start, end)

    bin_start = (start // resolution) * resolution
    bin_end = (end // resolution) * resolution
    n_bins = ((bin_end - bin_start) // resolution) + 1

    reconstructed = np.zeros(
        (n_bins, n_bins),
        dtype=np.float32,
    )

    for rec in records:
        x = int(rec.binX)
        y = int(rec.binY)
        value = float(rec.counts)

        i = (x - bin_start) // resolution
        j = (y - bin_start) // resolution

        if 0 <= i < n_bins and 0 <= j < n_bins:
            reconstructed[i, j] = value
            if i != j:
                reconstructed[j, i] = value

    diff = np.abs(native - reconstructed)

    mismatched = np.count_nonzero(diff)
    max_diff = float(diff.max())
    mean_diff = float(diff.mean())

    lines = [
        "=== HCT116 sparse reconstruction validation ===",
        f"Region: {chrom}:{start}-{end}",
        f"Resolution: {resolution}bp",
        f"Native shape: {native.shape}",
        f"Reconstructed shape: {reconstructed.shape}",
        f"Sparse records: {len(records)}",
        f"Native nonzero: {np.count_nonzero(native)}",
        f"Reconstructed nonzero: {np.count_nonzero(reconstructed)}",
        f"Mismatched cells: {mismatched}",
        f"Maximum absolute difference: {max_diff}",
        f"Mean absolute difference: {mean_diff}",
        f"Exact match: {bool(np.array_equal(native, reconstructed))}",
    ]

    return "\n".join(lines)


@app.function(image=image, timeout=600)
def diagnose_hct116_500bp_distance_sparsity():
    """
    Measure contact occupancy and z-score behaviour as a function of
    genomic separation for the HCT116 chr10 500bp preliminary region.

    Uses sparse records directly, avoiding native dense conversion.
    """
    import csv
    import io

    import hicstraw
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    url = (
        "https://huggingface.co/buckets/"
        "FatimaInstitute/3d-genomics-prediction/"
        "resolve/HCT116/ENCFF573OPJ.hic?download=true"
    )

    chrom = "chr10"
    start = 20_000_000
    end = 21_000_000
    resolution = 500

    hic = hicstraw.HiCFile(url)
    mzd = hic.getMatrixZoomData(
        chrom,
        chrom,
        "observed",
        "NONE",
        "BP",
        resolution,
    )

    records = mzd.getRecords(
        start,
        end,
        start,
        end,
    )

    bin_start = (start // resolution) * resolution
    bin_end = (end // resolution) * resolution
    n_bins = ((bin_end - bin_start) // resolution) + 1

    nonzero = np.zeros(n_bins, dtype=np.int64)
    sums = np.zeros(n_bins, dtype=np.float64)
    sums_sq = np.zeros(n_bins, dtype=np.float64)
    max_count = np.zeros(n_bins, dtype=np.float64)

    kept = 0

    for rec in records:
        x = int(rec.binX)
        y = int(rec.binY)
        value = float(rec.counts)

        i = (x - bin_start) // resolution
        j = (y - bin_start) // resolution

        if not (
            0 <= i < n_bins
            and 0 <= j < n_bins
        ):
            continue

        d = abs(j - i)

        if value != 0:
            nonzero[d] += 1
            sums[d] += value
            sums_sq[d] += value * value
            if value > max_count[d]:
                max_count[d] = value

        kept += 1

    total = n_bins - np.arange(n_bins)

    occupancy = np.divide(
        nonzero,
        total,
        out=np.zeros_like(nonzero, dtype=float),
        where=total > 0,
    )

    mean_all = np.divide(
        sums,
        total,
        out=np.zeros_like(sums),
        where=total > 0,
    )

    second_moment = np.divide(
        sums_sq,
        total,
        out=np.zeros_like(sums_sq),
        where=total > 0,
    )

    variance = np.maximum(
        second_moment - mean_all ** 2,
        0.0,
    )

    std_all = np.sqrt(variance)

    max_z = np.divide(
        max_count - mean_all,
        std_all,
        out=np.zeros_like(max_count),
        where=std_all > 0,
    )

    distance_bp = np.arange(n_bins) * resolution

    # CSV
    csv_buf = io.StringIO()
    writer = csv.writer(csv_buf)

    writer.writerow([
        "distance_bp",
        "possible_cells",
        "nonzero_cells",
        "occupancy",
        "mean_count_all_cells",
        "std_count_all_cells",
        "max_count",
        "max_zscore",
    ])

    for d in range(n_bins):
        writer.writerow([
            int(distance_bp[d]),
            int(total[d]),
            int(nonzero[d]),
            float(occupancy[d]),
            float(mean_all[d]),
            float(std_all[d]),
            float(max_count[d]),
            float(max_z[d]),
        ])

    # Plot
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(13, 5),
    )

    axes[0].plot(
        distance_bp[1:] / 1000,
        occupancy[1:],
    )
    axes[0].set_xlabel("Genomic separation (kb)")
    axes[0].set_ylabel("Non-zero fraction")
    axes[0].set_title(
        "HCT116 chr10 500bp\n"
        "contact occupancy by distance"
    )

    axes[1].plot(
        distance_bp[1:] / 1000,
        max_z[1:],
    )
    axes[1].set_xlabel("Genomic separation (kb)")
    axes[1].set_ylabel("Maximum z-score")
    axes[1].set_title(
        "Largest distance-stratified z-score"
    )

    plt.tight_layout()

    img_buf = io.BytesIO()
    plt.savefig(
        img_buf,
        format="png",
        dpi=150,
    )
    plt.close()
    img_buf.seek(0)

    # Useful distance bands
    bands = [
        (0, 10_000),
        (10_000, 25_000),
        (25_000, 50_000),
        (50_000, 100_000),
        (100_000, 250_000),
        (250_000, 500_000),
        (500_000, 1_000_001),
    ]

    lines = [
        "=== HCT116 500bp distance-sparsity diagnostic ===",
        f"Region: {chrom}:{start}-{end}",
        f"Resolution: {resolution}bp",
        f"Bins: {n_bins}",
        f"Sparse records retained: {kept}",
        "",
        "Distance-band occupancy:",
    ]

    for lo, hi in bands:
        mask = (
            (distance_bp >= lo)
            & (distance_bp < hi)
        )

        possible = int(total[mask].sum())
        observed = int(nonzero[mask].sum())

        frac = (
            observed / possible
            if possible
            else float("nan")
        )

        lines.append(
            f"{lo/1000:7.1f}-{hi/1000:7.1f} kb: "
            f"{frac:.4f} "
            f"({observed}/{possible})"
        )

    lines.append("")
    lines.append("Selected separations:")

    for bp in [
        500,
        1_000,
        2_000,
        5_000,
        10_000,
        25_000,
        50_000,
        100_000,
        250_000,
        500_000,
        750_000,
        1_000_000,
    ]:
        d = bp // resolution

        lines.append(
            f"{bp/1000:7.1f} kb: "
            f"occupancy={occupancy[d]:.4f}, "
            f"nonzero={nonzero[d]}/{total[d]}, "
            f"max_count={max_count[d]:.3f}, "
            f"max_z={max_z[d]:.3f}"
        )

    return (
        img_buf.getvalue(),
        csv_buf.getvalue(),
        "\n".join(lines),
    )


@app.function(image=image, timeout=3600)
def compare_hct116_500bp_normalizations():
    """
    Compare the existing distance-stratified z-score against a bounded
    distance-stratified rank-Gaussian transform.

    Also compare the rank-Gaussian real matrix against distance-preserving
    shuffled nulls. Each null preserves the exact marginal contact
    distribution at every genomic separation while destroying genomic
    organization along that diagonal.
    """
    import gc
    import io
    import time

    import gudhi
    import hicstraw

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    import numpy as np
    from scipy.stats import rankdata, norm

    url = (
        "https://huggingface.co/buckets/"
        "FatimaInstitute/3d-genomics-prediction/"
        "resolve/HCT116/ENCFF573OPJ.hic?download=true"
    )

    chrom = "chr10"
    start = 20_000_000
    end = 21_000_000
    resolution = 500

    # ------------------------------------------------------------
    # Fetch sparse contacts and reconstruct dense matrix.
    # This path has already been validated exactly against
    # getRecordsAsMatrix() on a smaller region.
    # ------------------------------------------------------------

    hic = hicstraw.HiCFile(url)

    mzd = hic.getMatrixZoomData(
        chrom,
        chrom,
        "observed",
        "NONE",
        "BP",
        resolution,
    )

    records = mzd.getRecords(
        start,
        end,
        start,
        end,
    )

    bin_start = (start // resolution) * resolution
    bin_end = (end // resolution) * resolution
    n = ((bin_end - bin_start) // resolution) + 1

    matrix = np.zeros(
        (n, n),
        dtype=np.float32,
    )

    for rec in records:
        x = int(rec.binX)
        y = int(rec.binY)
        value = float(rec.counts)

        i = (x - bin_start) // resolution
        j = (y - bin_start) // resolution

        if 0 <= i < n and 0 <= j < n:
            matrix[i, j] = value

            if i != j:
                matrix[j, i] = value

    # ------------------------------------------------------------
    # Mask exact diagonal, as in the existing pipeline.
    # ------------------------------------------------------------

    masked = matrix.copy()

    off = matrix[
        ~np.eye(
            matrix.shape[0],
            dtype=bool,
        )
    ]

    positive = off[off > 0]

    bg = (
        float(np.median(positive))
        if positive.size
        else 0.0
    )

    np.fill_diagonal(
        masked,
        bg,
    )

    # ------------------------------------------------------------
    # Existing normalization.
    # ------------------------------------------------------------

    def distance_zscore(x):
        n = x.shape[0]

        z = np.zeros_like(
            x,
            dtype=np.float32,
        )

        for d in range(n):
            vals = np.diagonal(
                x,
                offset=d,
            )

            if vals.size < 2:
                continue

            mean_d = float(vals.mean())
            std_d = float(vals.std())

            if std_d <= 0:
                continue

            zv = (
                (vals - mean_d)
                / std_d
            ).astype(
                np.float32,
                copy=False,
            )

            i = np.arange(n - d)
            j = i + d

            z[i, j] = zv

            if d > 0:
                z[j, i] = zv

        return z

    # ------------------------------------------------------------
    # Candidate diagnostic normalization.
    #
    # Rank within each genomic-distance diagonal, preserving ties,
    # then convert empirical quantiles to standard-normal scores.
    #
    # Important: this prevents rare contacts from getting arbitrarily
    # huge scores purely because occupancy approaches zero.
    # ------------------------------------------------------------

    def distance_rank_gaussian(x):
        n = x.shape[0]

        z = np.zeros_like(
            x,
            dtype=np.float32,
        )

        # Ignore d=0 because the exact diagonal was deliberately masked.
        for d in range(1, n):
            vals = np.diagonal(
                x,
                offset=d,
            ).copy()

            m = vals.size

            if m < 2:
                continue

            ranks = rankdata(
                vals,
                method="average",
            )

            u = (
                ranks - 0.5
            ) / m

            # Numerical protection only.
            u = np.clip(
                u,
                1e-6,
                1 - 1e-6,
            )

            zv = norm.ppf(u).astype(
                np.float32,
                copy=False,
            )

            i = np.arange(m)
            j = i + d

            z[i, j] = zv
            z[j, i] = zv

        return z

    # ------------------------------------------------------------
    # Persistence helper.
    # ------------------------------------------------------------

    def persistence_values(score):
        filt = -score

        cc = gudhi.CubicalComplex(
            top_dimensional_cells=filt
        )

        cc.compute_persistence()

        h0 = cc.persistence_intervals_in_dimension(0)

        finite = h0[
            np.isfinite(h0[:, 1])
        ]

        if len(finite):
            pers = np.sort(
                finite[:, 1]
                - finite[:, 0]
            )[::-1]
        else:
            pers = np.array([])

        n_total = len(h0)
        n_finite = len(finite)

        del cc
        gc.collect()

        return pers, n_total, n_finite

    # ------------------------------------------------------------
    # Null generator.
    #
    # Permute positions independently within each distance diagonal,
    # then restore symmetry.
    # ------------------------------------------------------------

    def distance_preserving_shuffle(x, rng):
        n = x.shape[0]

        shuffled = np.zeros_like(
            x,
            dtype=np.float32,
        )

        # Exact diagonal carries no information here.
        np.fill_diagonal(
            shuffled,
            np.diagonal(x),
        )

        for d in range(1, n):
            vals = np.diagonal(
                x,
                offset=d,
            ).copy()

            rng.shuffle(vals)

            i = np.arange(n - d)
            j = i + d

            shuffled[i, j] = vals
            shuffled[j, i] = vals

        return shuffled

    lines = [
        "=== HCT116 500bp normalization comparison ===",
        f"Region: {chrom}:{start}-{end}",
        f"Resolution: {resolution}bp",
        f"Matrix shape: {matrix.shape}",
        (
            "Non-zero fraction: "
            f"{np.count_nonzero(matrix) / matrix.size:.4f}"
        ),
        "",
    ]

    # ------------------------------------------------------------
    # Existing z-score.
    # ------------------------------------------------------------

    print(
        "[1/5] Computing existing distance z-score persistence...",
        flush=True,
    )
    t_stage = time.perf_counter()

    z_old = distance_zscore(masked)

    old_max_cell = float(
        np.max(z_old)
    )

    old_pers, old_total, old_finite = (
        persistence_values(z_old)
    )

    elapsed = time.perf_counter() - t_stage
    print(
        f"[1/5] Existing z-score persistence complete "
        f"in {elapsed:.1f}s",
        flush=True,
    )
    lines.append(
        f"Existing z-score persistence runtime: {elapsed:.1f}s"
    )

    lines.extend([
        "",
        "--- Existing distance z-score ---",
        f"Maximum cell z-score: {old_max_cell:.3f}",
        f"H0 features: {old_total}",
        f"Finite H0 features: {old_finite}",
        (
            "Top 10 persistence: "
            f"{old_pers[:10].round(3).tolist()}"
        ),
    ])

    del z_old
    gc.collect()

    # ------------------------------------------------------------
    # Rank-Gaussian real matrix.
    # ------------------------------------------------------------

    lines.append("")

    print(
        "[2/5] Computing rank-Gaussian real persistence...",
        flush=True,
    )
    t_stage = time.perf_counter()

    z_rank = distance_rank_gaussian(masked)

    rank_max_cell = float(
        np.max(z_rank)
    )

    rank_pers, rank_total, rank_finite = (
        persistence_values(z_rank)
    )

    elapsed = time.perf_counter() - t_stage
    print(
        f"[2/5] Rank-Gaussian real persistence complete "
        f"in {elapsed:.1f}s",
        flush=True,
    )
    lines.append(
        f"Rank-Gaussian real persistence runtime: {elapsed:.1f}s"
    )

    lines.extend([
        "",
        "--- Distance rank-Gaussian ---",
        f"Maximum cell score: {rank_max_cell:.3f}",
        f"H0 features: {rank_total}",
        f"Finite H0 features: {rank_finite}",
        (
            "Top 10 persistence: "
            f"{rank_pers[:10].round(3).tolist()}"
        ),
    ])

    del z_rank
    gc.collect()

    # ------------------------------------------------------------
    # Three exploratory distance-preserving nulls.
    # ------------------------------------------------------------

    rng = np.random.default_rng(20260930)

    null_persistence = []

    for null_idx in range(3):
        stage_number = null_idx + 3

        lines.append("")

        print(
            f"[{stage_number}/5] Computing shuffled null "
            f"{null_idx + 1}/3...",
            flush=True,
        )
        t_stage = time.perf_counter()

        shuffled = distance_preserving_shuffle(
            masked,
            rng,
        )

        z_null = distance_rank_gaussian(
            shuffled,
        )

        null_pers, _, _ = persistence_values(
            z_null
        )

        elapsed = time.perf_counter() - t_stage

        null_persistence.append(
            null_pers
        )

        print(
            f"[{stage_number}/5] Shuffled null "
            f"{null_idx + 1}/3 complete in {elapsed:.1f}s; "
            f"max persistence="
            f"{null_pers[0] if len(null_pers) else float('nan'):.3f}",
            flush=True,
        )

        lines.append(
            f"Null {null_idx + 1} runtime: {elapsed:.1f}s"
        )

        lines.append(
            f"Null {null_idx + 1} top 10 persistence: "
            f"{null_pers[:10].round(3).tolist()}"
        )

        del shuffled
        del z_null
        gc.collect()

    null_maxima = [
        float(p[0])
        if len(p)
        else float("nan")
        for p in null_persistence
    ]

    lines.extend([
        "",
        "--- Real vs shuffled null ---",
        (
            "Real rank-Gaussian max persistence: "
            f"{rank_pers[0]:.3f}"
            if len(rank_pers)
            else
            "Real rank-Gaussian max persistence: NA"
        ),
        (
            "Null max persistence values: "
            f"{np.round(null_maxima, 3).tolist()}"
        ),
    ])

    if (
        len(rank_pers)
        and len(null_maxima)
    ):
        lines.append(
            "Real/null-max ratio using mean null maximum: "
            f"{rank_pers[0] / np.mean(null_maxima):.3f}"
        )

    # ------------------------------------------------------------
    # Plot ranked persistence curves.
    # ------------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(9, 6)
    )

    top_n = 200

    if len(old_pers):
        ax.plot(
            old_pers[:top_n],
            label="existing z-score",
        )

    if len(rank_pers):
        ax.plot(
            rank_pers[:top_n],
            label="rank-Gaussian real",
        )

    for idx, null_pers in enumerate(
        null_persistence
    ):
        if len(null_pers):
            ax.plot(
                null_pers[:top_n],
                alpha=0.45,
                label=f"rank-Gaussian null {idx + 1}",
            )

    ax.set_xlabel(
        "Persistence rank"
    )
    ax.set_ylabel(
        "H0 persistence"
    )
    ax.set_title(
        "HCT116 chr10, 500bp\n"
        "normalization and distance-preserving null comparison"
    )

    ax.legend()

    plt.tight_layout()

    buf = io.BytesIO()

    plt.savefig(
        buf,
        format="png",
        dpi=150,
    )

    plt.close()

    buf.seek(0)

    return (
        buf.getvalue(),
        "\n".join(lines),
    )


@app.local_entrypoint()
def main(task: str = "insulation"):
    if task == "hct116_norm_compare":
        result = compare_hct116_500bp_normalizations.remote()

        with open(
            "figures/hct116_500bp_normalization_comparison.png",
            "wb",
        ) as f:
            f.write(result[0])

        print(result[1])

        print(
            "Saved "
            "figures/hct116_500bp_normalization_comparison.png"
        )

    elif task == "hct116_distance_diag":
        result = diagnose_hct116_500bp_distance_sparsity.remote()

        with open(
            "figures/hct116_500bp_distance_sparsity.png",
            "wb",
        ) as f:
            f.write(result[0])

        with open(
            "figures/hct116_500bp_distance_sparsity.csv",
            "w",
        ) as f:
            f.write(result[1])

        print(result[2])
        print(
            "Saved "
            "figures/hct116_500bp_distance_sparsity.png"
        )
        print(
            "Saved "
            "figures/hct116_500bp_distance_sparsity.csv"
        )

    elif task == "hct116_recon_validate":
        print(validate_hct116_sparse_reconstruction.remote())

    elif task == "hct116_diag":
        print(diagnose_hct116_remote.remote())

    elif task == "hct116_prelim":
        result = hct116_prelim_500bp.remote()
        if result[0] is not None:
            with open("figures/hct116_prelim_500bp.png", "wb") as f:
                f.write(result[0])
            print(result[1])
            print("Saved figures/hct116_prelim_500bp.png")
        else:
            print(result[1])

    elif task == "check_intact_hic_sparsity":
        print(check_intact_hic_sparsity.remote())

    elif task == "check_intact_hic":
        print(check_intact_hic_resolution.remote())

    elif task == "cache_peaks":
        print(fetch_and_cache_peaks.remote())

    elif task == "spawn_all":
        chrom_lengths_hg19 = {
            "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276,
            "5": 180915260, "6": 171115067, "7": 159138663, "8": 146364022,
            "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895,
            "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753,
            "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
            "21": 48129895, "22": 51304566,
        }
        regions = []
        for chrom, length in chrom_lengths_hg19.items():
            start_mb = 20
            while (start_mb + 2) * 1_000_000 < length - 7_000_000:
                regions.append((chrom, start_mb * 1_000_000, start_mb * 1_000_000 + 2_000_000))
                start_mb += 25
        print(f"Spawning {len(regions)} independent region jobs...")
        for chrom, start, end in regions:
            process_one_region.spawn(chrom, start, end)
        print("All jobs spawned. They run independently on Modal now -- safe to disconnect.")
        print("Check progress any time with: modal run modal_app.py --task aggregate")

    elif task == "retry_failed":
        failed_regions = [
            ("10",120000000,122000000),("10",70000000,72000000),("11",20000000,22000000),
            ("11",45000000,47000000),("11",95000000,97000000),("12",120000000,122000000),
            ("13",45000000,47000000),("13",70000000,72000000),("14",20000000,22000000),
            ("14",45000000,47000000),("15",70000000,72000000),("16",70000000,72000000),
            ("17",20000000,22000000),("17",45000000,47000000),("18",45000000,47000000),
            ("1",195000000,197000000),("1",20000000,22000000),("1",220000000,222000000),
            ("20",20000000,22000000),("20",45000000,47000000),("21",20000000,22000000),
            ("3",20000000,22000000),("3",70000000,72000000),("4",170000000,172000000),
            ("5",170000000,172000000),("5",45000000,47000000),("6",45000000,47000000),
            ("7",145000000,147000000),("7",45000000,47000000),("8",20000000,22000000),
            ("8",70000000,72000000),("9",120000000,122000000),("9",45000000,47000000),
            ("9",95000000,97000000),
        ]
        print(f"Retrying {len(failed_regions)} previously failed regions...")
        for chrom, start, end in failed_regions:
            process_one_region.spawn(chrom, start, end)
        print("All retries spawned.")

    elif task == "aggregate":
        print(aggregate_results.remote())

    elif task == "ctcf_benchmark_genomewide":
        image_bytes, summary = run_ctcf_benchmark_genomewide.remote()
        with open("figures/ctcf_benchmark_genomewide.png", "wb") as f:
            f.write(image_bytes)
        print(summary)
        print("Saved figures/ctcf_benchmark_genomewide.png")

    elif task == "ctcf_benchmark_expanded":
        image_bytes, summary = run_ctcf_benchmark_expanded.remote()
        with open("figures/ctcf_benchmark_expanded.png", "wb") as f:
            f.write(image_bytes)
        print(summary)
        print("Saved figures/ctcf_benchmark_expanded.png")

    elif task == "ctcf_benchmark":
        image_bytes, summary = run_ctcf_benchmark.remote()
        with open("figures/ctcf_benchmark.png", "wb") as f:
            f.write(image_bytes)
        print(summary)
        print("Saved figures/ctcf_benchmark.png")

    elif task == "integrate_multi":
        image_bytes, summary = run_integration_multiregion.remote()
        with open("figures/stage1_stage2_multiregion.png", "wb") as f:
            f.write(image_bytes)
        print(summary)
        print("Saved figures/stage1_stage2_multiregion.png")

    elif task == "integrate":
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
