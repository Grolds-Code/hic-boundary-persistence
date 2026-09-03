# Project Journal

A chronological record of what was tried, what broke, and why -- kept
separately from the paper and README so the reasoning trail stays
readable even as the code itself keeps changing. Nothing here gets
deleted or rewritten as the project progresses; new entries are added
to the bottom.

## 1. Synthetic validation -- a real bug, and the fix

Before trusting the Stage 2 persistent homology method on real data,
it was tested against a synthetic contact matrix with known, planted
TAD (Topologically Associating Domain) boundaries.

The first version had a real bug: a naive cubical filtration on a
full block-diagonal matrix trivially merges every planted block into
one connected component from birth. This happens because every
diagonal pixel is, by construction, part of its own block, and
consecutive diagonal pixels are always corner-adjacent in a cubical
complex -- which is enough to connect two cells regardless of block
structure. The fix, masking the exact diagonal before filtration,
turned out not to be a special-purpose patch: real Hi-C analyses
already exclude the exact diagonal as standard practice, since
self-contacts are a known assay artifact, not real biological signal.

After the fix, the method recovered all 4 planted boundaries in
synthetic tests, separated from noise by close to an order of
magnitude in persistence. See src/synthetic_validation.py.

## 2. Real data, Finding 1 -- raw counts are dominated by distance decay

Running the corrected method on real GM12878 data (chr21) surfaced a
new problem synthetic data could not reveal, since the synthetic
model had no analogue of genomic distance decay. Every top-ranked
"boundary" sat directly on the diagonal, tracking pure distance decay
(nearby DNA always touches more) rather than real domain structure.

Fix: switch from raw observed contact counts to observed/expected
(O/E) values, dividing out the expected distance-dependent background
at each genomic distance.

## 3. Real data, Finding 2 -- O/E at long range is dominated by noise

Switching to O/E moved the problem rather than solving it: top-ranked
features jumped to the far corners of the analysis window -- the
longest-range pairs within it, regardless of window size. At long
genomic distance, expected counts shrink toward zero, so a single
stray read produces a large, spurious O/E ratio.

Fix (at the time): restrict the persistence computation to a bounded
genomic distance band (~300kb) around the diagonal, independent of
the overall window size. See src/real_data_persistence.py.

## 4. Distance-band sensitivity sweep -- the 300kb cutoff does not hold up

Limitation 6 in the paper flagged the 300kb cutoff as chosen by direct
observation on one region, not validated. Testing it properly (first
locally, then extended via Modal to include lower cutoffs) swept
cutoffs from 20kb to 1000kb on the same region.

Result: top persistence rises monotonically at every single cutoff
tested -- 0.31 (20kb) / 0.38 (50kb) / 0.77 (100kb) / 1.09 (200kb) /
1.42 (300kb) / 2.78 (500kb) / 3.82 (1000kb) -- with no plateau
anywhere, not even near 20kb, close to the tightest possible band.
This is a stronger, more honest finding than "not yet validated": a
hard distance-band cutoff does not have a correct value to find. It
postpones the long-range noise problem rather than removing it. See
src/distance_band_sweep.py and modal_app.py (function run_sweep).

## 5. Distance-stratified z-score -- testing a more principled fix

Since a hard cutoff does not resolve the underlying issue (variance in
O/E ratios growing continuously with genomic distance), the next test
was a distance-stratified z-score: normalize each bin-pair by both the
mean AND standard deviation of values at that specific genomic
distance, rather than discarding data past a threshold. If this is
the right fix, top persistence should stay roughly flat as the
analysis window grows, unlike the cutoff approach, which kept
climbing.

Result: top persistence across window sizes 1Mb / 2Mb / 3Mb / 5Mb was
5.64 / 5.56 / 5.73 / 8.17. This is a real improvement over the cutoff
approach -- roughly flat from 1Mb to 3Mb, versus the cutoff approach's
steady climb (about 12x from 20kb to 1Mb) over a comparable range. It
is not perfectly scale-invariant, though: the jump at 5Mb is real and
not yet explained -- possibly genuine larger-scale structure, an edge
effect, or still-thin sample counts at the most extreme distances even
in a larger window. Worth investigating further before treating the
z-score approach as fully validated, but it is a clear improvement
over a hard cutoff and the more principled fix going forward. See
modal_app.py (function run_zscore_comparison).

## Infrastructure notes, kept here so they are not rediscovered the hard way

- Windows cannot compile cooltools or hic-straw from source (MSVC
  lacks a random() symbol these C extensions expect) -- this is why
  the project runs in WSL2, not native Windows.
- hic-straw needs libcurl4-openssl-dev (local builds) and
  additionally zlib1g-dev (Modal's minimal base image) as system
  prerequisites before it will compile.
- Modal compute is reserved for anything that would otherwise use
  significant local bandwidth or compute time; single small-region
  pulls are safe to run locally.

## 6. Stage 1 -- insulation score, and a fixed-threshold trap

src/insulation.py had been an empty stub since the start of the
project. Built the standard method (Crane et al. 2015): a sliding
window along the diagonal, boundaries called at local minima of the
log2 insulation score.

Validated on synthetic data first (same discipline as Stage 2): all 4
planted boundaries recovered correctly, ranked correctly by
prominence, no false positives ahead of them.

On real data, the first run found 0 boundaries. Diagnosis (via a
dedicated diagnostic Modal function, not guessing) showed why: the
prominence threshold (0.15) was tuned against the synthetic test's
deliberately dramatic 10:1 domain-to-background contrast, which
produced large prominence values. Real biological insulation dips are
far subtler -- the maximum possible prominence in the real data tested
was 0.0707, less than half the threshold, so nothing could ever pass.

Fix: switched from a fixed absolute prominence threshold to a
percentile-based one (keep the top half of detected dips, regardless
of their absolute scale). This adapts to whatever data it is run on,
rather than requiring hand-tuning per dataset. After the fix: 27
candidate boundaries found on the same chr21 2Mb region. See
src/insulation.py and modal_app.py (function run_insulation).

General lesson: any threshold validated only on synthetic data should
be treated as unvalidated for real data until checked directly --
this is the second time in this project a synthetic-only parameter
choice failed silently on real data (the first being the 300kb
distance-band cutoff).

## 7. First real attempt at Claim 1 -- inconclusive, and why that's honest

Built the actual integration this project has been building toward:
for each pair of consecutive Stage 1 boundaries, treat the region
between them as a candidate domain, run Stage 2 persistence within it,
and compare the resulting ranking against Stage 1's own insulation-
score prominence ranking for the same domains. This is Claim 1's
actual question, tested for the first time.

First pass (28 domains, one 2Mb region): Spearman rho=-0.205, p=0.296.
After dropping domains under 5 bins (too small for persistence to mean
much, likely just noise): 15 domains, rho=0.064, p=0.819.

Neither result is statistically distinguishable from zero correlation.
This is not evidence that persistence and insulation score agree or
disagree -- it means one 2Mb window is not enough data to say anything
yet. 15-28 domains is a pilot, not a sample. The honest conclusion at
this point is that Claim 1 needs a real sample: multiple regions,
likely multiple chromosomes, before the correlation (or lack of one)
means anything. See modal_app.py (function run_integration).

## 8. Multi-region Claim 1 test -- first statistically significant result

Pooled the Stage1-vs-Stage2 domain comparison across 8 regions (5
windows on chr21, 3 on chr20) instead of one, to get real statistical
power. 107 domains total.

Result: rho=0.234, p=0.0155 -- statistically significant at
alpha=0.05, unlike either single-region attempt. The correlation is
real but weak: rho^2 suggests only about 5% of the variance in one
method is explained by the other. Read together, this means the two
methods are not independent (makes sense, both track real boundary
strength to some extent) but are also far from redundant -- a strong
signal that persistence is capturing something insulation-score
prominence does not, which is the core premise Claim 1 is built on.

Important caveat: this correlation tells us how similar the two
methods are, not which one is more biologically accurate. That
question needs the actual CTCF/cohesin enrichment benchmark (paper
Section 4.5), not yet built. See modal_app.py (function
run_integration_multiregion).

## 9. First real test against biological ground truth -- promising, not yet proven

Built the actual CTCF/cohesin benchmark from the paper's Section 4.5:
fetched real ChIP-seq peaks for GM12878 (CTCF: ENCFF833FTF, RAD21:
ENCFF753RGL, SMC3: ENCFF572RPI, all confirmed hg19 to match the Hi-C
data -- one candidate CTCF file, ENCFF797SDL, was correctly rejected
after checking, since it turned out to be hg38). For each of the 107
domains from the multi-region test, checked whether a real peak sits
near the boundary, then tested whether Stage 1 (insulation prominence)
or Stage 2 (persistence) ranking better predicts that.

Point estimates favor persistence on both markers: CTCF AUC 0.646
(Stage 2) vs 0.561 (Stage 1); cohesin AUC 0.572 (Stage 2) vs 0.473
(Stage 1) -- notably, Stage 1 performed at or below random chance
(AUC 0.5) for predicting cohesin binding.

However, a paired bootstrap 95% confidence interval on the AUC
difference includes zero for both markers (CTCF: [-0.063, 0.224];
cohesin: [-0.064, 0.260]). This means the apparent advantage for
persistence, while a real and encouraging direction, cannot yet be
distinguished from chance at n=107 domains from 8 regions on 2
chromosomes. This is not a negative result -- it is the honest,
correct read of a promising but underpowered first test. More regions
and more chromosomes are needed before this becomes a claim rather
than a direction. See modal_app.py (function run_ctcf_benchmark).
