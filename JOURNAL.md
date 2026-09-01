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
