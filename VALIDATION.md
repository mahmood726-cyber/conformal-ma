# Validation record

Truth-first: every number here was produced by the commands shown; nothing is
hand-entered. This documents (1) the cross-checks of the core math against R,
and (2) the honest refutation of the original "conformal 92% vs standard 70%"
headline.

## 1. Core pooling validated against R `metafor` (gold standard)

`src/conformal_core.dl_pool` vs `metafor::rma(yi, vi, method="DL")`, 60 datasets,
k in {4,5,8,12,20,40} (seed 12345). Script: `tools/metafor_check.R`.

| quantity | max abs diff | verdict |
|---|---|---|
| mu (pooled estimate) | 3.33e-16 | PASS |
| se(mu)               | 5.55e-17 | PASS |
| tau2 (DerSimonian-Laird) | 8.67e-18 | PASS |

metafor version 5.0.1. Threshold 1e-6 (R-validation rule). Also independently
re-derived in a from-scratch R DL implementation (mu/se/tau2/standard-PI/HKSJ all
match to <=1.2e-15) as a cross-language correctness check.

## 2. Unit tests

`python -m pytest tests/ -q` -> 14 passed. Covers DL correctness vs a manual
implementation, tau2=0 on homogeneous data, I^2, the standard/HKSJ multipliers,
the HKSJ floor, the conformal quantile level, conformal's marginal-coverage
property, DGP reproducibility, and each random-effects law's mean/variance.

## 3. The original headline was a circular-evaluation artifact

The committed pipeline built ONE interval from all k studies, then counted how
many of those SAME k studies fell inside it. Conformal is calibrated on exactly
those leave-one-out residuals, so it is graded on the quantity it was fit to and
scores near-perfect by construction; the parametric intervals are not. Reproducing
the committed pipeline on the current 595-file corpus gives n=365 reviews (the
committed summary said 307 -- a different/older snapshot), conformal 0.919 /
standard 0.706 / hksj 0.669, i.e. the in-sample artifact reproduces.

### 3a. Honest leave-one-out on the SAME real data (`pipeline.py`)

Each held-out study is predicted from the OTHER k-1; no method sees its test
point. n=365, nominal 95%:

| method | mean LOO coverage | median |
|---|---|---|
| standard  | 0.917 | 0.938 |
| hksj      | 0.910 | 0.923 |
| conformal | **0.857** | 0.875 |

Median conformal/standard width ratio = 0.917 (conformal NARROWER, not 3x wider).
By heterogeneity: low I^2 std 0.949 / conf 0.841; high I^2 std 0.888 / conf 0.888.

**The honest ranking is the OPPOSITE of the published claim**: out-of-sample,
the standard PI covers best and conformal under-covers. The manuscript's
"conformal 93% vs standard 58% at high heterogeneity" does not hold -- standard
never collapses to 58% out-of-sample.

### 3b. Known-truth simulation (`honest_coverage.py`)

Draw k published studies from a known (mu, tau2) + publication selection + a
random-effects law; build the three PIs; draw a genuinely held-out new study from
the same truth and check coverage. Identical predictive scale across methods, so
only the multiplier (parametric t vs conformal empirical quantile) differs.
(Headline numbers from the high-replication run are in
`sim_output/honest_coverage_main.json`; see README for the table.)

Pattern (tau2=0.05, y_new estimand): standard ~0.94-0.96, conformal ~0.89-0.91
across normal/t3/skew/mixture; conformal under-covers most at k=5 (finite-sample
conformal calibration), approaching nominal by k>=10. This corroborates the
real-data finding.

## How to reproduce
```
python -m pytest tests/ -q
python pipeline.py --data <Pairwise70/data> --out data/output
python honest_coverage.py --reps 2000
Rscript tools/metafor_check.R
```
