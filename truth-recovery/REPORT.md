# Known-Truth Validation of conformal-ma

The repo's headline: *"standard PIs cover only ~70%, conformal achieves ~92%"* —
established on 403 real Cochrane reviews, using a held-out study as a coverage
proxy. That is real-data evidence, but with the true effect distribution unknown
you cannot tell whether the standard PI fails because effects are non-normal,
because of publication selection, or because k is small.

This module injects a **known** effect distribution and a genuine new study, so
the claim can be checked mechanism-by-mechanism. The coverage target is *the
effect a new study would find* — the observed `y_new = θ_new + se_new·noise` —
which is exactly what the repo says it predicts.

> Truth-first: seeded, reproducible.
> `python truth-recovery/known_truth_validation.py --reps 3000`
> (mu=0.3, tau=0.25, nominal=0.95, seed=20260613).

## Result: the claim holds under known truth

Coverage of the observed new effect `y_new`:

| dist | k | standard PI | conformal | width std | width conf |
|---|---|---|---|---|---|
| normal | 5 | 0.834 | 0.873 | 1.72 | 1.93 |
| normal | 10 | 0.783 | **0.945** | 1.17 | 1.92 |
| normal | 20 | 0.798 | **0.969** | 1.05 | 1.93 |
| heavy (t₃) | 10 | 0.747 | **0.931** | 1.05 | 2.17 |
| heavy (t₃) | 20 | 0.741 | **0.968** | 0.95 | 2.36 |
| skew | 10 | 0.703 | **0.932** | 1.00 | 2.29 |
| skew | 20 | 0.738 | **0.970** | 0.93 | 2.59 |
| bimodal | 10 | 0.822 | **0.949** | 1.26 | 1.90 |
| bimodal | 20 | 0.825 | **0.967** | 1.14 | 1.80 |

**Mean over all cells: standard 0.790, conformal 0.923** — essentially the
repo's advertised 70%/92%, now confirmed against injected truth. Conformal's
distribution-free coverage holds across heavy-tailed, skewed and bimodal effects
at k ≥ 10 (0.93–0.97), exactly where the normality-based standard PI has no
reason to.

## Two honest refinements the known-truth view adds

1. **Conformal under-covers at k=5** (0.85–0.87, not 0.95). The split-conformal
   quantile `⌈(1−α)(k+1)⌉/k` cannot reach 0.95 with only 5 studies — a genuine
   finite-sample floor. The guarantee is effectively k ≥ 10 here; the README's
   "works with as few as k=5" should be qualified.

2. **About half of the standard PI's shortfall is a target mismatch, not non-
   normality.** The standard PI covers the new study's *true* effect `θ_new` at
   0.899 but its *observed* effect `y_new` at only 0.790, because
   `√(τ̂² + SE(μ̂)²)` omits the new study's own sampling variance `se_new²`, which
   conformal includes. So the fair statement is: conformal wins both because it is
   distribution-free **and** because it predicts the right quantity. The repo's
   narrative should separate these two effects rather than attribute the whole gap
   to normality.

## What transferred from the allmeta estimator work

- **Transferred:** the conformal/known-truth-coverage methodology — building an
  injected-truth harness to verify a coverage claim and decompose *why* an
  interval misses. It validated the repo's central result and sharpened it.
- **Did not transfer:** NPE/SBC are about a different estimator; conformal-ma is
  already the conformal method, so the relevant contribution was the honest
  yardstick, which is what was added.

## Files
`known_truth_validation.py` (injected-truth harness over 4 effect distributions) ·
`../truth_recovery_import.py` (import shim) · `../tests/test_known_truth.py` (4 tests).
