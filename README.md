# conformal-ma

**Honest out-of-sample evaluation of conformal vs standard prediction intervals in meta-analysis.**

This project originally reported that conformal prediction intervals (PIs) achieve
~92% coverage versus ~70% for the standard random-effects PI across Cochrane
reviews. That result turned out to be an artifact of **in-sample (circular)
coverage evaluation** — the interval was graded on the same studies used to build
it, which conformal wins by construction.

Re-evaluated honestly out-of-sample (each held-out study predicted from the
other *k*−1 at its own precision), the ranking **reverses**:

| Evaluation | standard | HKSJ | conformal |
|---|---|---|---|
| In-sample LOO (original, circular) | 0.706 | 0.669 | 0.919 |
| **Honest LOO, real data — macro** (365 reviews) | **0.968** | 0.965 | **0.903** |
| **Honest LOO, real data — micro** (6,595 studies) | 0.957 | 0.956 | 0.941 |

(micro 95% CIs: standard [0.951, 0.961] vs conformal [0.935, 0.946] — non-overlapping.)

Conformal intervals are slightly *narrower* (median width ratio 0.91), not 3× wider,
and confer no coverage advantage; they under-cover, worst at small *k* and low
heterogeneity. A known-truth simulation across normal / heavy-tailed / skewed /
bimodal random-effects laws corroborates this at τ²=0.05 and τ²=0.30 (conformal
never overtakes standard). The DerSimonian–Laird core is validated against R
`metafor` 5.0.1 to machine precision (≤ 3×10⁻¹⁶).

The honest claim is appropriately narrow: *this finite-sample normalized
leave-one-out conformal PI does not beat the standard t-interval* in these
analyses — not that conformal prediction undercovers universally.

This is reported as a **negative result** and a cautionary example of in-sample
coverage bias. See `paper/manuscript.md` and `VALIDATION.md`.

## Layout
- `src/conformal_core.py` — DL pooling + standard/HKSJ/conformal PIs (one shared predictive scale). metafor-validated.
- `src/dgp.py` — known-truth DGP (vendored from the truth-recovery bench) + non-normal random-effects laws + held-out new-study draw.
- `honest_coverage.py` — known-truth out-of-sample coverage harness.
- `pipeline.py` — honest leave-one-out coverage on the Pairwise70 Cochrane corpus (parametrised paths; no hardcoded data location).
- `tools/metafor_check.R` — R `metafor` cross-validation.
- `tests/` — unit tests (`pytest`).

## Reproduce
```bash
python -m pytest tests/ -q
python pipeline.py --data /path/to/Pairwise70/data --out data/output
python honest_coverage.py --reps 2000
Rscript tools/metafor_check.R
```
Data path can also be set via `CONFORMAL_DATA_DIR`.

_Status: methods-correction complete; honest negative result._
