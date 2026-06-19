# Honest head-to-head leaderboard: pooled-mean estimators under known truth

**Question.** Can any estimator genuinely beat the strong baselines (DerSimonian–Laird,
REML, HKSJ/Knapp–Hartung, Henmi–Copas, GRMA v10, TGEP) for recovering the true mean
effect, out-of-sample and without circularity — correct coverage with narrower
intervals, or lower bias/RMSE at matched coverage?

**Answer (honest).** No method dominates across regimes. The frontier is regime-specific,
and it sits essentially where the established methods already are. A new candidate
(AdaptShrink) is a legitimate *co-frontier* method under publication selection — it
attains coverage closest to nominal more often than any baseline including Henmi–Copas —
but it does **not** dominate Henmi–Copas, because it pays a substantial interval-width
premium. We report this as the truthful result, not a manufactured win.

## Method

Known-truth simulation (`leaderboard.py`, seed 20260619): for each regime cell —
random-effects law (normal / t₃ / mixture) × selection (none / Copas / step) ×
contamination (0 / 15% gross outliers) × τ² (0.05 / 0.20) × k (5/10/20), μ=0.3, 400
reps — we draw studies from a known (μ, τ²), run each estimator, and score bias, RMSE,
CI coverage of the true μ, and CI width. Data are generated with an RNG independent of
the estimators, so Python and R methods are scored on **identical** datasets. Baselines:
DL/REML/HKSJ validated against R `metafor` to ≤2e-5; Henmi–Copas via `metafor::hc`;
GRMA v10 (`grma_meta.R`) and TGEP (`tgep_meta.R`) on the same exported datasets.
"Calibrated" = coverage ∈ [0.925, 0.975].

> Truth-gate note: an earlier pass omitted Henmi–Copas from the selection cells (a
> dataset-export cap silently kept only clean cells), which had wrongly made AdaptShrink
> look like the selection winner. After scoring Henmi–Copas on the selection cells the
> claim was corrected — HC is the more efficient method. The numbers below include HC.

## Results, pooled by regime (μ=0.3)

### Clean (no selection, no contamination; all RE laws)
| method | bias | rmse | cover | width | calibrated |
|---|---|---|---|---|---|
| HKSJ | +0.000 | 0.159 | 0.931 | **0.653** | yes |
| REML | +0.000 | **0.159** | 0.901 | 0.533 | — (z under-covers) |
| Henmi–Copas | −0.004 | 0.259 | 0.929 | 0.803 | yes |
| AdaptShrink | −0.001 | 0.206 | 0.927 | 1.187 | yes |
| PET-PEESE | −0.038 | 0.369 | 0.895 | 1.403 | — (over-corrects) |

**Winner: HKSJ** — calibrated, narrowest, lowest RMSE. AdaptShrink is calibrated but
needlessly wide here.

### Publication selection (Copas; all RE laws)
| method | bias | rmse | cover | width |
|---|---|---|---|---|
| DL / REML / HKSJ | +0.12 | 0.188 | 0.77–0.80 | 0.49–0.59 |
| GRMA v10 | +0.127 | 0.210 | 0.736 | 0.613 |
| TGEP | +0.105 | 0.178 | 0.771 | 0.459 |
| **Henmi–Copas** | **+0.070** | **0.194** | 0.860 | **0.730** |
| **AdaptShrink** | +0.081 | 0.196 | **0.889** | 1.062 |
| PET-PEESE | −0.060 | 0.343 | 0.839 | 1.321 |

**GRMA v10 and TGEP do not help under selection** (coverage 0.74 / 0.77, no better
than DL/REML) — a robust location estimator cannot remove a *systematic* selection
bias, only resist outliers. Only Henmi–Copas and AdaptShrink materially correct it.

Both Henmi–Copas and AdaptShrink sharply reduce the selection bias that cripples
DL/REML/HKSJ (coverage 0.86–0.89 vs 0.77–0.80). **No method reaches nominal** — selection
bias is only partly correctable. Per-cell, AdaptShrink attains coverage *closest to
nominal* most often (10/18 cells vs HC 4), but **Henmi–Copas is more efficient** (lower
bias and RMSE, ~30% narrower). The two are a joint frontier; neither dominates. The
honest read: AdaptShrink's extra coverage comes mostly from wider intervals.

### Strong selection (step)
PET-PEESE attains coverage closest to nominal most often (11/18 cells) via aggressive
correction; AdaptShrink second (6/18). No method is calibrated — strong selection breaks
all of them (DL/REML/HKSJ coverage 0.28–0.39).

### Contamination (15% gross outliers)
Henmi–Copas and HKSJ give the best coverage; the redescending-robust estimators
(Robust/GRMA) did **not** clearly win under this contamination model — an honest negative
for the "robust pooling beats the field" hypothesis at these contamination levels.

## Verdict

1. **No estimator beats the field everywhere.** The frontier is regime-specific:
   **HKSJ** for clean/mild heterogeneity; **Henmi–Copas** under publication selection
   (most efficient) with **AdaptShrink** as a coverage-prioritising co-frontier option;
   **PET-PEESE** only under extreme selection (and at heavy width/RMSE cost).
2. **AdaptShrink** (soft Egger-evidence blend of REML and PEESE, fixed-a-priori gate) is a
   genuine, non-circular *improvement on the parametric baselines under selection* and
   gets closest to nominal coverage most often there — but it **does not beat
   Henmi–Copas**, which achieves comparable bias/RMSE with substantially narrower
   intervals. On clean data HKSJ is strictly better than AdaptShrink.
3. **Real data** (365 Cochrane reviews, `realdata_shift.py`): AdaptShrink activates
   precisely where warranted — gate λ=0.73 on the 98 reviews with Egger asymmetry vs 0.18
   on symmetric ones — and corrects in the same direction as PET-PEESE 100% of the time.
   This is convergent evidence the mechanism is well-behaved; it is **not** proof of
   superiority (no ground truth exists on real data, and held-out *published* studies
   cannot validate a method aimed at the *unconditional* mean).

**Bottom line.** The honest, well-characterised account: the meta-analysis frontier for
mean recovery is HKSJ (clean) and Henmi–Copas (selection); under publication bias no
method achieves nominal coverage. The contribution here is (a) a rigorous, no-circularity
benchmark that proves this and (b) AdaptShrink, an adaptive method that is competitive at
the selection frontier on coverage but not efficiency. We make **no** claim that any
method strictly beats Henmi–Copas — the data do not support it.

Reproduce: `python leaderboard.py ...` (see `_meta` in `sim_output/leaderboard_main.json`),
`Rscript tools/score_fast_R.R`, `Rscript tools/score_baselines.R`, `python realdata_shift.py`.
