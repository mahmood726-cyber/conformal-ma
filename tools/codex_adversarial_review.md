# Codex Review: conformal-ma refutation

Line references are to the current working tree under `C:\Projects\conformal-ma`.

## 1. Methodological verdict

Mostly sound, not yet referee-proof.

The original evaluation is genuinely circular if it built/calibrated an interval on all `k` studies and then counted those same studies as covered. The current files describe that correctly: `honest_coverage.py:6-11` and `pipeline.py:8-18`. That metric is a calibration/in-sample diagnostic, not evidence of future-study prediction. It especially favors conformal because its multiplier is an empirical order statistic of those same residuals.

Your two replacements remove direct leakage. The simulation builds PIs, then draws a new unseen study (`honest_coverage.py:46-76`). The real-data pipeline deletes the held-out study before fitting (`pipeline.py:128-137`). That is the right refutation target.

Sharing `sigma_pred` is not unfair for the specific question "parametric multiplier vs empirical conformal multiplier." The shared scale is explicit in `conformal_core.py:15-24` and implemented in `predictive_sd` at `conformal_core.py:99-102`. It is a controlled comparison, not a handicap.

Caveat: this is a particular normalized jackknife/LOO conformal PI, not exact full conformal. Scores use LOO fits (`conformal_core.py:163-169`), while the final interval uses the full/training fit scale (`conformal_core.py:171-174`). That is defensible and probably not the main source of undercoverage, but do not claim an exact finite-sample conformal guarantee from this implementation.

## 2. Bugs / implementation checks

`dl_pool`: no obvious DL bug. The fixed-effect weights, Q, DL denominator C, truncation at zero, and random-effects refit are standard (`conformal_core.py:72-83`). Add input validation for nonpositive/nonfinite `vi`, but callers mostly guard this.

`standard_pi`: correct for the observable `y_new` estimand because it uses `sqrt(tau2_hat + se_mu^2 + se_new^2)` (`conformal_core.py:109-118`). It is not the conventional latent `theta_new` PI unless `se_new^2` is removed.

`hksj_pi`: the q floor is implemented (`conformal_core.py:134-140`). The formula is HK-style. But with `q=1`, HKSJ can still be narrower than the standard PI because it uses `t_{k-1}` (`conformal_core.py:139`) while standard uses `t_{k-2}` (`conformal_core.py:114`). If you claim "floored HKSJ cannot narrow below DL/standard," that claim is false for this PI implementation.

`conformal_pi`: it does exclude study `i` from its own calibration score: `np.delete` is used before refitting (`conformal_core.py:163-168`). The external real-data LOO also excludes the test point before calling `all_pis` (`pipeline.py:128-137`). No leakage bug found there.

Conformal quantile: the intended order statistic is `m = ceil((1-alpha)(k+1))` (`conformal_core.py:171`). For alpha 0.05 and calibration `k < 19`, `m > k`; exact finite-sample split conformal would require an infinite interval. This code clamps to the maximum score (`conformal_core.py:171-172`), so 95% coverage is impossible at small k. Best-case iid score coverage is about `k/(k+1)`: 4/5 for LOO folds with 4 calibration studies, 5/6 for k=5, 10/11 for k=10, 15/16 for k=15. Your 0.89-0.91 conformal coverage is therefore expected, not a mysterious bug.

NumPy detail: `np.quantile(scores, m/k, method="higher")` does not always return the m-th empirical order statistic. For larger k it is often one order statistic too high, hence conservative. Use `np.sort(scores)[m-1]` when `m <= k`, and `inf` or an explicitly disclosed max-score approximation when `m > k`.

Material real-data issue: `pipeline.py` predicts held-out `yi[i]` but calls `C.all_pis(y_lo, v_lo)` without `se_new=sei[i]` (`pipeline.py:131`). That defaults to the median training SE (`conformal_core.py:184-191`). So the real LOO checks actual held-out studies at their actual SEs against a median-precision PI. Either use `se_new=sei[i]` for observed-study LOO, or state that the estimand is not conditional on the held-out study precision.

Aggregation issue: the reported real-data coverage is a mean of per-review coverages (`pipeline.py:191-194`). Thus `n=365` is reviews, not held-out studies. Report both macro-average and micro-average over all held-out studies, with binomial/clustered uncertainty intervals.

## 3. `y_new` vs `theta_new`

Targeting `y_new` is legitimate if the estimand is "what effect estimate will a future study report at anticipated precision?" That is also the natural scale for real-data LOO because the held-out data are observed estimates, not latent true effects.

But conventional meta-analysis PIs usually target the latent true effect `theta_new`. Your code acknowledges both (`dgp.py:144-164`, `honest_coverage.py:67-71`), but the comments conflict: `honest_coverage.py:57-62` says `y_new` is primary, `honest_coverage.py:120` calls `theta_new` primary, and `honest_coverage.py:132-134` switches back to `y_new`. Fix this before submission.

Do not call the `theta_new` results "classic PI" if the interval still includes `se_new^2`; that makes the theta check conservative by construction (`honest_coverage.py:57-62`).

## 4. Is small-k conformal undercoverage real?

Yes, for this finite-sample max-score conformal implementation. At 95%, when calibration k is below 19, a finite max-score interval cannot achieve exact 95% marginal coverage. This is a finite-sample conformal limitation plus an implementation choice to clamp rather than return infinity.

It is not caused by failure to exclude the test point. It is also not caused by sharing `sigma_pred`; the same scale is applied to all methods. Larger-k behavior should be analyzed separately because the NumPy quantile currently tends to be conservative there.

## 5. Referee rejection risks

Highest risk: real-data LOO uses median training SE rather than held-out SE (`pipeline.py:131`, `conformal_core.py:184-191`). That can be attacked as an estimand mismatch.

High risk: simulation training studies may be publication-selected (`dgp.py:104-124`), but the future study is explicitly unselected (`dgp.py:150-152`). If the paper claims prediction of future published Cochrane-like studies, add a sensitivity where the new study is drawn from the same selected/published distribution.

High risk: overgeneralizing from this implementation to "conformal prediction undercovers." The fair claim is narrower: this finite-sample normalized LOO conformal PI does not beat the standard t PI in these simulations/LOO analyses.

Medium risk: `n=365` is reviews and coverage is macro-averaged (`pipeline.py:186-194`). Add total held-out study count, micro-average, and confidence intervals.

Medium risk: the primary-outcome extraction is heuristic (`pipeline.py:73-107`). A referee may ask whether outcome selection, effect-scale detection, or CI-derived SEs drive the real-data result.

Medium risk: standard/HKSJ are evaluated on the widened `y_new` scale. This is coherent, but it is not the usual latent-theta PI unless explicitly framed.

Bottom line: the circularity refutation is real. The conformal "win" from the original in-sample evaluation should not be trusted. But fix the held-out-SE issue, clarify the estimand, report micro/macro uncertainty, and add a selected-future sensitivity before making a strong manuscript claim.
