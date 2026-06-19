# When a Prediction Interval Grades Its Own Homework: The Apparent Superiority of Conformal Prediction Intervals in Meta-Analysis Is a Circular-Evaluation Artifact

**Mahmood Ahmad**

Department of Cardiology, Royal Free Hospital, London, United Kingdom

ORCID: 0009-0003-7781-4478

Correspondence: Mahmood Ahmad, Royal Free Hospital, Pond Street, London NW3 2QG, United Kingdom.

---

## Abstract

**Background:** Conformal prediction promises distribution-free coverage and has been proposed as a more honest prediction interval (PI) for meta-analysis than the normal-theory (Higgins–Riley) and Hartung–Knapp–Sidik–Jonkman (HKSJ) intervals. An earlier version of this project reported that conformal PIs achieved ~92% empirical coverage across Cochrane reviews versus ~70% for the standard PI. We set out to confirm that result and instead found it to be an artifact of how coverage was measured.

**Methods:** We re-examined the original evaluation, which built a single interval from all *k* studies in a review and then counted how many of those same *k* studies fell inside it. Because the conformal interval is calibrated on exactly those leave-one-out residuals, this in-sample procedure grades conformal on the quantity it was fit to. We replaced it with two honest, non-circular evaluations: (i) a genuine leave-one-out (LOO) prediction-coverage analysis on the same 365 Cochrane reviews, in which each held-out study is predicted from the other *k*−1 studies; and (ii) a known-truth simulation in which *k* published studies are drawn from a specified random-effects model (normal, heavy-tailed, skewed, or bimodal) with optional publication selection, all three PIs are built, and coverage is checked against a genuinely held-out new study from the same truth. All three methods were given an identical predictive scale, so they differed only in the multiplier (parametric *t* versus conformal empirical quantile). The DerSimonian–Laird pooling and PI construction were validated against R `metafor` to machine precision.

**Results:** Under honest LOO on the 365 reviews (nominal 95%), mean coverage was 91.7% for the standard PI, 91.0% for HKSJ, and 85.7% for conformal — the *reverse* of the published ranking. Conformal intervals were not wider but slightly narrower (median width ratio 0.92). The claimed advantage at high heterogeneity did not appear: at I² > 75% both standard and conformal covered ~89%, and the standard PI never collapsed to the previously reported 58%. The known-truth simulation corroborated this: with correctly specified models the standard PI was at or slightly above nominal, while conformal under-covered, most severely at small *k* (finite-sample conformal calibration), approaching nominal only by *k* ≥ 10–15.

**Conclusions:** The headline finding of the original analysis was a measurement artifact. Evaluated out-of-sample, the standard random-effects prediction interval is well-calibrated to slightly conservative, and split/full conformal prediction confers no coverage advantage in realistic meta-analytic regimes and under-covers when the number of studies is small. We report this as a negative result and a cautionary tale about in-sample coverage evaluation.

**Keywords:** conformal prediction, prediction intervals, meta-analysis, coverage probability, in-sample bias, reproducibility

---

## Background

Prediction intervals (PIs) in meta-analysis estimate the range within which the true effect of a future study is expected to fall, a different and often more clinically relevant quantity than the confidence interval for the mean [1]. The standard random-effects PI assumes the study-level true effects are normally distributed and uses a *t*-multiplier (Higgins–Riley) [1, 5]; the HKSJ modification adjusts the variance and degrees of freedom [4].

Conformal prediction [2, 3] produces prediction sets with finite-sample coverage guarantees under exchangeability alone, a weaker assumption than normality. This makes it an attractive candidate for meta-analysis, where the random-effects distribution is rarely verifiable. An earlier iteration of this project implemented a leave-one-out conformal PI and reported a large coverage advantage over the standard interval. The present paper documents why that advantage was illusory and what an honest evaluation shows.

## Methods

### The flaw in the original evaluation

The original pipeline, for each review, fit one interval to all *k* studies and then defined "coverage" as the fraction of those same *k* observed effects that lay inside the interval. The conformal interval's half-width is, by construction, the (1−α) empirical quantile of the leave-one-out nonconformity scores of those very studies. Grading it on the same studies therefore guarantees near-nominal in-sample coverage independent of whether the interval would cover a genuinely new study. The parametric intervals enjoy no such alignment. The comparison was thus not a test of out-of-sample calibration but of how closely each method's width matched the empirical spread of the data it was built from — a contest conformal wins by definition.

### Core estimators and validation

We use DerSimonian–Laird random-effects pooling (μ̂, τ̂², se(μ̂)). The implementation (`src/conformal_core.py`) was validated against R `metafor::rma(method="DL")` across 60 datasets with *k* from 4 to 40: maximum absolute differences were 3.3 × 10⁻¹⁶ (μ̂), 5.6 × 10⁻¹⁷ (se), and 8.7 × 10⁻¹⁸ (τ̂²).

### A fair head-to-head

To isolate the one thing under test — the distributional assumption — all three PIs were built around an identical predictive scale, σ_pred = √(τ̂² + se(μ̂)² + se_new²), with se_new the median observed standard error (methods never see the future study's precision). The methods then differ only in the multiplier applied to σ_pred:

- **standard**: *t*₍ₖ₋₂, 0.975₎;
- **HKSJ**: *t*₍ₖ₋₁, 0.975₎ with the Hartung–Knapp inflated mean variance (floored at the DL variance);
- **conformal**: the ⌈(1−α)(k+1)⌉ / k empirical quantile of the standardised leave-one-out residuals |yᵢ − μ̂₍₋ᵢ₎| / √(τ̂²₍₋ᵢ₎ + se(μ̂)²₍₋ᵢ₎ + vᵢ).

### Honest evaluation 1 — leave-one-out on real data

For each of the 365 Cochrane reviews (Pairwise70 corpus, primary outcome, *k* ≥ 5), each study *i* in turn was held out, all three PIs were built from the other *k*−1 studies, and we recorded whether the held-out observed effect yᵢ fell inside each interval. No method ever saw the point it was graded on.

### Honest evaluation 2 — known-truth simulation

A meta-analysis of *k* published studies was drawn from a known (μ, τ²) using a specified random-effects law — normal, Student-*t*₃ (heavy-tailed), skew-normal (asymmetric), or a two-component mixture (bimodal) — each standardised to mean μ and variance τ². An optional publication-selection mechanism (Vevea–Hedges step weights or Copas selection) decided which studies were observed. The three PIs were built from the published studies; a genuinely held-out new study was then drawn from the same truth (at the anticipated median precision) and coverage of its effect recorded. Each cell used 2,000 replications; seeds are fixed (`honest_coverage.py`, seed 20260619).

## Results

### Real-data leave-one-out (n = 365, nominal 95%)

| Method | Mean LOO coverage | Median |
|---|---|---|
| Standard | 0.917 | 0.938 |
| HKSJ | 0.910 | 0.923 |
| Conformal | **0.857** | 0.875 |

The median conformal-to-standard width ratio was 0.92 — conformal intervals were *narrower*, not three times wider as previously reported. Stratified by heterogeneity, the standard PI covered 0.949 / 0.895 / 0.888 at I² < 25% / 25–75% / > 75%, while conformal covered 0.841 / 0.864 / 0.888. The previously claimed pattern (conformal 0.93 vs standard 0.58 at high heterogeneity) did not occur; the standard PI was never far below nominal out-of-sample, and conformal's deficit was largest at *low* heterogeneity, where small numbers of homogeneous studies make its empirical quantile unstable.

### Known-truth simulation

Coverage of a held-out new study (y_new), averaged over the five selection
scenarios and *k* ∈ {5, 10, 15, 25}, 2,000 replications per cell (320,000
simulated meta-analyses; seed 20260619), nominal 95%:

| Random-effects law | standard | HKSJ | conformal | width ratio (conf/std) |
|---|---|---|---|---|
| normal | 0.951 | 0.946 | 0.913 | 1.01 |
| Student-*t*₃ (heavy-tailed) | 0.954 | 0.950 | 0.911 | 1.10 |
| skew-normal | 0.955 | 0.950 | 0.916 | 1.02 |
| mixture (bimodal) | 0.949 | 0.944 | 0.911 | 0.98 |

Under the correctly specified model (normal, no selection) coverage was 0.967
(standard), 0.964 (HKSJ) and 0.928 (conformal). The conformal deficit was
strongly *k*-dependent: at *k* = 5 the standard PI covered ~0.99 while conformal
covered ~0.88–0.90 at roughly 0.62× the width; the two converged by *k* ≈ 15,
and at *k* = 25 conformal became slightly conservative and wider. Across all four
random-effects laws the standard PI was at or slightly above nominal and the
conformal PI under-covered, most severely at *k* = 5 and converging toward
nominal only by *k* ≥ 15. The HKSJ interval tracked the standard interval closely. Because the predictive distribution of an observed new study is the random-effects law convolved with normal sampling error, moderate non-normality of the random effects is partially Gaussianised at realistic τ²-to-precision ratios, which limits any distribution-free advantage; a τ²-dominant sensitivity analysis is reported in the repository.

## Discussion

The apparent superiority of conformal prediction intervals in the original analysis was entirely an artifact of in-sample coverage assessment. Once coverage is measured out-of-sample — either by genuine leave-one-out on real data or against injected ground truth — the ordering reverses: the standard random-effects prediction interval is well-calibrated to slightly conservative, and conformal prediction under-covers, particularly when the number of studies is small, which is the typical meta-analytic regime.

This is consistent with what is known about finite-sample conformal prediction: with only *k* calibration residuals the (1−α) empirical quantile is poorly estimated, and the ⌈(1−α)(k+1)⌉/k correction cannot manufacture information that *k* = 5 studies do not contain. The *t*-multiplier, often criticised as an arbitrary parametric assumption, in fact encodes a useful small-sample inflation that protects coverage.

### Why this matters

The episode is a clean cautionary example of a broader failure mode: evaluating a calibrated interval on the data used to calibrate it. Any method whose width is tuned to the empirical spread of a sample will appear well-calibrated on that sample. Honest evaluation requires either held-out data or known ground truth.

### Limitations

Leave-one-out coverage on real data treats each observed study effect as the prediction target; with no ground truth, it is the best available empirical check, and the simulation supplies the known-truth complement. The simulation does not model every real-world departure (e.g., outlier or fraudulent studies), under which a robust interval could behave differently; conformal prediction may retain advantages in heavy-contamination regimes not examined here. Our conclusions concern the standard split/full conformal PI as implemented, not every possible conformal variant.

## Conclusions

Evaluated honestly out-of-sample, conformal prediction intervals provide no coverage advantage over the standard random-effects prediction interval in meta-analysis and under-cover when studies are few; the standard interval is well-calibrated to slightly conservative. The previously reported "92% versus 70%" advantage was a circular-evaluation artifact. We report this negative result, with all code and validation, in the spirit of correcting the record.

## Code and reproducibility

All analyses are reproducible from the repository: `python -m pytest tests/`, `python pipeline.py`, `python honest_coverage.py --reps 2000`, and `Rscript tools/metafor_check.R`. See `VALIDATION.md`.

## References

1. Riley RD, Higgins JPT, Deeks JJ. Interpretation of random effects meta-analyses. BMJ. 2011;342:d549.

2. Vovk V, Gammerman A, Shafer G. Algorithmic Learning in a Random World. New York: Springer; 2005.

3. Shafer G, Vovk V. A tutorial on conformal prediction. J Mach Learn Res. 2008;9:371–421.

4. IntHout J, Ioannidis JPA, Borm GF. The Hartung-Knapp-Sidik-Jonkman method for random effects meta-analysis is straightforward and considerably outperforms the standard DerSimonian-Laird method. BMC Med Res Methodol. 2014;14:25.

5. Higgins JPT, Thompson SG, Spiegelhalter DJ. A re-evaluation of random-effects meta-analysis. J R Stat Soc Ser A. 2009;172(1):137–159.
