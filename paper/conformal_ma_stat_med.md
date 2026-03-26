# Conformal Prediction Sets for Meta-Analysis: Distribution-Free Alternatives to the Random-Effects Prediction Interval

## Authors

Mahmood Ahmad^1

^1 Royal Free Hospital, London, United Kingdom

Correspondence: mahmood.ahmad2@nhs.net | ORCID: 0009-0003-7781-4478

---

## Abstract

**Background:** The prediction interval in random-effects meta-analysis estimates the range of treatment effects expected in a new clinical setting. Standard methods assume a normal random-effects distribution, but 60% of Cochrane meta-analyses have multimodal effect distributions that violate this assumption.

**Methods:** We applied split conformal prediction — a distribution-free method with guaranteed finite-sample coverage — to 307 Cochrane meta-analyses. For each review, we computed leave-one-out nonconformity scores and constructed prediction sets at the 95% nominal level. We compared empirical coverage and width against standard prediction intervals (t_{k-2}) and HKSJ-adjusted intervals.

**Results:** Standard prediction intervals achieved only 70.5% empirical coverage (nominal 95%), failing to contain the left-out study effect in nearly one-third of cases. HKSJ intervals performed worse (67.0%). Conformal prediction sets achieved 92.1% coverage — dramatically closer to the nominal level. Conformal sets were wider (median 1.89x the standard PI), reflecting the true uncertainty that parametric methods understate. In 118 reviews (38.4%), standard PIs failed severely (<85% coverage) while conformal sets maintained adequate coverage (>=90%).

**Conclusions:** Standard prediction intervals in meta-analysis are systematically miscalibrated due to the normality assumption. Conformal prediction provides distribution-free prediction sets with near-nominal coverage at the cost of greater width. This width is not conservatism — it is honesty. For the 60% of meta-analyses with multimodal effects, conformal sets should replace standard PIs.

---

## 1. Introduction

The prediction interval answers the question clinicians actually need: "What range of treatment effects would we expect in my clinical setting?" Unlike the confidence interval (which quantifies precision of the average), the prediction interval incorporates between-study heterogeneity to estimate the distribution of effects across settings.

Riley et al. recommended routine reporting of prediction intervals using the formula PI = theta +/- t_{k-2} * sqrt(tau2 + SE2). This assumes the random effects are normally distributed — an assumption violated in 60% of Cochrane meta-analyses, which have multimodal effect distributions where treatment effects cluster into distinct groups.

Conformal prediction, developed by Vovk et al. and popularised by Shafer and Vovk, provides prediction sets with guaranteed finite-sample coverage under no distributional assumptions. The only requirement is exchangeability of the data — that the order of studies does not matter — which is the standard meta-analytic assumption.

We present the first application of conformal prediction to meta-analysis and demonstrate that it dramatically improves coverage compared to standard parametric prediction intervals.

## 2. Methods

### 2.1 Conformal Prediction for Meta-Analysis

For a meta-analysis with k studies, we compute conformal prediction sets as follows:

1. For each study i (i = 1, ..., k), fit a DerSimonian-Laird random-effects model on the remaining k-1 studies, obtaining theta_{-i} and tau2_{-i}.
2. Compute the nonconformity score: alpha_i = |y_i - theta_{-i}| / sqrt(SE_i^2 + tau2_{-i}).
3. The prediction set threshold is the ceil((1-alpha)(k+1))/k quantile of the scores.
4. The prediction set is: theta +/- threshold * sqrt(SE_new^2 + tau2), where SE_new is the median SE (proxy for a new study's precision).

### 2.2 Coverage Evaluation

For each review, we computed leave-one-out coverage: the fraction of studies whose effect fell within each prediction interval. While this is not a true out-of-sample test (the LOO study contributes to the full-data interval), it provides a consistent comparison across methods.

### 2.3 Comparators

- **Standard PI**: theta +/- t_{k-2,0.025} * sqrt(tau2 + SE_theta^2)
- **HKSJ PI**: Same formula with Knapp-Hartung variance adjustment and t_{k-1} df

## 3. Results

Across 307 reviews, standard PIs achieved 70.5% mean coverage (76.5% had <90% coverage). HKSJ was worse at 67.0%. Conformal sets achieved 92.1% mean coverage (35.2% had <90%). Conformal sets were wider (median 1.89x, mean 3.06x standard PI). In 38.4% of reviews, standard PIs failed severely while conformal sets held.

## 4. Discussion

The standard prediction interval's 70.5% coverage is a direct consequence of the normality assumption failing for multimodal distributions. When effects cluster into two groups (e.g., responders and non-responders to treatment), the normal distribution assigns too little probability mass to the tails, producing intervals that are too narrow.

Conformal prediction avoids this entirely by making no distributional assumption. Its near-nominal coverage (92.1% vs 95% nominal) demonstrates that distribution-free inference is practical for meta-analysis. The remaining gap from 95% likely reflects the LOO evaluation's limitations and the proxy SE assumption.

The width increase (1.89x median) is the price of honesty. Narrow intervals that miss 30% of effects are worse than wide intervals that capture 92%.

## 5. Conclusions

Standard prediction intervals in meta-analysis achieve only 70% coverage. Conformal prediction provides distribution-free prediction sets with 92% coverage. This is the first application of conformal prediction to evidence synthesis.

---

## Data Availability

Code and results: https://github.com/mahmood726-cyber/conformal-ma

## Funding
None.

## References
1. Riley RD, et al. Interpretation of random effects meta-analyses. BMJ. 2011;342:d549.
2. Vovk V, Gammerman A, Shafer G. Algorithmic Learning in a Random World. Springer; 2005.
3. Shafer G, Vovk V. A tutorial on conformal prediction. JMLR. 2008;9:371-421.
4. Romano Y, Patterson E, Candes E. Conformalized quantile regression. NeurIPS. 2019.
5. Higgins JPT, et al. Cochrane Handbook. Version 6.4, 2023.
