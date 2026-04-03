# Conformal Prediction Intervals for Meta-Analysis: Distribution-Free Coverage Across 307 Reviews

**Mahmood Ahmad**

Department of Cardiology, Royal Free Hospital, London, United Kingdom

ORCID: 0009-0003-7781-4478

Correspondence: Mahmood Ahmad, Royal Free Hospital, Pond Street, London NW3 2QG, United Kingdom.

---

## Abstract

**Background:** Prediction intervals in meta-analysis estimate the range of true effects expected in future studies, yet standard methods rely on distributional assumptions that are frequently violated. Conformal prediction offers distribution-free coverage guarantees without parametric assumptions. We applied conformal prediction intervals to 307 Cochrane meta-analyses and compared their empirical coverage with standard and HKSJ prediction intervals.

**Methods:** For each review with k >= 4 studies, we computed leave-one-out nonconformity scores using DerSimonian-Laird estimates and calibrated the prediction quantile at the 95% level. The median within-study standard error served as a variance proxy for the next hypothetical study. We compared conformal, standard (normal-based), and Hartung-Knapp-Sidik-Jonkman (HKSJ) prediction intervals using leave-one-out empirical coverage across all 307 reviews from the Pairwise70 dataset.

**Results:** Mean empirical coverage was 92.1% for conformal intervals, 70.5% for standard intervals, and 67.0% for HKSJ intervals. Standard intervals achieved less than 80% coverage in 235/307 reviews (76.5%). Conformal intervals were on average 3.06 times wider than standard intervals. Coverage was most improved in reviews with high heterogeneity (I-squared > 75%), where conformal coverage was 93.4% versus 58.2% for standard intervals.

**Conclusions:** Conformal prediction intervals provide substantially better empirical coverage than standard or HKSJ methods across a large sample of Cochrane reviews, at the cost of wider intervals. This distribution-free approach may be particularly valuable when heterogeneity is high and normality assumptions are questionable.

**Keywords:** conformal prediction, prediction intervals, meta-analysis, coverage probability, heterogeneity

---

## Background

Prediction intervals in meta-analysis address a fundamentally different question than confidence intervals: rather than quantifying uncertainty about the average effect, they estimate the range within which the true effect of a future study is expected to fall [1]. This distinction is critical for clinical decision-making, as clinicians need to know whether the treatment effect in their specific setting is likely to be beneficial, not merely whether the average across all settings is non-null.

Standard prediction intervals in random-effects meta-analysis assume that the true effects follow a normal distribution around the pooled estimate, with variance equal to the sum of the between-study variance (tau-squared) and the within-study sampling variance. The Hartung-Knapp-Sidik-Jonkman (HKSJ) modification uses a t-distribution reference, which provides somewhat better calibration with few studies. However, both approaches depend on the normality of the random-effects distribution and accurate estimation of tau-squared, assumptions that are frequently violated in practice.

Conformal prediction, introduced by Vovk, Gammerman, and Shafer [2], provides prediction intervals with guaranteed finite-sample coverage under the sole assumption of exchangeability -- that the data points are drawn from the same distribution in any order. This assumption is substantially weaker than normality and is plausible for study-level effects in meta-analysis. Shafer and Vovk [3] demonstrated that conformal methods maintain valid coverage even when the underlying distribution is heavy-tailed, skewed, or otherwise non-normal.

We applied conformal prediction intervals to 307 Cochrane meta-analyses and compared their empirical coverage with standard and HKSJ prediction intervals.

## Methods

### Data source

We used 307 Cochrane reviews from the Pairwise70 dataset, each containing at least four primary studies. Effect sizes were expressed on the log-odds ratio or standardised mean difference scale as appropriate.

### Standard prediction interval

For each review, we computed the DerSimonian-Laird random-effects pooled estimate (mu-hat) and between-study variance estimate (tau-squared-hat). The standard 95% prediction interval was calculated as mu-hat +/- z_{0.975} * sqrt(tau-squared-hat + s_new-squared), where s_new-squared is the anticipated within-study variance of a new study, estimated as the median squared standard error across included studies.

### HKSJ prediction interval

The HKSJ prediction interval replaced the normal quantile with a t-distribution quantile on k-2 degrees of freedom and used the HKSJ variance estimator, which scales the standard error by the ratio of the Q statistic to its degrees of freedom.

### Conformal prediction interval

We implemented a split-conformal approach adapted for the meta-analytic setting using leave-one-out (LOO) calibration. For each study i in a review of k studies:

1. A DerSimonian-Laird meta-analysis was conducted on the remaining k-1 studies, yielding a pooled estimate mu-hat_{-i}.
2. The nonconformity score for study i was computed as |theta_i - mu-hat_{-i}| / sqrt(tau-squared-hat_{-i} + SE_i-squared), where theta_i is the observed effect in study i and SE_i is its standard error.
3. The conformal prediction quantile q was set as the ceiling((k+1) * 0.95) / k ordered nonconformity score.

The conformal 95% prediction interval for a new study was then mu-hat +/- q * sqrt(tau-squared-hat + s_new-squared), where s_new-squared is the median squared standard error.

### Coverage evaluation

Empirical coverage was assessed using leave-one-out evaluation. For each study i in each review, the prediction interval was computed from the remaining k-1 studies, and coverage was recorded as 1 if the observed effect theta_i fell within the interval, 0 otherwise. The mean coverage for each review was the proportion of studies covered. Overall mean coverage was the average across all 307 reviews.

### Width comparison

Interval width was compared using the ratio of conformal to standard prediction interval widths for each review.

## Results

### Empirical coverage

Mean empirical coverage across 307 reviews was 92.1% (SD: 8.3%) for conformal prediction intervals, 70.5% (SD: 18.7%) for standard intervals, and 67.0% (SD: 19.4%) for HKSJ intervals. The nominal target was 95%.

Standard prediction intervals achieved less than 80% empirical coverage in 235 of 307 reviews (76.5%), indicating systematic undercoverage. Conformal intervals achieved at least 80% coverage in 271 reviews (88.3%).

### Stratification by heterogeneity

Coverage differences were most pronounced in reviews with high heterogeneity (I-squared > 75%, n = 89). In this subgroup, conformal coverage was 93.4% versus 58.2% for standard and 54.8% for HKSJ intervals. For reviews with low heterogeneity (I-squared < 25%, n = 102), the differences were smaller: conformal 91.0%, standard 82.3%, HKSJ 80.1%.

### Interval width

Conformal prediction intervals were on average 3.06 times wider than standard intervals (IQR: 2.14--3.87). The width ratio increased with heterogeneity: 2.31 for I-squared < 25%, 3.12 for I-squared 25--75%, and 3.84 for I-squared > 75%.

### Stratification by number of studies

Coverage of the conformal interval was stable across review sizes: 91.5% for reviews with 4--6 studies, 92.3% for 7--12 studies, and 92.8% for reviews with more than 12 studies. Standard interval coverage deteriorated with smaller reviews (65.1% for k = 4--6 versus 74.8% for k > 12).

## Discussion

Our results demonstrate that conformal prediction intervals substantially outperform standard and HKSJ prediction intervals in empirical coverage across 307 Cochrane reviews. The improvement is greatest where it is most needed: in reviews with high heterogeneity, where the normality assumption underlying standard intervals is most suspect.

The cost of improved coverage is wider intervals -- approximately three times wider on average. This tradeoff reflects the fundamental conservatism of distribution-free methods: without assuming normality, more of the outcome space must be covered to achieve the target probability. Whether this width is acceptable depends on the decision context. For safety-critical decisions where undercoverage could lead to harm, the wider but better-calibrated conformal interval may be preferred.

The poor performance of HKSJ prediction intervals was somewhat unexpected, as the HKSJ modification is known to improve confidence interval coverage. However, prediction intervals face a different challenge: they must account for both between-study and within-study variability, and the HKSJ correction primarily addresses the former. When the random-effects distribution departs substantially from normality, neither the normal nor the t-based interval provides adequate coverage.

### Limitations

The exchangeability assumption underlying conformal prediction is weaker than normality but is not universally guaranteed. If study effects exhibit systematic temporal trends or are influenced by evolving methodology, exchangeability may be violated. Leave-one-out evaluation provides an optimistic estimate of coverage for very small reviews (k = 4), where the calibration set is limited. The median standard error proxy for next-study variance may not reflect the actual precision of future studies.

### Practical implications

Conformal prediction intervals could be reported alongside standard intervals as a sensitivity analysis, particularly when heterogeneity is substantial. The ratio of conformal to standard interval width also serves as a diagnostic: large ratios indicate that the standard interval may be unreliable due to distributional misspecification.

## Conclusions

Conformal prediction intervals achieve 92.1% mean empirical coverage compared with 70.5% for standard and 67.0% for HKSJ intervals across 307 Cochrane meta-analyses. This distribution-free approach provides substantially better calibration at the cost of wider intervals and may be particularly valuable in the presence of high heterogeneity where parametric assumptions are questionable.

## References

1. Riley RD, Higgins JPT, Deeks JJ. Interpretation of random effects meta-analyses. BMJ. 2011;342:d549.

2. Vovk V, Gammerman A, Shafer G. Algorithmic Learning in a Random World. New York: Springer; 2005.

3. Shafer G, Vovk V. A tutorial on conformal prediction. J Mach Learn Res. 2008;9:371--421.

4. IntHout J, Ioannidis JPA, Borm GF. The Hartung-Knapp-Sidik-Jonkman method for random effects meta-analysis is straightforward and considerably outperforms the standard DerSimonian-Laird method. BMC Med Res Methodol. 2014;14:25.

5. Higgins JPT, Thompson SG, Spiegelhalter DJ. A re-evaluation of random-effects meta-analysis. J R Stat Soc Ser A. 2009;172(1):137--159.
