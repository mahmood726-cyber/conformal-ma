Mahmood Ahmad
Tahir Heart Institute
author@example.com

Conformal Prediction Intervals for Meta-Analysis: Distribution-Free Coverage Across 307 Cochrane Reviews

Can distribution-free prediction intervals provide guaranteed coverage for the next study effect in meta-analysis without assuming normality? Conformal prediction was adapted for random-effects meta-analysis using leave-one-out nonconformity scores from DerSimonian-Laird estimates, then applied to 307 Cochrane reviews with at least four studies each, comparing coverage against standard and HKSJ prediction intervals. The method computes standardized residuals for each left-out study, takes the calibrated quantile, and projects the interval using median standard error as proxy for next-study variance. Conformal prediction intervals achieved 92.1% mean empirical coverage compared with 70.5% for standard and 67.0% for HKSJ intervals across 307 reviews. Standard intervals exhibited undercoverage in 76.5% of reviews, while conformal intervals were wider by a factor of 3.06 but maintained finite-sample guarantees. These findings suggest conventional prediction intervals systematically understate uncertainty about future study effects in heterogeneous meta-analyses. However, the limitation of wider intervals means conformal sets may be too conservative for clinical decisions requiring precise effect boundaries.

Outside Notes

Type: methods
Primary estimand: Prediction interval coverage probability
App: ConformalMA Pipeline v1.0
Data: 307 Cochrane systematic reviews (Pairwise70 dataset)
Code: https://github.com/mahmood726-cyber/conformal-ma
Version: 1.0
Validation: DRAFT

References

1. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. 2nd ed. Wiley; 2021.
2. Higgins JPT, Thompson SG, Deeks JJ, Altman DG. Measuring inconsistency in meta-analyses. BMJ. 2003;327(7414):557-560.
3. Cochrane Handbook for Systematic Reviews of Interventions. Version 6.4. Cochrane; 2023.

AI Disclosure

This work represents a compiler-generated evidence micro-publication (i.e., a structured, pipeline-based synthesis output). AI (Claude, Anthropic) was used as a constrained synthesis engine operating on structured inputs and predefined rules for infrastructure generation, not as an autonomous author. The 156-word body was written and verified by the author, who takes full responsibility for the content. This disclosure follows ICMJE recommendations (2023) that AI tools do not meet authorship criteria, COPE guidance on transparency in AI-assisted research, and WAME recommendations requiring disclosure of AI use. All analysis code, data, and versioned evidence capsules (TruthCert) are archived for independent verification.
