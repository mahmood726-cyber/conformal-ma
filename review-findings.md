## REVIEW CLEAN
## Code Review Audit: ConformalMA (pipeline.py)
### Date: 2026-04-03
### Summary: 0 P0, 1 P1, 3 P2

---

#### P0 -- Critical

None.

#### P1 -- Important

- **P1-1** [Robustness]: `load_review()` imports pandas inside the function (line 168). This is fine functionally but adds startup cost per call. Minor -- does not affect correctness.

#### P2 -- Minor / Enhancement

- **P2-1** [Statistics]: DerSimonian-Laird formula is correct in both LOO (line 72: `Q - (k-2)` for k-1 studies, df=k-2) and full-data (line 94: `Q - (k-1)`, df=k-1). The `C` denominator handles the standard DL correction.

- **P2-2** [Statistics]: Conformal quantile level (line 84): `ceil((1-alpha)*(k+1))/k` is the correct finite-sample conformal quantile formula from the split conformal prediction literature.

- **P2-3** [Statistics]: Standard PI uses `t_{k-2}` (line 123) -- correct for DL-based prediction intervals. HKSJ PI uses `t_{k-1}` (line 142) with HKSJ variance adjustment -- correct.

- **P2-4** [Security]: CSV output uses `csv.DictWriter` without formula injection guards. Since all values are numeric or review IDs from filenames (not user input), this is low risk. Pipeline output is not user-facing.

#### Checklist

- [x] DL formula correct (LOO and full-data)
- [x] Conformal quantile formula correct
- [x] Standard PI uses t_{k-2} -- correct
- [x] HKSJ PI uses t_{k-1} with adjusted variance -- correct
- [x] SE computation from CI: `(CI.end - CI.start) / (2 * 1.96)` -- correct
- [x] Handles k < 4 edge case (returns None)
- [x] Handles k < 3 for standard/HKSJ PI (returns None)
- [x] Division by zero guarded: `C > 0` check, `se > 0` filter
- [x] Coverage comparison uses LOO approach (valid empirical check)
- [x] Output directory created with `mkdir(parents=True, exist_ok=True)`
