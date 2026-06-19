"""pipeline.py -- real-data companion to the known-truth honest harness.

Applies the standard / HKSJ / conformal prediction intervals to the Pairwise70
Cochrane corpus and measures EMPIRICAL leave-one-out prediction coverage.

What changed (and why)
----------------------
The earlier version of this file built ONE interval from all k studies and then
counted how many of those SAME k studies fell inside it. That is in-sample and
circular -- the conformal interval is calibrated on exactly those residuals, so
it scores near-perfect by construction while the parametric intervals do not.
The "92% vs 70%" headline that produced was an artifact of that circularity, not
an out-of-sample fact (see honest_coverage.py for the known-truth refutation).

This version does HONEST leave-one-out prediction coverage, which is what the
manuscript actually describes: for each study i, the interval is built from the
OTHER k-1 studies and we check whether the held-out y_i falls inside it. No
method ever sees the point it is graded on, and all three are graded identically.

There is no ground truth on real data, so LOO prediction coverage is the best
available empirical measure; the simulation in honest_coverage.py supplies the
known-truth check.

Paths are parameters, not hardcoded. Pass --data DIR / --out DIR, or set
CONFORMAL_DATA_DIR / CONFORMAL_OUT_DIR. Default data dir is the Pairwise70 clone
under the user's home if present.
"""

import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
import conformal_core as C  # noqa: E402

try:
    import pyreadr
except ImportError:  # pragma: no cover - dependency guard
    pyreadr = None

METHODS = ("standard", "hksj", "conformal")


def default_data_dir():
    env = os.environ.get("CONFORMAL_DATA_DIR")
    if env:
        return Path(env)
    home = Path(os.path.expanduser("~"))
    for cand in (home / "Pairwise70" / "data",
                 Path(r"C:\Models\Pairwise70\data")):
        if cand.is_dir():
            return cand
    return home / "Pairwise70" / "data"


# ---------------------------------------------------------------------------
# Data loading (unchanged extraction logic; loads the Cochrane primary outcome)
# ---------------------------------------------------------------------------

def load_review(rda_path):
    import pandas as pd
    result = pyreadr.read_r(str(rda_path))
    df = list(result.values())[0].copy()
    df.columns = df.columns.str.replace(' ', '.', regex=False)
    review_id = rda_path.stem.split('_')[0]

    groups = []
    for (grp, num), sub in df.groupby(['Analysis.group', 'Analysis.number']):
        has_binary = (sub['Experimental.cases'].notna() & (sub['Experimental.cases'] > 0)).any()
        groups.append({'grp': grp, 'num': num, 'k': len(sub), 'binary': has_binary})
    if not groups:
        return None
    gdf = pd.DataFrame(groups)
    binary = gdf[gdf['binary']]
    best = binary.loc[binary['k'].idxmax()] if len(binary) > 0 else gdf.loc[gdf['k'].idxmax()]
    primary = df[(df['Analysis.group'] == best['grp']) & (df['Analysis.number'] == best['num'])]

    has_binary = (primary['Experimental.cases'].notna() & (primary['Experimental.cases'] > 0)).any()
    scale = 'ratio' if has_binary else ('ratio' if (primary['Mean'].dropna() > 0).all() else 'difference')

    if scale == 'ratio':
        v = (primary['Mean'].notna() & (primary['Mean'] > 0) & primary['CI.start'].notna() &
             (primary['CI.start'] > 0) & primary['CI.end'].notna() & (primary['CI.end'] > 0))
        sub = primary[v]
        if len(sub) < 4:
            return None
        yi = np.log(sub['Mean'].values.astype(float))
        sei = (np.log(sub['CI.end'].values.astype(float)) - np.log(sub['CI.start'].values.astype(float))) / (2 * 1.96)
    else:
        v = primary['Mean'].notna() & primary['CI.start'].notna() & primary['CI.end'].notna()
        sub = primary[v]
        if len(sub) < 4:
            return None
        yi = sub['Mean'].values.astype(float)
        sei = (sub['CI.end'].values.astype(float) - sub['CI.start'].values.astype(float)) / (2 * 1.96)

    ok = (sei > 0) & np.isfinite(yi) & np.isfinite(sei)
    yi, sei = yi[ok], sei[ok]
    if len(yi) < 4:
        return None
    return {'review_id': review_id, 'yi': yi, 'sei': sei, 'k': len(yi), 'scale': scale}


# ---------------------------------------------------------------------------
# Honest leave-one-out prediction coverage
# ---------------------------------------------------------------------------

def loo_coverage(yi, sei, alpha=0.05):
    """For each study i, build each PI from the other k-1 studies and test
    whether the held-out y_i is covered. Returns per-method coverage and the
    full-data tau2 / I2 / conformal-vs-standard width ratio (for description).

    Requires k >= 5 so that each LOO fold has k-1 >= 4 (conformal needs 4)."""
    yi = np.asarray(yi, float)
    vi = np.asarray(sei, float) ** 2
    k = len(yi)
    if k < 5:
        return None

    hits = {m: 0 for m in METHODS}
    folds = {m: 0 for m in METHODS}
    for i in range(k):
        y_lo = np.delete(yi, i)
        v_lo = np.delete(vi, i)
        # Predict the held-out study at ITS OWN precision (se_new = sei[i]), so the
        # interval targets exactly the quantity we test (the observed y_i). Using
        # the median training se instead would create an estimand mismatch.
        pis = C.all_pis(y_lo, v_lo, alpha=alpha, se_new=float(np.sqrt(vi[i])))
        for m in METHODS:
            pi = pis[m]
            if pi is None:
                continue
            folds[m] += 1
            hits[m] += int(pi["lo"] <= yi[i] <= pi["hi"])
    cov = {m: (hits[m] / folds[m] if folds[m] else float("nan")) for m in METHODS}

    full = C.all_pis(yi, vi, alpha=alpha)
    wr = (full["conformal"]["width"] / full["standard"]["width"]
          if full["standard"]["width"] > 0 else float("nan"))
    return {"cov": cov, "hits": hits, "folds": folds, "tau2": full["fit"]["tau2"],
            "I2": C.i_squared(full["fit"]["Q"], k), "width_ratio_conf_std": wr}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, default=default_data_dir())
    ap.add_argument("--out", type=Path,
                    default=Path(os.environ.get("CONFORMAL_OUT_DIR",
                                                "data/output")))
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--min-k", type=int, default=5)
    args = ap.parse_args()

    if pyreadr is None:
        sys.exit("pyreadr is required to read .rda files: pip install pyreadr")
    if not args.data.is_dir():
        sys.exit(f"data dir not found: {args.data}\n"
                 f"set --data or CONFORMAL_DATA_DIR to the Pairwise70 .rda folder")
    args.out.mkdir(parents=True, exist_ok=True)

    print("Conformal Meta-Analysis -- honest leave-one-out coverage")
    print("=" * 56)
    print(f"  data: {args.data}")

    t0 = time.time()
    rda_files = sorted(args.data.glob('*.rda'))
    results = []
    micro = {m: [0, 0] for m in METHODS}  # [hits, folds] pooled over ALL held-out studies
    for rda in rda_files:
        review = load_review(rda)
        if review is None or review['k'] < args.min_k:
            continue
        r = loo_coverage(review['yi'], review['sei'], alpha=args.alpha)
        if r is None:
            continue
        for m in METHODS:
            micro[m][0] += r["hits"][m]
            micro[m][1] += r["folds"][m]
        results.append({"review_id": review['review_id'], "k": review['k'],
                        "scale": review['scale'], "tau2": round(r["tau2"], 4),
                        "I2": round(r["I2"], 1),
                        "cov_standard": round(r["cov"]["standard"], 4),
                        "cov_hksj": round(r["cov"]["hksj"], 4),
                        "cov_conformal": round(r["cov"]["conformal"], 4),
                        "width_ratio_conf_std": round(r["width_ratio_conf_std"], 3)})

    n = len(results)
    elapsed = time.time() - t0
    if n == 0:
        sys.exit("no reviews passed filtering")

    def col(name):
        return np.array([row[name] for row in results], float)
    means = {m: float(np.nanmean(col(f"cov_{m}"))) for m in METHODS}
    meds = {m: float(np.nanmedian(col(f"cov_{m}"))) for m in METHODS}
    wr = col("width_ratio_conf_std")

    def wilson(h, nn, z=1.96):
        if nn == 0:
            return (float("nan"), float("nan"))
        p = h / nn
        d = 1 + z * z / nn
        c = p + z * z / (2 * nn)
        hw = z * ((p * (1 - p) / nn + z * z / (4 * nn * nn)) ** 0.5)
        return ((c - hw) / d, (c + hw) / d)

    micro_cov = {m: (micro[m][0] / micro[m][1] if micro[m][1] else float("nan")) for m in METHODS}
    micro_ci = {m: wilson(micro[m][0], micro[m][1]) for m in METHODS}

    print(f"  reviews: {n}   held-out studies: {micro['standard'][1]}   time: {elapsed:.1f}s")
    print(f"\n  {'Method':12s}{'macro mean':>12s}{'macro med':>11s}"
          f"{'micro':>9s}{'micro 95% CI':>20s}")
    for m in METHODS:
        lo, hi = micro_ci[m]
        print(f"  {m:12s}{means[m]:>12.4f}{meds[m]:>11.4f}{micro_cov[m]:>9.4f}"
              f"   [{lo:.4f}, {hi:.4f}]")
    print(f"\n  median conformal/standard width ratio: {np.nanmedian(wr):.3f}")
    print(f"  conformal wider than standard: "
          f"{int(np.sum(wr > 1))}/{n} ({100*np.mean(wr > 1):.1f}%)")

    # I2 stratification (honest, computed -- not asserted)
    I2 = col("I2")
    strata = {"low (<25)": I2 < 25, "mod (25-75)": (I2 >= 25) & (I2 <= 75),
              "high (>75)": I2 > 75}
    strat_out = {}
    print(f"\n  LOO coverage by heterogeneity stratum:")
    for label, mask in strata.items():
        if mask.sum() == 0:
            continue
        row = {m: round(float(np.nanmean(col(f"cov_{m}")[mask])), 4) for m in METHODS}
        row["n"] = int(mask.sum())
        strat_out[label] = row
        print(f"    {label:14s} n={row['n']:4d}  "
              f"std={row['standard']:.3f} hksj={row['hksj']:.3f} conf={row['conformal']:.3f}")

    fields = list(results[0].keys())
    with open(args.out / 'conformal_loo_results.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(results)
    summary = {"n_reviews": n,
               "n_heldout_studies": micro["standard"][1], "alpha": args.alpha,
               "loo_coverage_macro_mean": {m: round(means[m], 4) for m in METHODS},
               "loo_coverage_macro_median": {m: round(meds[m], 4) for m in METHODS},
               "loo_coverage_micro": {m: round(micro_cov[m], 4) for m in METHODS},
               "loo_coverage_micro_95CI": {m: [round(micro_ci[m][0], 4), round(micro_ci[m][1], 4)]
                                           for m in METHODS},
               "median_width_ratio_conf_std": round(float(np.nanmedian(wr)), 3),
               "by_heterogeneity": strat_out,
               "elapsed_seconds": round(elapsed, 1),
               "note": "Honest out-of-sample LOO: each held-out study predicted from the "
                       "other k-1 at its own SE; no method sees its test point. macro = mean "
                       "of per-review coverages; micro = pooled over all held-out studies."}
    with open(args.out / 'conformal_loo_summary.json', 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2)
    print(f"\n  saved -> {args.out}/")


if __name__ == '__main__':
    main()
