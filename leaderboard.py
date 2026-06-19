"""leaderboard.py -- honest head-to-head of pooled-mean estimators on known truth.

For every regime cell (selection mechanism x random-effects law x contamination
x mu x tau2 x k) we simulate `reps` meta-analyses from a KNOWN (mu, tau2),
optionally add publication selection and/or gross-outlier contamination, run each
estimator, and score against the truth:

    bias    = mean(mu_hat - mu)
    rmse    = sqrt(mean((mu_hat - mu)^2))
    cover   = P(CI covers true mu)            (target 1 - alpha)
    width   = mean CI width

A method is "calibrated" in a cell if cover in [1-alpha-0.025, 1-alpha+0.025].
Among calibrated methods, lower RMSE / narrower width is better. The honest win
condition: a method that is calibrated AND has lower RMSE (or narrower width at
matched coverage) than every baseline, averaged across the regime grid, WITHOUT
losing badly in any single regime.

Truth-first: all numbers are produced here from seeded simulation. Baselines
(DL/REML/HKSJ) are validated against R metafor; GRMA v10, Henmi-Copas and TGEP
are scored on the SAME exported datasets by the R lane and merged in.
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
import dgp            # noqa: E402
import methods as M   # noqa: E402

BASE_SEED = 20260619

# estimator name -> callable(yi, vi, alpha, rng) ; rng used only by bootstrap ones
ESTIMATORS = {
    "DL":          lambda y, v, a, r: M.dl(y, v, a),
    "REML":        lambda y, v, a, r: M.reml(y, v, a),
    "HKSJ":        lambda y, v, a, r: M.hksj(y, v, a),
    "PET-PEESE":   lambda y, v, a, r: M.pet_peese(y, v, a),
    "Robust":      lambda y, v, a, r: M.robust_redescend(y, v, a, n_boot=300, rng=r),
    "Egger-gated": lambda y, v, a, r: M.egger_gated(y, v, a),
    "AdaptShrink": lambda y, v, a, r: M.adaptive_shrink(y, v, a, n_boot=400, rng=r),
    "Ensemble":    lambda y, v, a, r: M.ensemble_loo(y, v, a, n_boot=150, rng=r),
}


def run(reps, ks, mus, tau2s, scenarios, re_dists, contams, alpha, methods,
        export_path=None, verbose=True):
    # Two independent RNG streams: one for DATA (so the simulated datasets are
    # IDENTICAL for a given grid/reps regardless of which methods are run, enabling
    # clean merging across separate method-subset runs), one for the bootstrap
    # methods. Data rng is reseeded per cell so cells are independent of method set.
    acc = {}   # cell-key -> method -> dict of running sums
    export = []  # optional: dump datasets for the R lane
    cell_seed = 0
    for re_dist in re_dists:
        for scen in scenarios:
            for contam in contams:
                for tau2 in tau2s:
                    for k in ks:
                        for mu in mus:
                            cell = f"{re_dist}|{scen}|c{contam}|t{tau2}|k{k}|m{mu}"
                            cell_seed += 1
                            data_rng = np.random.default_rng(BASE_SEED + cell_seed)
                            method_rng = np.random.default_rng(BASE_SEED + 100000 + cell_seed)
                            for m in methods:
                                acc.setdefault(cell, {}).setdefault(
                                    m, {"bias": 0.0, "sq": 0.0, "cov": 0, "w": 0.0, "n": 0})
                            for _ in range(reps):
                                y, v, info = dgp.generate(mu, tau2, k, scen, data_rng, re_dist=re_dist)
                                if contam > 0:
                                    y, v = dgp.contaminate(y, v, data_rng, contam)
                                if not (np.all(np.isfinite(y)) and np.all(v > 0)):
                                    continue
                                if export_path is not None and len(export) < 4000:
                                    export.append({"cell": cell, "mu": mu, "tau2": tau2,
                                                   "y": y.round(6).tolist(),
                                                   "v": v.round(8).tolist()})
                                for m in methods:
                                    try:
                                        res = ESTIMATORS[m](y, v, alpha, method_rng)
                                    except Exception:
                                        continue
                                    mh = res["mu"]
                                    if not np.isfinite(mh):
                                        continue
                                    a = acc[cell][m]
                                    a["bias"] += mh - mu
                                    a["sq"] += (mh - mu) ** 2
                                    a["cov"] += int(res["ci_lo"] <= mu <= res["ci_hi"])
                                    a["w"] += res["ci_hi"] - res["ci_lo"]
                                    a["n"] += 1
            if verbose:
                print(f"  done re_dist={re_dist}", flush=True)
    if export_path is not None:
        Path(export_path).write_text(json.dumps(export), encoding="utf-8")
        if verbose:
            print(f"  exported {len(export)} datasets -> {export_path}", flush=True)
    return acc


def summarise(acc, alpha):
    nominal = 1 - alpha
    band = (nominal - 0.025, nominal + 0.025)
    cells = {}
    for cell, md in acc.items():
        cells[cell] = {}
        for m, a in md.items():
            n = a["n"]
            if n == 0:
                continue
            cells[cell][m] = {
                "bias": round(a["bias"] / n, 4),
                "rmse": round((a["sq"] / n) ** 0.5, 4),
                "cover": round(a["cov"] / n, 4),
                "width": round(a["w"] / n, 4),
                "calibrated": band[0] <= a["cov"] / n <= band[1],
                "n": n,
            }
    return {"nominal": nominal, "band": list(band), "cells": cells}


def aggregate(summary):
    """Across-grid aggregate per method: mean RMSE, mean coverage, mean width,
    fraction of cells calibrated, and a 'calibrated-RMSE' (mean RMSE over the
    cells where that method is calibrated)."""
    by_method = {}
    for cell, md in summary["cells"].items():
        for m, s in md.items():
            d = by_method.setdefault(m, {"rmse": [], "cover": [], "width": [],
                                         "calib": 0, "cells": 0, "calib_rmse": []})
            d["rmse"].append(s["rmse"]); d["cover"].append(s["cover"])
            d["width"].append(s["width"]); d["cells"] += 1
            if s["calibrated"]:
                d["calib"] += 1
                d["calib_rmse"].append(s["rmse"])
    out = {}
    for m, d in by_method.items():
        out[m] = {
            "mean_rmse": round(float(np.mean(d["rmse"])), 4),
            "mean_cover": round(float(np.mean(d["cover"])), 4),
            "mean_width": round(float(np.mean(d["width"])), 4),
            "frac_calibrated": round(d["calib"] / d["cells"], 3),
            "calib_rmse": round(float(np.mean(d["calib_rmse"])), 4) if d["calib_rmse"] else None,
            "n_cells": d["cells"],
        }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=400)
    ap.add_argument("--ks", type=int, nargs="+", default=[5, 10, 20])
    ap.add_argument("--mus", type=float, nargs="+", default=[0.0, 0.3])
    ap.add_argument("--tau2s", type=float, nargs="+", default=[0.05, 0.20])
    ap.add_argument("--scenarios", nargs="+", default=["none", "step_strong", "copas_strong"])
    ap.add_argument("--re-dists", nargs="+", default=["normal", "t3", "mixture"])
    ap.add_argument("--contams", type=float, nargs="+", default=[0.0, 0.15])
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--methods", nargs="+", default=list(ESTIMATORS.keys()))
    ap.add_argument("--export", default=None, help="dump datasets JSON for the R lane")
    ap.add_argument("--out", default="sim_output/leaderboard.json")
    args = ap.parse_args()

    t0 = time.time()
    print(f"[leaderboard] reps={args.reps} methods={args.methods}", flush=True)
    acc = run(args.reps, args.ks, args.mus, args.tau2s, args.scenarios,
              args.re_dists, args.contams, args.alpha, args.methods, args.export)
    summary = summarise(acc, args.alpha)
    summary["aggregate"] = aggregate(summary)
    summary["_meta"] = {"reps": args.reps, "ks": args.ks, "mus": args.mus,
                        "tau2s": args.tau2s, "scenarios": args.scenarios,
                        "re_dists": args.re_dists, "contams": args.contams,
                        "alpha": args.alpha, "methods": args.methods,
                        "seed": BASE_SEED, "seconds": round(time.time() - t0, 1)}
    outp = Path(args.out); outp.parent.mkdir(parents=True, exist_ok=True)
    outp.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("\nAGGREGATE (across all regime cells):")
    print(f"  {'method':12s}{'mean_rmse':>10s}{'mean_cov':>10s}{'mean_w':>9s}"
          f"{'%calib':>8s}{'calRMSE':>9s}")
    for m, s in sorted(summary["aggregate"].items(), key=lambda kv: kv[1]["mean_rmse"]):
        cr = s["calib_rmse"] if s["calib_rmse"] is not None else float("nan")
        print(f"  {m:12s}{s['mean_rmse']:>10.4f}{s['mean_cover']:>10.4f}"
              f"{s['mean_width']:>9.3f}{100*s['frac_calibrated']:>7.0f}%{cr:>9.4f}")
    print(f"\n[leaderboard] {summary['_meta']['seconds']}s -> {outp}", flush=True)


if __name__ == "__main__":
    main()
