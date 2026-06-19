"""honest_coverage.py -- KNOWN-TRUTH, OUT-OF-SAMPLE validation of conformal vs
parametric meta-analytic prediction intervals.

Why this file exists
--------------------
The original pipeline graded a prediction interval by checking how many of the
SAME studies used to build it fell inside it (in-sample leave-one-out). That is
circular: the conformal interval is calibrated on exactly those residuals, so it
is graded on the quantity it was fit to and looks near-perfect by construction,
while the parametric intervals are not. The reported "conformal 92% vs standard
70%" is an artifact of that circularity, not an out-of-sample fact.

This harness fixes it. For each simulated meta-analysis we:
  1. draw k PUBLISHED studies from a known (mu, tau2) + selection mechanism +
     random-effects law (normal / heavy-tailed / skewed / bimodal),
  2. build the standard, HKSJ and conformal 95% PIs from those k studies,
  3. draw a GENUINELY HELD-OUT new study (theta_new, y_new) from the SAME true
     model -- never seen by any method,
  4. record whether each PI covers theta_new (the classic PI estimand) and
     y_new (the observable next-study effect).

All three methods share an identical predictive scale and target; they differ
only in the multiplier (parametric t vs distribution-free empirical quantile).
So any coverage gap is attributable to the distributional assumption alone.

Truth-first: every number is produced by seeded simulation here; nothing is
hand-entered. Reproduce with `python honest_coverage.py --reps 400`.
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
import dgp  # noqa: E402
import conformal_core as C  # noqa: E402

BASE_SEED = 20260619
METHODS = ("standard", "hksj", "conformal")


def one_rep(mu, tau2, k, scenario, re_dist, rng, alpha):
    """Return per-method coverage of theta_new and y_new, plus widths, for one
    simulated meta-analysis. None if the meta-analysis is undefined."""
    y, v, info = dgp.generate(mu, tau2, k, scenario, rng, re_dist=re_dist)
    # guard: all-equal or degenerate variance can make DL undefined
    if not np.all(np.isfinite(y)) or not np.all(v > 0):
        return None
    pis = C.all_pis(y, v, alpha=alpha)
    if pis["conformal"] is None or pis["standard"] is None or pis["hksj"] is None:
        return None

    # The held-out study is observed at the SAME precision the methods anticipate
    # (the median observed se), so coverage reflects the distributional question
    # alone, not a precision-forecast mismatch. y_new is the primary estimand;
    # theta_new (its latent true effect) is reported as a secondary, and is
    # covered conservatively because the same y-scale interval is wider than the
    # tau2+se_mu^2 scale a theta_new-only interval would use.
    se_med = pis["se_new"]
    theta_new, y_new, _ = dgp.draw_new_study(rng, mu, tau2, re_dist=re_dist,
                                             se_new=se_med)

    rec = {"theta_new": {}, "y_new": {}, "width": {}, "mult": {}}
    for m in METHODS:
        pi = pis[m]
        rec["theta_new"][m] = int(pi["lo"] <= theta_new <= pi["hi"])
        rec["y_new"][m] = int(pi["lo"] <= y_new <= pi["hi"])
        rec["width"][m] = pi["width"]
        rec["mult"][m] = pi["mult"]
    rec["tau2_hat"] = pis["fit"]["tau2"]
    rec["I2"] = C.i_squared(pis["fit"]["Q"], k)
    return rec


def run(reps, ks, scenarios, re_dists, mus, tau2, alpha, verbose=True):
    rng = np.random.default_rng(BASE_SEED)
    rows = []
    for re_dist in re_dists:
        for scen in scenarios:
            for k in ks:
                for mu in mus:
                    for _ in range(reps):
                        rec = one_rep(mu, tau2, k, scen, re_dist, rng, alpha)
                        if rec is None:
                            continue
                        rows.append({"re_dist": re_dist, "scen": scen,
                                     "k": k, "mu": mu, **rec})
            if verbose:
                print(f"  done re_dist={re_dist:8s} scen={scen}", flush=True)
    return rows


def _agg(rows, key_fields, target):
    """Mean coverage of `target` (theta_new|y_new) per method, grouped by the
    tuple of key_fields. Returns {groupkey: {method: cov, 'n': n, ...}}."""
    groups = {}
    for r in rows:
        gk = tuple(r[f] for f in key_fields)
        groups.setdefault(gk, []).append(r)
    out = {}
    for gk, rs in groups.items():
        n = len(rs)
        entry = {"n": n}
        for m in METHODS:
            entry[m] = round(float(np.mean([r[target][m] for r in rs])), 4)
        entry["mean_width_ratio_conf_std"] = round(float(np.mean(
            [r["width"]["conformal"] / r["width"]["standard"]
             for r in rs if r["width"]["standard"] > 0])), 3)
        out["|".join(map(str, gk))] = entry
    return out


def analyse(rows, alpha):
    nominal = 1 - alpha
    out = {"nominal": nominal, "n_total": len(rows)}
    # primary estimand: theta_new (the classic meta-analysis PI target)
    out["theta_new_by_re_dist"] = _agg(rows, ["re_dist"], "theta_new")
    out["theta_new_by_re_dist_scenario"] = _agg(rows, ["re_dist", "scen"], "theta_new")
    out["theta_new_by_re_dist_k"] = _agg(rows, ["re_dist", "k"], "theta_new")
    # secondary: y_new (observable next study)
    out["y_new_by_re_dist"] = _agg(rows, ["re_dist"], "y_new")
    # the honest summary: how far is each method's mean coverage from nominal,
    # under a correctly specified model (normal, no selection) vs the rest
    correct = [r for r in rows if r["re_dist"] == "normal" and r["scen"] == "none"]
    misspec = [r for r in rows if not (r["re_dist"] == "normal" and r["scen"] == "none")]
    def mc(rs, tgt):
        return {m: round(float(np.mean([r[tgt][m] for r in rs])), 4) for m in METHODS} if rs else {}
    # headline on the PRIMARY estimand y_new (matches the shared predictive scale)
    out["headline"] = {
        "estimand": "y_new (next study at median precision)",
        "correct_model_normal_none": {"n": len(correct),
                                      "y_new": mc(correct, "y_new")},
        "misspecified_any": {"n": len(misspec),
                             "y_new": mc(misspec, "y_new")},
    }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=400)
    ap.add_argument("--ks", type=int, nargs="+", default=[5, 10, 15, 25])
    ap.add_argument("--mus", type=float, nargs="+", default=[0.0, 0.3])
    ap.add_argument("--tau2", type=float, default=0.05)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--re-dists", nargs="+", default=dgp.RE_DISTS)
    ap.add_argument("--scenarios", nargs="+", default=dgp.SCENARIOS)
    ap.add_argument("--out", default="sim_output/honest_coverage.json")
    args = ap.parse_args()

    t0 = time.time()
    print(f"[honest] reps={args.reps} ks={args.ks} re_dists={args.re_dists} "
          f"scenarios={args.scenarios} tau2={args.tau2}", flush=True)
    rows = run(args.reps, args.ks, args.scenarios, args.re_dists, args.mus,
               args.tau2, args.alpha)
    out = analyse(rows, args.alpha)
    out["_meta"] = {"reps": args.reps, "ks": args.ks, "mus": args.mus,
                    "tau2": args.tau2, "alpha": args.alpha,
                    "re_dists": args.re_dists, "scenarios": args.scenarios,
                    "seed": BASE_SEED, "seconds": round(time.time() - t0, 1)}
    outpath = Path(args.out)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    outpath.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out["headline"], indent=2))
    print("y_new coverage by random-effects law:")
    print(json.dumps(out["y_new_by_re_dist"], indent=2))
    print(f"[honest] {len(rows)} sims in {out['_meta']['seconds']}s -> {outpath}",
          flush=True)


if __name__ == "__main__":
    main()
