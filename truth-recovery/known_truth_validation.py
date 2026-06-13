"""
known_truth_validation.py -- Does conformal prediction actually deliver its
promised coverage, and does it beat the standard PI, under KNOWN TRUTH?

The repo's headline claim ("standard PIs cover only ~70%, conformal achieves
~92%") comes from 403 real Cochrane reviews, using a held-out study as a coverage
proxy. That is real-data evidence but the true effect distribution is unknown, so
you cannot separate "the standard PI is wrong because effects are non-normal"
from "...because of publication selection" from "...because k is small".

This harness injects a KNOWN effect distribution and a genuine new study, so the
claim can be checked mechanism-by-mechanism. The coverage target is the effect a
NEW study would FIND -- the observed y_new = theta_new + se_new * noise -- which
is exactly the quantity the repo says it predicts.

Effect distributions for theta_i (all with mean mu, between-study SD tau):
  normal   : theta ~ N(mu, tau)                  (standard-PI assumptions hold)
  heavy    : scaled t_3                           (heavy tails)
  skew     : centered/scaled lognormal           (right-skew)
  bimodal  : 0.5 N(mu-d, .) + 0.5 N(mu+d, .)      (multimodal -> conformal's case)

Truth-first: seeded, reproducible. `python truth-recovery/known_truth_validation.py`
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from pipeline import conformal_prediction_set, standard_prediction_interval

BASE_SEED = 20260613
DISTS = ["normal", "heavy", "skew", "bimodal"]


def _draw_theta(rng, dist, mu, tau, size):
    if dist == "normal":
        return rng.normal(mu, tau, size)
    if dist == "heavy":                       # t_3 scaled to SD tau
        t = rng.standard_t(3, size)
        return mu + tau * t / np.sqrt(3.0)
    if dist == "skew":                        # lognormal centered, scaled to SD tau
        s = 0.9
        z = rng.lognormal(0.0, s, size)
        z = (z - np.exp(s * s / 2)) / np.sqrt((np.exp(s * s) - 1) * np.exp(s * s))
        return mu + tau * z
    if dist == "bimodal":                     # two symmetric modes, total SD ~ tau
        d = tau                               # mode offset
        comp = rng.random(size) < 0.5
        spread = tau * 0.35
        return np.where(comp, mu - d, mu + d) + rng.normal(0, spread, size)
    raise ValueError(dist)


def run_dist(dist, mu, tau, k, reps, alpha, seed):
    rng = np.random.default_rng(seed)
    cov_std_ynew = cov_conf_ynew = 0
    cov_std_theta = cov_conf_theta = 0
    w_std = w_conf = 0.0
    n = 0
    for _ in range(reps):
        sei = np.exp(rng.uniform(np.log(0.10), np.log(0.70), size=k))
        theta = _draw_theta(rng, dist, mu, tau, k)
        yi = theta + sei * rng.normal(0, 1, size=k)

        conf = conformal_prediction_set(yi, sei, alpha=alpha)
        if conf is None:
            continue
        std = standard_prediction_interval(conf["theta"], conf["se_theta"],
                                           conf["tau2"], k, alpha=alpha)
        if std is None:
            continue

        # a genuine NEW study: its true effect and what it would observe
        theta_new = _draw_theta(rng, dist, mu, tau, 1)[0]
        se_new = float(np.median(sei))
        y_new = theta_new + se_new * rng.normal(0, 1)

        n += 1
        # coverage of the OBSERVED new effect (the repo's stated target)
        if std["standard_lo"] <= y_new <= std["standard_hi"]:
            cov_std_ynew += 1
        if conf["conformal_lo"] <= y_new <= conf["conformal_hi"]:
            cov_conf_ynew += 1
        # coverage of the new study's TRUE effect (decomposition)
        if std["standard_lo"] <= theta_new <= std["standard_hi"]:
            cov_std_theta += 1
        if conf["conformal_lo"] <= theta_new <= conf["conformal_hi"]:
            cov_conf_theta += 1
        w_std += std["standard_width"]
        w_conf += conf["conformal_width"]

    if n == 0:
        return None
    return {
        "n": n,
        "cov_standard_ynew": round(cov_std_ynew / n, 3),
        "cov_conformal_ynew": round(cov_conf_ynew / n, 3),
        "cov_standard_theta_new": round(cov_std_theta / n, 3),
        "cov_conformal_theta_new": round(cov_conf_theta / n, 3),
        "width_standard": round(w_std / n, 3),
        "width_conformal": round(w_conf / n, 3),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=3000)
    ap.add_argument("--ks", type=int, nargs="+", default=[5, 10, 20])
    ap.add_argument("--mu", type=float, default=0.3)
    ap.add_argument("--tau", type=float, default=0.25)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--out", default="truth-recovery/known_truth_results.json")
    args = ap.parse_args()

    t0 = time.time()
    grid = {}
    seed = BASE_SEED
    for dist in DISTS:
        grid[dist] = {}
        for k in args.ks:
            seed += 1
            grid[dist][k] = run_dist(dist, args.mu, args.tau, k, args.reps,
                                     args.alpha, seed)
    result = {"mu": args.mu, "tau": args.tau, "alpha": args.alpha,
              "reps": args.reps, "ks": args.ks, "seed": BASE_SEED,
              "seconds": round(time.time() - t0, 1), "grid": grid}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(result, f, indent=2)

    nominal = 1 - args.alpha
    print(f"\n# conformal-ma known-truth validation  (mu={args.mu}, tau={args.tau}, "
          f"nominal={nominal}, reps={args.reps}/cell, seed={BASE_SEED})\n")
    print("Coverage of the OBSERVED new study effect y_new (the repo's target):\n")
    print(f"{'dist':9s} {'k':>3s} {'standard':>9s} {'conformal':>10s}   {'w_std':>6s} {'w_conf':>6s}")
    for dist in DISTS:
        for k in args.ks:
            c = grid[dist][k]
            if not c:
                continue
            print(f"{dist:9s} {k:>3d} {c['cov_standard_ynew']:>9.3f} "
                  f"{c['cov_conformal_ynew']:>10.3f}   "
                  f"{c['width_standard']:>6.3f} {c['width_conformal']:>6.3f}")
    # means
    def m(key):
        vals = [grid[d][k][key] for d in DISTS for k in args.ks if grid[d][k]]
        return sum(vals) / len(vals)
    print(f"\nMean over all cells -> standard y_new: {m('cov_standard_ynew'):.3f}  "
          f"conformal y_new: {m('cov_conformal_ynew'):.3f}  (nominal {nominal})")
    print(f"Mean theta_new      -> standard: {m('cov_standard_theta_new'):.3f}  "
          f"conformal: {m('cov_conformal_theta_new'):.3f}")
    print(f"\nWrote {args.out}  ({result['seconds']}s)")


if __name__ == "__main__":
    main()
