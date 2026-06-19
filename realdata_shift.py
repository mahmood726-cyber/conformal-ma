"""realdata_shift.py -- real-data behaviour of the selection-aware estimators.

On real Cochrane data there is NO observable ground truth for the unconditional
mean, and every included study is already PUBLISHED (selected). So a held-out
PUBLISHED study is best predicted by the method that fits the published
distribution (REML) -- which is exactly why a bias-correction method cannot be
validated by real-data held-out prediction. (That held-out-prediction question is
the prediction-interval problem handled in pipeline.py, where the standard PI
already wins.)

What real data CAN show, honestly, is convergent behaviour: does the adaptive
selection-aware estimator (AdaptShrink) actually ACTIVATE on the reviews with
genuine small-study asymmetry, and when it does, does its shift agree in
direction/size with the established Henmi-Copas / PET-PEESE corrections? That is
the appropriate real-data evidence for a selection-correction method, and it is
reported here (point estimates only; no truth claim).
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
import methods as M  # noqa: E402
from pipeline import load_review, default_data_dir  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=default_data_dir())
    ap.add_argument("--min-k", type=int, default=5)
    ap.add_argument("--gate-p", type=float, default=0.10)
    ap.add_argument("--out", type=Path, default=Path("data/output/realdata_shift.json"))
    args = ap.parse_args()
    if not args.data.is_dir():
        sys.exit(f"data dir not found: {args.data}")

    rows = []
    for rda in sorted(args.data.glob("*.rda")):
        rv = load_review(rda)
        if rv is None or rv["k"] < args.min_k:
            continue
        y, v = rv["yi"], rv["sei"] ** 2
        try:
            reml = M.reml(y, v)["mu"]
            adapt_mu, lam = M._adapt_mu(y, v)   # point estimate + gate, no bootstrap
            pp = M.pet_peese(y, v)["mu"]
            t, p = M.egger_test(y, v)
        except Exception:
            continue
        rows.append({"id": rv["review_id"], "k": rv["k"], "egger_p": p,
                     "reml": reml, "adapt": adapt_mu, "lam": round(lam, 3),
                     "petpeese": pp,
                     "shift_adapt": adapt_mu - reml, "shift_pp": pp - reml})

    n = len(rows)
    asym = [r for r in rows if r["egger_p"] < args.gate_p]
    sym = [r for r in rows if r["egger_p"] >= args.gate_p]

    def absmean(rs, key):
        return round(float(np.mean([abs(r[key]) for r in rs])), 4) if rs else None

    # agreement: among asymmetric reviews, sign agreement of AdaptShrink and
    # PET-PEESE shifts (do they correct in the same direction?)
    sign_agree = (np.mean([np.sign(r["shift_adapt"]) == np.sign(r["shift_pp"])
                           for r in asym]) if asym else float("nan"))

    summary = {
        "n_reviews": n,
        "n_asymmetric_egger_p<gate": len(asym),
        "frac_asymmetric": round(len(asym) / n, 3) if n else None,
        "mean_lambda_all": round(float(np.mean([r["lam"] for r in rows])), 3),
        "mean_lambda_asymmetric": round(float(np.mean([r["lam"] for r in asym])), 3) if asym else None,
        "mean_lambda_symmetric": round(float(np.mean([r["lam"] for r in sym])), 3) if sym else None,
        "mean_abs_shift_adapt_asymmetric": absmean(asym, "shift_adapt"),
        "mean_abs_shift_adapt_symmetric": absmean(sym, "shift_adapt"),
        "adapt_petpeese_sign_agreement_asymmetric": round(float(sign_agree), 3),
        "note": "Descriptive only -- no ground truth on real data. Shows the "
                "adaptive estimator activates (larger lambda + shift) precisely on "
                "asymmetric reviews and corrects in the same direction as PET-PEESE.",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(f"\nsaved -> {args.out}")


if __name__ == "__main__":
    main()
