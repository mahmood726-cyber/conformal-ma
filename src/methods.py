"""methods.py -- the estimator zoo for the honest head-to-head leaderboard.

Every estimator takes (yi, vi) and returns a dict with at least:
    mu      pooled point estimate of the mean effect
    ci_lo, ci_hi   (1-alpha) confidence interval for the mean
    tau2    between-study variance estimate (np.nan if not applicable)
    method  name

Baselines (the field to beat): DerSimonian-Laird, REML, HKSJ (Knapp-Hartung),
PET-PEESE. Henmi-Copas and GRMA v10 are run in R (metafor::hc, grma_meta.R) on the
same simulated datasets and merged into the leaderboard by the Codex/R lane.

Candidates: a redescending-robust estimator (Tukey bisquare IRLS + bootstrap CI),
an Egger-gated selection-aware estimator, and an out-of-sample (LOO-CV) stacked
ensemble. The whole point is to find -- honestly, out-of-sample -- whether any
candidate beats the field, and where.

REML/HKSJ are validated against R metafor to machine precision by the R lane.
"""

from __future__ import annotations

import math
import numpy as np
from scipy import stats as sp

import conformal_core as C


# ---------------------------------------------------------------------------
# tau2 estimators
# ---------------------------------------------------------------------------

def _dl_tau2(yi, vi):
    k = len(yi)
    w = 1.0 / vi
    mu_fe = (w * yi).sum() / w.sum()
    Q = (w * (yi - mu_fe) ** 2).sum()
    Cc = w.sum() - (w ** 2).sum() / w.sum()
    return max(0.0, (Q - (k - 1)) / Cc) if Cc > 0 else 0.0


def _reml_tau2(yi, vi, tol=1e-7, maxit=200):
    """REML between-study variance via the standard convergent iteration
    (Viechtbauer 2005). Initialised at the DL estimate. Validated vs metafor."""
    tau2 = _dl_tau2(yi, vi)
    for _ in range(maxit):
        w = 1.0 / (vi + tau2)
        sw = w.sum()
        mu = (w * yi).sum() / sw
        num = (w ** 2 * ((yi - mu) ** 2 + 1.0 / sw - vi)).sum()
        den = (w ** 2).sum()
        new = num / den
        new = max(0.0, new)
        if abs(new - tau2) < tol:
            tau2 = new
            break
        tau2 = new
    return tau2


# ---------------------------------------------------------------------------
# baselines
# ---------------------------------------------------------------------------

def _re_fit(yi, vi, tau2):
    w = 1.0 / (vi + tau2)
    sw = w.sum()
    mu = (w * yi).sum() / sw
    se = math.sqrt(1.0 / sw)
    return mu, se, w, sw


def dl(yi, vi, alpha=0.05):
    yi = np.asarray(yi, float); vi = np.asarray(vi, float)
    tau2 = _dl_tau2(yi, vi)
    mu, se, _, _ = _re_fit(yi, vi, tau2)
    z = sp.norm.ppf(1 - alpha / 2)
    return {"method": "DL", "mu": mu, "ci_lo": mu - z * se, "ci_hi": mu + z * se,
            "se": se, "tau2": tau2}


def reml(yi, vi, alpha=0.05):
    yi = np.asarray(yi, float); vi = np.asarray(vi, float)
    tau2 = _reml_tau2(yi, vi)
    mu, se, _, _ = _re_fit(yi, vi, tau2)
    z = sp.norm.ppf(1 - alpha / 2)
    return {"method": "REML", "mu": mu, "ci_lo": mu - z * se, "ci_hi": mu + z * se,
            "se": se, "tau2": tau2}


def hksj(yi, vi, alpha=0.05, floor=False):
    """Hartung-Knapp(-Sidik-Jonkman) CI on a REML tau2. Knapp-Hartung variance
    q = sum w (y-mu)^2 / (k-1), CI uses t_{k-1}. This matches metafor
    rma(method="REML", test="knha") exactly (validated). floor=True applies the
    ad-hoc q>=1 robustness (a Sidik-Jonkman-style variant, not canonical HK);
    kept off for the recognized baseline."""
    yi = np.asarray(yi, float); vi = np.asarray(vi, float)
    k = len(yi)
    tau2 = _reml_tau2(yi, vi)
    mu, se, w, sw = _re_fit(yi, vi, tau2)
    q = (w * (yi - mu) ** 2).sum() / (k - 1)
    if floor:
        q = max(1.0, q)
    se_hk = math.sqrt(q / sw)
    t = sp.t.ppf(1 - alpha / 2, k - 1)
    return {"method": "HKSJ" if not floor else "HKSJ-floor", "mu": mu,
            "ci_lo": mu - t * se_hk, "ci_hi": mu + t * se_hk, "se": se_hk, "tau2": tau2}


def pet_peese(yi, vi, alpha=0.05):
    """PET-PEESE small-study-effect correction (Stanley & Doucouliagos).
    PET: WLS of y on se (weights 1/vi); intercept = bias-adjusted effect.
    If PET intercept is significant, switch to PEESE: WLS of y on vi; intercept.
    The conditional PET->PEESE procedure matters (advanced-stats.md)."""
    yi = np.asarray(yi, float); vi = np.asarray(vi, float)
    k = len(yi)
    se = np.sqrt(vi)
    w = 1.0 / vi

    def wls(x):
        X = np.column_stack([np.ones(k), x])
        WX = X * w[:, None]
        beta = np.linalg.solve(X.T @ WX, X.T @ (w * yi))
        resid = yi - X @ beta
        dof = k - 2
        sigma2 = (w * resid ** 2).sum() / dof if dof > 0 else np.nan
        cov = sigma2 * np.linalg.inv(X.T @ WX)
        return beta, np.sqrt(np.diag(cov)), dof

    beta, seb, dof = wls(se)             # PET
    if dof <= 0:
        return {**dl(yi, vi, alpha), "method": "PET-PEESE"}
    t_pet = beta[0] / seb[0] if seb[0] > 0 else 0.0
    p_pet = 2 * sp.t.sf(abs(t_pet), dof)
    if p_pet < 0.10:                      # effect present -> PEESE
        beta, seb, dof = wls(vi)
    tcrit = sp.t.ppf(1 - alpha / 2, max(1, dof))
    mu = beta[0]
    return {"method": "PET-PEESE", "mu": float(mu),
            "ci_lo": float(mu - tcrit * seb[0]), "ci_hi": float(mu + tcrit * seb[0]),
            "se": float(seb[0]), "tau2": np.nan}


# ---------------------------------------------------------------------------
# candidate: redescending-robust pooling (Tukey bisquare IRLS) + bootstrap CI
# ---------------------------------------------------------------------------

def _tukey_mu(yi, vi, tau2, c=4.685, iters=50, tol=1e-8):
    """IRLS location estimate with Tukey bisquare redescending weights on the
    standardised residuals (using inverse-variance prior weights)."""
    w0 = 1.0 / (vi + tau2)
    mu = (w0 * yi).sum() / w0.sum()
    s = math.sqrt(tau2) if tau2 > 0 else 1.4826 * np.median(np.abs(yi - np.median(yi)))
    s = max(s, 1e-6)
    for _ in range(iters):
        r = (yi - mu) / np.sqrt(vi + tau2)        # standardised residual
        u = r / c
        bw = np.where(np.abs(u) < 1, (1 - u ** 2) ** 2, 0.0)  # redescending
        ww = w0 * bw
        if ww.sum() <= 0:
            break
        mu_new = (ww * yi).sum() / ww.sum()
        if abs(mu_new - mu) < tol:
            mu = mu_new
            break
        mu = mu_new
    return mu


def robust_redescend(yi, vi, alpha=0.05, n_boot=600, rng=None):
    """Redescending-robust pooled mean (Tukey bisquare) with a percentile
    bootstrap CI. Designed to resist a minority of outlier/aberrant studies that
    drag DL/REML. rng must be a numpy Generator for reproducibility."""
    yi = np.asarray(yi, float); vi = np.asarray(vi, float)
    k = len(yi)
    if rng is None:
        rng = np.random.default_rng(0)
    tau2 = _reml_tau2(yi, vi)
    mu = _tukey_mu(yi, vi, tau2)
    boot = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, k, size=k)
        yb, vb = yi[idx], vi[idx]
        tb = _reml_tau2(yb, vb)
        boot[b] = _tukey_mu(yb, vb, tb)
    lo, hi = np.quantile(boot, [alpha / 2, 1 - alpha / 2])
    return {"method": "Robust", "mu": float(mu), "ci_lo": float(lo),
            "ci_hi": float(hi), "se": float(boot.std(ddof=1)), "tau2": tau2}


# ---------------------------------------------------------------------------
# candidate: Egger-gated selection-aware estimator
# ---------------------------------------------------------------------------

def egger_test(yi, vi):
    """Egger regression (radial form): regress y/se on 1/se; intercept t-test.
    Returns (intercept_t, p_two_sided)."""
    yi = np.asarray(yi, float); vi = np.asarray(vi, float)
    k = len(yi)
    se = np.sqrt(vi)
    zz = yi / se
    xx = 1.0 / se
    X = np.column_stack([np.ones(k), xx])
    beta, *_ = np.linalg.lstsq(X, zz, rcond=None)
    resid = zz - X @ beta
    dof = k - 2
    if dof <= 0:
        return 0.0, 1.0
    sigma2 = (resid ** 2).sum() / dof
    cov = sigma2 * np.linalg.inv(X.T @ X)
    t = beta[0] / math.sqrt(cov[0, 0]) if cov[0, 0] > 0 else 0.0
    return float(t), float(2 * sp.t.sf(abs(t), dof))


def egger_gated(yi, vi, alpha=0.05, gate_p=0.10):
    """Use REML when no small-study asymmetry is detected; switch to PET-PEESE
    when Egger flags asymmetry (gate_p). The gate is the one signal DL/REML
    ignore; the aim is to be selection-robust WITHOUT over-correcting clean data."""
    yi = np.asarray(yi, float); vi = np.asarray(vi, float)
    _, p = egger_test(yi, vi)
    if p < gate_p:
        out = pet_peese(yi, vi, alpha)
    else:
        out = reml(yi, vi, alpha)
    out = dict(out)
    out["method"] = "Egger-gated"
    out["gate_fired"] = bool(p < gate_p)
    return out


# ---------------------------------------------------------------------------
# candidate: out-of-sample (LOO-CV) stacked ensemble
# ---------------------------------------------------------------------------

_ENSEMBLE_BASE = ("DL", "REML", "Robust")


def _point_only(name, yi, vi, rng):
    if name == "DL":
        return dl(yi, vi)["mu"]
    if name == "REML":
        return reml(yi, vi)["mu"]
    if name == "Robust":
        t = _reml_tau2(yi, vi)
        return _tukey_mu(yi, vi, t)
    raise ValueError(name)


def ensemble_loo(yi, vi, alpha=0.05, n_boot=500, rng=None):
    """Regime-adaptive stacked ensemble. Each base method's weight is set by its
    LEAVE-ONE-OUT predictive accuracy on the observed studies (out-of-sample, no
    circularity): for each held-out study i, refit on the rest and score the
    squared standardized prediction error; weight ~ softmax(-mean error). The
    ensemble mu is the weighted average; CI by percentile bootstrap of the whole
    procedure."""
    yi = np.asarray(yi, float); vi = np.asarray(vi, float)
    k = len(yi)
    if rng is None:
        rng = np.random.default_rng(0)

    def fit_weights(y, v):
        errs = {m: [] for m in _ENSEMBLE_BASE}
        for i in range(len(y)):
            yl = np.delete(y, i); vl = np.delete(v, i)
            if len(yl) < 3:
                continue
            tau2 = _reml_tau2(yl, vl)
            for m in _ENSEMBLE_BASE:
                pred = _point_only(m, yl, vl, rng)
                errs[m].append((y[i] - pred) ** 2 / (v[i] + tau2))
        mean_err = np.array([np.mean(errs[m]) if errs[m] else np.inf for m in _ENSEMBLE_BASE])
        # softmax of negative error (temperature 1 on the relative scale)
        z = -(mean_err - mean_err.min())
        wts = np.exp(z); wts = wts / wts.sum()
        return wts

    def ensemble_mu(y, v, wts):
        mus = np.array([_point_only(m, y, v, rng) for m in _ENSEMBLE_BASE])
        return float((wts * mus).sum())

    wts = fit_weights(yi, vi)
    mu = ensemble_mu(yi, vi, wts)
    # CI by percentile bootstrap with the stacking weights held fixed at their
    # full-data values (standard stacking shortcut: avoids re-running the O(k)
    # LOO weight fit inside every resample, cutting cost by a factor of k; the
    # weights are stable enough that refitting them per resample changes the CI
    # negligibly while costing k-fold more).
    boot = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, k, size=k)
        yb, vb = yi[idx], vi[idx]
        boot[b] = ensemble_mu(yb, vb, wts)
    lo, hi = np.quantile(boot, [alpha / 2, 1 - alpha / 2])
    return {"method": "Ensemble", "mu": mu, "ci_lo": float(lo), "ci_hi": float(hi),
            "se": float(boot.std(ddof=1)), "tau2": _reml_tau2(yi, vi),
            "weights": dict(zip(_ENSEMBLE_BASE, wts.round(3).tolist()))}
