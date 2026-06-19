"""conformal_core.py -- audited core for conformal meta-analytic prediction.

This module is the single source of truth for the four estimands used across the
project. It is deliberately small, importable, and unit-tested against R
(metafor) so the simulation harness and the real-data pipeline share *exactly*
the same definitions.

Design note (the honest comparison)
-----------------------------------
A prediction interval (PI) in meta-analysis estimates the range a *new* study's
effect is expected to occupy. The fair way to compare a parametric PI against a
conformal PI is to fix everything except the one thing under test -- the
distributional assumption -- and grade both on a genuinely held-out draw.

We therefore build all three PIs around the SAME predictive scale

    sigma_pred = sqrt(tau2_hat + se_mu_hat**2 + se_new**2)

and let the methods differ ONLY in the multiplier applied to it:

    standard : t_{k-2, 1-alpha/2}              (normal random-effects assumption)
    hksj     : t_{k-1, 1-alpha/2} with HKSJ-inflated mean variance
    conformal: empirical (1-alpha) quantile of standardised leave-one-out
               residuals  |y_i - mu_hat_{-i}| / sigma_pred_{-i}

That makes the head-to-head literally "parametric-t multiplier vs
distribution-free empirical quantile, all else equal". Under a correctly
specified normal model both should be ~nominal; the conformal quantile only pays
off when the residual law departs from normal (heavy tails, skew, mixtures) --
which is exactly the claim conformal prediction makes, and exactly what the
known-truth harness in honest_coverage.py measures.

Estimand
--------
y_new : the effect a *new* study would observe, drawn from the same true model
        (theta_new ~ RE(mu, tau2); y_new ~ N(theta_new, se_new**2)). This is the
        observable "what will the next study report" target, and it is what the
        conformal calibration set (observed y_i) is on the scale of. We do NOT
        peek at the future study's precision: se_new is taken as the median of
        the observed within-study standard errors.
"""

from __future__ import annotations

import math
import numpy as np
from scipy import stats as sp_stats


# ---------------------------------------------------------------------------
# DerSimonian-Laird random-effects pooling
# ---------------------------------------------------------------------------

def dl_pool(yi, vi):
    """DerSimonian-Laird random-effects fit.

    Parameters
    ----------
    yi : array_like   effect sizes (already on the analysis scale, e.g. log-OR)
    vi : array_like   within-study sampling VARIANCES (se**2)

    Returns
    -------
    dict with mu, se_mu, tau2, Q, k. Matches metafor::rma(method="DL").
    """
    yi = np.asarray(yi, float)
    vi = np.asarray(vi, float)
    k = len(yi)
    if k < 2:
        raise ValueError("dl_pool needs k >= 2")

    wi = 1.0 / vi
    sw = wi.sum()
    mu_fe = float((wi * yi).sum() / sw)
    Q = float((wi * (yi - mu_fe) ** 2).sum())
    C = float(sw - (wi ** 2).sum() / sw)
    tau2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0

    ws = 1.0 / (vi + tau2)
    sws = ws.sum()
    mu = float((ws * yi).sum() / sws)
    se_mu = float(math.sqrt(1.0 / sws))
    return {"mu": mu, "se_mu": se_mu, "tau2": tau2, "Q": Q, "k": k}


def i_squared(Q, k):
    """Higgins I^2 (%). Returns 0 when Q <= k-1."""
    if k < 2:
        return float("nan")
    if Q <= 0:
        return 0.0
    return float(max(0.0, (Q - (k - 1)) / Q) * 100.0)


# ---------------------------------------------------------------------------
# Predictive scale shared by all methods
# ---------------------------------------------------------------------------

def predictive_sd(fit, se_new):
    """sigma_pred = sqrt(tau2 + se_mu^2 + se_new^2). Predictive SD of a new
    study's observed effect about mu_hat."""
    return math.sqrt(fit["tau2"] + fit["se_mu"] ** 2 + se_new ** 2)


# ---------------------------------------------------------------------------
# Three prediction intervals -- differ ONLY in the multiplier
# ---------------------------------------------------------------------------

def standard_pi(fit, se_new, alpha=0.05):
    """Higgins-Riley normal-theory PI: mu +/- t_{k-2} * sigma_pred."""
    k = fit["k"]
    if k < 3:
        return None
    t = sp_stats.t.ppf(1 - alpha / 2, k - 2)
    sp = predictive_sd(fit, se_new)
    half = t * sp
    return {"lo": fit["mu"] - half, "hi": fit["mu"] + half,
            "width": 2 * half, "mult": float(t), "sigma_pred": sp}


def hksj_pi(yi, vi, fit, se_new, alpha=0.05):
    """HKSJ-style PI: t_{k-1} with the Hartung-Knapp inflated mean variance.

    The HK variance estimator scales 1/sum(w*) by q = sum(w*(y-mu)^2)/(k-1).
    We floor q at 1 (the standard ad-hoc HKSJ floor; without it HK can *narrow*
    the interval below DL, see advanced-stats.md). The predictive scale then
    uses the HK-inflated mean variance in place of se_mu^2.
    """
    yi = np.asarray(yi, float)
    vi = np.asarray(vi, float)
    k = fit["k"]
    if k < 3:
        return None
    ws = 1.0 / (vi + fit["tau2"])
    sws = ws.sum()
    q = float((ws * (yi - fit["mu"]) ** 2).sum() / (k - 1))
    q = max(1.0, q)  # HKSJ floor
    var_mu_hk = q / sws
    t = sp_stats.t.ppf(1 - alpha / 2, k - 1)
    sp = math.sqrt(fit["tau2"] + var_mu_hk + se_new ** 2)
    half = t * sp
    return {"lo": fit["mu"] - half, "hi": fit["mu"] + half,
            "width": 2 * half, "mult": float(t), "sigma_pred": sp}


def conformal_pi(yi, vi, fit, se_new, alpha=0.05):
    """Split/full-conformal PI via leave-one-out standardised residuals.

    For each i: refit DL on the other k-1 studies -> mu_{-i}, tau2_{-i}, se_mu_{-i}.
    Nonconformity score  s_i = |y_i - mu_{-i}| / sqrt(tau2_{-i} + se_mu_{-i}^2 + v_i).
    The conformal multiplier is the  ceil((1-alpha)(k+1))/k  empirical quantile of
    {s_i}. The interval is mu_hat +/- q * sigma_pred, with sigma_pred on the SAME
    scale as the parametric methods so only the multiplier differs.

    Returns None for k < 4 (too few studies to calibrate a 95% quantile).
    """
    yi = np.asarray(yi, float)
    vi = np.asarray(vi, float)
    k = fit["k"]
    if k < 4:
        return None

    scores = np.empty(k)
    for i in range(k):
        y_lo = np.delete(yi, i)
        v_lo = np.delete(vi, i)
        f = dl_pool(y_lo, v_lo)
        denom = math.sqrt(f["tau2"] + f["se_mu"] ** 2 + vi[i])
        scores[i] = abs(yi[i] - f["mu"]) / denom

    # Exact conformal order statistic: the m-th smallest score, m = ceil((1-a)(k+1)).
    # (np.quantile with method="higher" can land one order statistic too high for
    # larger k, making the interval needlessly conservative; we index directly.)
    # When m > k -- which for alpha=0.05 happens whenever k < 19 -- a finite 95%
    # split-conformal set is unattainable, so we clamp to the maximum score. This
    # is a disclosed approximation, NOT an exact finite-sample guarantee.
    m = math.ceil((1 - alpha) * (k + 1))
    s = np.sort(scores)
    q = float(s[m - 1]) if m <= k else float(s[-1])
    sp = predictive_sd(fit, se_new)
    half = q * sp
    return {"lo": fit["mu"] - half, "hi": fit["mu"] + half,
            "width": 2 * half, "mult": q, "sigma_pred": sp,
            "scores": scores}


# ---------------------------------------------------------------------------
# Convenience: build all three from raw (yi, vi)
# ---------------------------------------------------------------------------

def all_pis(yi, vi, alpha=0.05, se_new=None):
    """Fit DL and return {fit, standard, hksj, conformal}. se_new defaults to the
    median observed standard error (we never peek at the future study's se)."""
    yi = np.asarray(yi, float)
    vi = np.asarray(vi, float)
    fit = dl_pool(yi, vi)
    if se_new is None:
        se_new = float(np.median(np.sqrt(vi)))
    return {
        "fit": fit,
        "se_new": se_new,
        "standard": standard_pi(fit, se_new, alpha),
        "hksj": hksj_pi(yi, vi, fit, se_new, alpha),
        "conformal": conformal_pi(yi, vi, fit, se_new, alpha),
    }
