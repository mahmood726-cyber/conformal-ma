"""dgp.py -- Known-truth data-generating process for meta-analyses.

Vendored from the allmeta truth-recovery bench (and the Pairwise70
truth-recovery-validation module) and EXTENDED here with non-normal
random-effects laws, because the whole point of a distribution-free prediction
interval is to be tested when the random-effects distribution is NOT normal.

A "meta-analysis" is the set of PUBLISHED studies an analyst observes. Studies
are drawn from a true random-effects model; a selection rule decides which get
published. k is the number of *published* studies; the true (mu, tau2) and the
random-effects law are what an honest method must recover / predict.

Selection mechanisms (unchanged from the vendored harness)
    none, step_weak, step_strong, copas_weak, copas_strong

Random-effects laws (NEW; all standardised to mean mu and variance tau2)
    normal   : N(mu, tau2)                       -- the assumption standard PIs make
    t3       : Student-t(df=3), scaled to var tau2 -- heavy tails
    skew     : skew-normal(a=6), scaled            -- asymmetric
    mixture  : 0.5 N(-d) + 0.5 N(+d), scaled       -- bimodal / latent subgroups

Everything is driven by an explicit numpy Generator -> fully reproducible.
"""

import numpy as np
from scipy import stats

_STEP_CUTS = np.array([0.025, 0.05])
_STEP_WEIGHTS = {
    "weak":   np.array([1.0, 0.75, 0.55]),
    "strong": np.array([1.0, 0.35, 0.10]),
}
_COPAS = {
    "weak":   {"g0": -0.10, "g1": 0.12, "rho": 0.50},
    "strong": {"g0": -0.20, "g1": 0.12, "rho": 0.90},
}

SCENARIOS = ["none", "step_weak", "step_strong", "copas_weak", "copas_strong"]
RE_DISTS = ["normal", "t3", "skew", "mixture"]

# skew-normal shape used for the "skew" law
_SKEW_A = 6.0


def _draw_se(rng, k, se_lo, se_hi):
    """Log-uniform standard errors -> realistic spread of study precisions."""
    return np.exp(rng.uniform(np.log(se_lo), np.log(se_hi), size=k))


def _standardised_re(rng, size, re_dist):
    """Draw `size` values with mean 0 and variance 1 from the chosen law."""
    if re_dist == "normal":
        return rng.standard_normal(size)
    if re_dist == "t3":
        df = 3.0
        x = rng.standard_t(df, size=size)
        return x / np.sqrt(df / (df - 2.0))  # var(t_df) = df/(df-2)
    if re_dist == "skew":
        a = _SKEW_A
        delta = a / np.sqrt(1.0 + a * a)
        mean = delta * np.sqrt(2.0 / np.pi)
        var = 1.0 - 2.0 * delta * delta / np.pi
        x = stats.skewnorm.rvs(a, size=size, random_state=rng)
        return (x - mean) / np.sqrt(var)
    if re_dist == "mixture":
        d = 1.0
        comp_var = 0.25
        signs = rng.integers(0, 2, size=size) * 2 - 1  # +/-1
        x = signs * d + rng.standard_normal(size) * np.sqrt(comp_var)
        var = d * d + comp_var  # mean 0 by symmetry
        return x / np.sqrt(var)
    raise ValueError(f"unknown re_dist {re_dist!r}")


def _draw_theta(rng, mu, tau2, size, re_dist):
    return mu + np.sqrt(tau2) * _standardised_re(rng, size, re_dist)


def generate(mu, tau2, k, scenario, rng, re_dist="normal",
             se_lo=0.10, se_hi=0.70, max_factor=400):
    """Return (y, v, info) for one published meta-analysis of observed size k.

    y    : observed effect sizes (length k)
    v    : their sampling variances (se**2)
    info : dict(n_generated, k, sel_frac, degenerate, re_dist)
    """
    if scenario == "none":
        se = _draw_se(rng, k, se_lo, se_hi)
        theta = _draw_theta(rng, mu, tau2, k, re_dist)
        y = rng.normal(theta, se)
        return y, se ** 2, {"n_generated": k, "k": k, "sel_frac": 1.0,
                            "degenerate": False, "re_dist": re_dist}

    kind = "weak" if scenario.endswith("weak") else "strong"
    is_step = scenario.startswith("step")
    if is_step:
        weights = _STEP_WEIGHTS[kind]
    else:
        cp = _COPAS[kind]

    keep_y, keep_se = [], []
    n_examined = 0
    cap = max_factor * k
    while len(keep_y) < k and n_examined < cap:
        b = max(k, 64)
        se = _draw_se(rng, b, se_lo, se_hi)
        theta = _draw_theta(rng, mu, tau2, b, re_dist)
        eps = rng.normal(0.0, 1.0, size=b)
        y = theta + se * eps
        if is_step:
            p_one = stats.norm.sf(y / se)
            idx = np.searchsorted(_STEP_CUTS, p_one, side="right")
            w = weights[idx]
            published = rng.random(b) < w
        else:
            d = cp["rho"] * eps + np.sqrt(1 - cp["rho"] ** 2) * rng.normal(0, 1, size=b)
            z = cp["g0"] + cp["g1"] / se + d
            published = z > 0
        n_examined += b
        for yi, sei, pub in zip(y, se, published):
            if pub:
                keep_y.append(yi)
                keep_se.append(sei)
                if len(keep_y) >= k:
                    break

    degenerate = len(keep_y) < k
    if degenerate:
        need = k - len(keep_y)
        se = _draw_se(rng, need, se_lo, se_hi)
        theta = _draw_theta(rng, mu, tau2, need, re_dist)
        y = rng.normal(theta, se)
        keep_y.extend(list(y))
        keep_se.extend(list(se))

    yy = np.array(keep_y[:k])
    ss = np.array(keep_se[:k])
    sel_frac = k / max(1, n_examined)
    return yy, ss ** 2, {"n_generated": n_examined, "k": k,
                         "sel_frac": sel_frac, "degenerate": degenerate,
                         "re_dist": re_dist}


def draw_new_study(rng, mu, tau2, re_dist="normal", se_new=None,
                   se_lo=0.10, se_hi=0.70):
    """Draw a genuinely held-out NEW study from the same true model.

    Returns (theta_new, y_new, se_new). theta_new is the new study's true effect
    (the classic latent PI estimand); y_new is what it would observe. The new
    study is UNCONDITIONAL (not subject to publication selection) -- a prediction
    interval should cover the population of future studies, not the published
    subset.

    se_new : if given, the new study is observed at that (anticipated) precision.
             The harness passes the median observed se so that the target
             precision matches what the methods anticipate, isolating the
             distributional question from precision-forecast noise. If None, a
             fresh log-uniform se is drawn.
    """
    theta_new = float(_draw_theta(rng, mu, tau2, 1, re_dist)[0])
    if se_new is None:
        se_new = float(_draw_se(rng, 1, se_lo, se_hi)[0])
    y_new = float(rng.normal(theta_new, se_new))
    return theta_new, y_new, se_new
