"""Unit tests for the audited conformal core.

These cover: DL pooling correctness on a hand-checkable dataset, the conformal
multiplier/quantile logic, the HKSJ floor, the in-sample marginal-coverage
property of conformal (a sanity property distinct from the out-of-sample
honest harness), and reproducibility of the DGP.
"""

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import conformal_core as C  # noqa: E402
import dgp  # noqa: E402


def _manual_dl(yi, vi):
    yi = np.asarray(yi, float); vi = np.asarray(vi, float); k = len(yi)
    wi = 1 / vi; sw = wi.sum()
    mu_fe = (wi * yi).sum() / sw
    Q = (wi * (yi - mu_fe) ** 2).sum()
    Cc = sw - (wi ** 2).sum() / sw
    tau2 = max(0.0, (Q - (k - 1)) / Cc)
    ws = 1 / (vi + tau2); sws = ws.sum()
    mu = (ws * yi).sum() / sws
    se = math.sqrt(1 / sws)
    return mu, se, tau2, Q


def test_dl_pool_matches_manual():
    yi = [0.10, 0.22, 0.35, 0.05, 0.28]
    vi = [0.01, 0.02, 0.015, 0.03, 0.012]
    f = C.dl_pool(yi, vi)
    mu, se, tau2, Q = _manual_dl(yi, vi)
    assert f["mu"] == pytest.approx(mu, abs=1e-12)
    assert f["se_mu"] == pytest.approx(se, abs=1e-12)
    assert f["tau2"] == pytest.approx(tau2, abs=1e-12)
    assert f["Q"] == pytest.approx(Q, abs=1e-12)


def test_dl_pool_homogeneous_tau2_zero():
    # identical effects -> Q ~ 0 -> tau2 = 0
    f = C.dl_pool([0.2, 0.2, 0.2, 0.2], [0.01, 0.02, 0.015, 0.03])
    assert f["tau2"] == pytest.approx(0.0, abs=1e-12)
    assert f["mu"] == pytest.approx(0.2, abs=1e-9)


def test_i_squared_zero_when_Q_small():
    assert C.i_squared(2.0, 5) == 0.0          # Q <= k-1 -> 0
    assert C.i_squared(40.0, 5) > 50.0          # large Q -> high I2


def test_standard_pi_symmetric_and_widens_with_t():
    f = C.dl_pool([0.1, 0.2, 0.3, 0.15, 0.25], [0.02] * 5)
    pi = C.standard_pi(f, se_new=0.1)
    mid = (pi["lo"] + pi["hi"]) / 2
    assert mid == pytest.approx(f["mu"], abs=1e-12)
    assert pi["width"] > 0
    # t_{k-2} multiplier for k=5 is t_3
    from scipy import stats as sp
    assert pi["mult"] == pytest.approx(sp.t.ppf(0.975, 3), abs=1e-12)


def test_hksj_floor_applied():
    # near-homogeneous -> q would be <1 without the floor; floored mult uses t_{k-1}
    yi = [0.20, 0.205, 0.198, 0.202, 0.201]
    vi = [0.02] * 5
    f = C.dl_pool(yi, vi)
    pi = C.hksj_pi(yi, vi, f, se_new=0.1)
    from scipy import stats as sp
    assert pi["mult"] == pytest.approx(sp.t.ppf(0.975, 4), abs=1e-12)
    assert pi["sigma_pred"] > 0


def test_conformal_needs_k_at_least_4():
    f = C.dl_pool([0.1, 0.2, 0.3], [0.02] * 3)
    assert C.conformal_pi([0.1, 0.2, 0.3], [0.02] * 3, f, se_new=0.1) is None


def test_conformal_quantile_level():
    # for k=10, alpha=0.05: ceil(0.95*11)/10 = ceil(10.45)/10 = 11/10 -> capped 1.0
    yi = list(np.linspace(-0.3, 0.4, 10))
    vi = [0.05] * 10
    f = C.dl_pool(yi, vi)
    pi = C.conformal_pi(yi, vi, f, se_new=0.2)
    assert pi is not None
    assert pi["mult"] > 0
    assert len(pi["scores"]) == 10


def test_conformal_marginal_coverage_property():
    """Distribution-free in-sample sanity: across many exchangeable datasets the
    conformal interval should cover a held-out exchangeable point at ~>= nominal.
    This is a coarse smoke check (not the honest out-of-sample harness)."""
    rng = np.random.default_rng(7)
    hits, n = 0, 600
    for _ in range(n):
        k = 12
        theta = rng.normal(0.2, math.sqrt(0.05), size=k + 1)
        se = rng.uniform(0.1, 0.5, size=k + 1)
        y = rng.normal(theta, se)
        f = C.dl_pool(y[:k], se[:k] ** 2)
        pi = C.conformal_pi(y[:k], se[:k] ** 2, f, se_new=float(np.median(se[:k])))
        # held-out observed effect at median precision
        y_new = rng.normal(theta[k], float(np.median(se[:k])))
        hits += pi["lo"] <= y_new <= pi["hi"]
    cov = hits / n
    assert 0.85 <= cov <= 0.99, cov  # near nominal, generous band for k=12


def test_dgp_reproducible():
    r1 = np.random.default_rng(123)
    r2 = np.random.default_rng(123)
    a = dgp.generate(0.3, 0.05, 10, "step_strong", r1, re_dist="t3")
    b = dgp.generate(0.3, 0.05, 10, "step_strong", r2, re_dist="t3")
    assert np.allclose(a[0], b[0]) and np.allclose(a[1], b[1])


@pytest.mark.parametrize("re_dist", dgp.RE_DISTS)
def test_re_dist_mean_variance(re_dist):
    """Each random-effects law is standardised to the requested mean/variance."""
    rng = np.random.default_rng(99)
    x = dgp._draw_theta(rng, 0.3, 0.05, 200_000, re_dist)
    assert abs(x.mean() - 0.3) < 0.01
    assert abs(x.var() - 0.05) < 0.01


def test_all_pis_smoke():
    rng = np.random.default_rng(1)
    y, v, _ = dgp.generate(0.3, 0.05, 12, "none", rng)
    out = C.all_pis(y, v)
    for m in ("standard", "hksj", "conformal"):
        assert out[m]["lo"] < out[m]["hi"]
