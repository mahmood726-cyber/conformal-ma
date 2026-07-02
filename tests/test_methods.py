"""Unit tests for the estimator zoo in methods.py.

These pin the structural contract every estimator must satisfy (finite point
estimate, ordered CI, correct multiplier / distribution), cross-check the DL
implementation against the audited conformal_core, and lock the reproducibility
of the bootstrap/stochastic candidates under a fixed numpy Generator.

They deliberately do NOT assert specific point-estimate *values* for the
candidate estimators (those are empirical claims graded by the honest harness),
only relationships that follow from the estimators' definitions.
"""

import math
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import stats as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import conformal_core as C  # noqa: E402
import methods as M  # noqa: E402


# A moderate, heterogeneous fixture (k=7) reused across tests.
YI = np.array([0.10, 0.22, 0.35, 0.05, 0.28, 0.15, 0.19])
VI = np.array([0.010, 0.020, 0.015, 0.030, 0.012, 0.020, 0.018])

# Estimators that are deterministic given (yi, vi).
_DETERMINISTIC = ("dl", "reml", "hksj", "pet_peese", "egger_gated")


@pytest.mark.parametrize("name", _DETERMINISTIC)
def test_deterministic_methods_contract(name):
    out = getattr(M, name)(YI, VI)
    assert math.isfinite(out["mu"])
    assert math.isfinite(out["ci_lo"]) and math.isfinite(out["ci_hi"])
    assert out["ci_lo"] < out["ci_hi"]
    assert out["method"]


def test_dl_matches_conformal_core():
    """methods.dl and conformal_core.dl_pool must agree on mu and tau2 -- they are
    two implementations of the same DerSimonian-Laird fit and are used
    interchangeably across the harness."""
    d = M.dl(YI, VI)
    f = C.dl_pool(YI, VI)
    assert d["mu"] == pytest.approx(f["mu"], abs=1e-12)
    assert d["tau2"] == pytest.approx(f["tau2"], abs=1e-12)
    # DL normal-theory CI half-width is z_{0.975} * se.
    z = sp.norm.ppf(0.975)
    assert (d["ci_hi"] - d["ci_lo"]) / 2 == pytest.approx(z * d["se"], rel=1e-12)


def test_dl_and_reml_ci_symmetric_about_mu():
    for out in (M.dl(YI, VI), M.reml(YI, VI)):
        assert (out["ci_lo"] + out["ci_hi"]) / 2 == pytest.approx(out["mu"], abs=1e-12)


def test_reml_tau2_nonnegative_and_converges():
    # Homogeneous data -> tau2 should collapse toward 0.
    yi = np.array([0.20, 0.20, 0.20, 0.20, 0.20])
    vi = np.array([0.01, 0.02, 0.015, 0.03, 0.012])
    assert M._reml_tau2(yi, vi) == pytest.approx(0.0, abs=1e-8)
    # tau2 is a variance -> never negative for any input.
    assert M._reml_tau2(YI, VI) >= 0.0
    # Clearly heterogeneous data (spread >> sampling variance) -> strictly positive.
    het_yi = np.array([0.0, 0.5, 1.0, -0.3, 0.8, 0.4])
    het_vi = np.array([0.01] * 6)
    assert M._reml_tau2(het_yi, het_vi) > 0.0


def test_hksj_multiplier_is_t_k_minus_1():
    """HKSJ CI half-width / se_hk must equal t_{k-1, 0.975}."""
    k = len(YI)
    out = M.hksj(YI, VI)
    mult = ((out["ci_hi"] - out["ci_lo"]) / 2) / out["se"]
    assert mult == pytest.approx(sp.t.ppf(0.975, k - 1), rel=1e-12)


def test_hksj_floor_widens_when_q_below_one():
    """On near-homogeneous data the Knapp-Hartung q < 1; the floor=True variant
    clamps q>=1, so its se (hence CI) is at least as wide as the unfloored HKSJ."""
    yi = np.array([0.200, 0.205, 0.198, 0.202, 0.201, 0.199])
    vi = np.array([0.02] * 6)
    plain = M.hksj(yi, vi)
    floored = M.hksj(yi, vi, floor=True)
    assert floored["se"] >= plain["se"]
    assert floored["method"] == "HKSJ-floor"
    assert (floored["ci_hi"] - floored["ci_lo"]) >= (plain["ci_hi"] - plain["ci_lo"])


def test_egger_test_symmetric_data_small_intercept():
    """With no funnel asymmetry the Egger intercept t should be modest and the
    two-sided p large."""
    rng = np.random.default_rng(11)
    se = rng.uniform(0.1, 0.5, size=40)
    yi = rng.normal(0.2, se)  # effect independent of se -> no small-study effect
    t, p = M.egger_test(yi, se ** 2)
    assert math.isfinite(t)
    assert 0.0 <= p <= 1.0
    assert p > 0.05  # should not flag asymmetry on clean data


def test_egger_test_degenerate_dof():
    # k=2 -> dof=0 -> documented (0.0, 1.0) fallback, never a crash.
    assert M.egger_test([0.1, 0.2], [0.01, 0.02]) == (0.0, 1.0)


def test_egger_gated_reports_gate_flag():
    out = M.egger_gated(YI, VI)
    assert out["method"] == "Egger-gated"
    assert isinstance(out["gate_fired"], bool)


def test_pet_peese_switches_under_strong_asymmetry():
    """Inject a strong small-study effect (imprecise studies biased upward);
    PET-PEESE must still return an ordered, finite CI."""
    rng = np.random.default_rng(3)
    se = rng.uniform(0.1, 0.6, size=25)
    yi = 0.1 + 1.5 * se + rng.normal(0, 0.05, size=25)  # bias grows with se
    out = M.pet_peese(yi, se ** 2)
    assert math.isfinite(out["mu"])
    assert out["ci_lo"] < out["ci_hi"]


def test_robust_reproducible_and_resists_outlier():
    """Same Generator seed -> identical bootstrap CI (reproducibility contract),
    and a single gross outlier drags DL further than the redescending estimator."""
    yi = np.array([0.20, 0.18, 0.22, 0.19, 0.21, 0.20, 2.50])  # last = outlier
    vi = np.array([0.02] * 7)
    r1 = M.robust_redescend(yi, vi, rng=np.random.default_rng(5))
    r2 = M.robust_redescend(yi, vi, rng=np.random.default_rng(5))
    assert r1["ci_lo"] == r2["ci_lo"] and r1["ci_hi"] == r2["ci_hi"]
    clean_mu = float(np.mean(yi[:-1]))
    assert abs(r1["mu"] - clean_mu) < abs(M.dl(yi, vi)["mu"] - clean_mu)


def test_adaptive_shrink_lambda_bounded_and_reproducible():
    a1 = M.adaptive_shrink(YI, VI, rng=np.random.default_rng(2))
    a2 = M.adaptive_shrink(YI, VI, rng=np.random.default_rng(2))
    assert 0.0 <= a1["lambda"] <= 1.0
    assert a1["ci_lo"] == a2["ci_lo"] and a1["ci_hi"] == a2["ci_hi"]
    assert a1["ci_lo"] < a1["ci_hi"]


def test_ensemble_weights_normalised_and_reproducible():
    e1 = M.ensemble_loo(YI, VI, rng=np.random.default_rng(3))
    e2 = M.ensemble_loo(YI, VI, rng=np.random.default_rng(3))
    w = e1["weights"]
    assert set(w) == set(M._ENSEMBLE_BASE)
    assert all(0.0 <= v <= 1.0 for v in w.values())
    assert sum(w.values()) == pytest.approx(1.0, abs=1e-2)  # rounded to 3 dp
    assert e1["mu"] == pytest.approx(e2["mu"], abs=1e-12)
    assert e1["ci_lo"] < e1["ci_hi"]
