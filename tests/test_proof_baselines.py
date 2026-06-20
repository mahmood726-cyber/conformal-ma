"""Oracle-locked numerical baselines — auto-generated 2026-06-20.

Each test pins a conformal-ma estimator to a value independently cross-checked
against a from-scratch scipy/closed-form oracle (engine NOT used on the
oracle side). Proof workflow: prove-methods-repos. Deterministic methods
pinned to <=1e-9; iterative to the engine stopping tol; Monte-Carlo to a
3-sigma band with a pinned seed (per the Monte-Carlo testing rule).
"""


# --- conformal_prediction_set ---
import numpy as np
import pipeline as P


def test_conformal_prediction_set_baseline():
    yi = np.array([0.20, -0.10, 0.35, 0.05, -0.25, 0.15])
    sei = np.array([0.10, 0.15, 0.12, 0.20, 0.18, 0.14])
    res = P.conformal_prediction_set(yi, sei, alpha=0.05)
    assert abs(res['theta'] - 0.09555811282017351) < 1e-9
    assert abs(res['tau2'] - 0.022860981687184877) < 1e-9
    assert abs(res['threshold'] - 2.0018731629401025) < 1e-9
    assert abs(res['conformal_lo'] - (-0.3238139214300354)) < 1e-9
    assert abs(res['conformal_hi'] - 0.5149301470703824) < 1e-9


def test_conformal_heldout_is_out_of_sample():
    # Held-out outlier must be genuinely EXCLUDED (non-circular coverage).
    yi = np.array([3.00, -0.10, 0.35, 0.05, -0.25, 0.15, 0.10])
    sei = np.array([0.10, 0.15, 0.12, 0.20, 0.18, 0.14, 0.11])
    conf = P.conformal_prediction_set(np.delete(yi, 0), np.delete(sei, 0), 0.05)
    assert not (conf['conformal_lo'] <= 3.0 <= conf['conformal_hi'])
    cov = P.heldout_interval_coverage(yi, sei)['conformal']
    assert abs(cov - 6.0 / 7.0) < 1e-9

# --- standard_prediction_interval (random-effects prediction interval) ---
def test_standard_prediction_interval_baseline():
    """Locks pipeline.standard_prediction_interval to its verified t_{k-2}
    (IntHout 2016) behaviour. NOTE: this is the t_{k-2} convention, which is
    DELIBERATELY NOT Cochrane v6.5 (t_{k-1}); see pipeline.py:146 docstring.
    Verified by independent scipy oracle (max abs diff 0.0 vs t_{k-2})."""
    import math
    import numpy as np
    from scipy import stats
    import pipeline

    theta, se_theta, tau2, alpha = 0.35, 0.08, 0.045, 0.05

    # Exact verified values for k=5 (engine == IntHout t_{k-2} oracle)
    r = pipeline.standard_prediction_interval(theta, se_theta, tau2, 5, alpha)
    assert abs(r['standard_lo'] - (-0.3715104812690635)) < 1e-9
    assert abs(r['standard_hi'] - (1.0715104812690635)) < 1e-9
    assert abs(r['standard_width'] - 1.443020962538127) < 1e-9

    # Confirm convention is t_{k-2} across a range of k via independent oracle
    pi_se = math.sqrt(tau2 + se_theta ** 2)
    for k in (4, 5, 6, 8, 10):
        eng = pipeline.standard_prediction_interval(theta, se_theta, tau2, k, alpha)
        t_km2 = stats.t.ppf(1 - alpha / 2, k - 2)
        assert abs(eng['standard_width'] - 2 * t_km2 * pi_se) < 1e-9
        # And it is NOT the Cochrane t_{k-1} interval
        t_km1 = stats.t.ppf(1 - alpha / 2, k - 1)
        assert abs(eng['standard_width'] - 2 * t_km1 * pi_se) > 1e-6

    # Guard: PI undefined for k < 3
    assert pipeline.standard_prediction_interval(theta, se_theta, tau2, 2, alpha) is None

# --- hksj_prediction_interval ---
def test_hksj_prediction_interval_baseline():
    import numpy as np
    from pipeline import hksj_prediction_interval

    yi = np.array([0.20, 0.45, -0.10, 0.33, 0.05, 0.28, 0.50])
    sei = np.array([0.15, 0.22, 0.18, 0.25, 0.30, 0.20, 0.12])
    k = len(yi)
    # Reconstruct the DL tau2 / RE theta the pipeline feeds in
    wi = 1.0 / sei**2
    sw = wi.sum()
    theta_fe = (wi * yi).sum() / sw
    Q = (wi * (yi - theta_fe) ** 2).sum()
    C = sw - (wi**2).sum() / sw
    tau2 = max(0.0, (Q - (k - 1)) / C)
    w = 1.0 / (sei**2 + tau2)
    theta = float((w * yi).sum() / w.sum())

    res = hksj_prediction_interval(yi, sei, theta, tau2, k, 0.05)

    assert abs(float(res['hksj_lo']) - (-0.13380069391068705)) < 1e-9
    assert abs(float(res['hksj_hi']) - 0.666949378681398) < 1e-9
    assert abs(float(res['hksj_width']) - 0.8007500725920851) < 1e-9
    # q* floor engages (raw < 1 -> floored to 1.0)
    assert abs(float(res['hksj_scale_raw']) - 0.8910121290934109) < 1e-9
    assert res['hksj_scale'] == 1.0
    # k < 3 is undefined
    assert hksj_prediction_interval(yi[:2], sei[:2], 0.0, 0.0, 2, 0.05) is None

# --- heldout_interval_coverage (leave-one-out out-of-sample empirical coverage of conformal / standard-PI / HKSJ-PI intervals) ---
import numpy as np
import pytest


def test_heldout_interval_coverage_baseline():
    """Locks heldout_interval_coverage to its verified LOO out-of-sample values.

    Independently cross-checked vs a from-scratch numpy/scipy reimplementation
    (exact match, maxAbsDiff=0.0) and proven non-circular: perturbing a held-out
    point leaves its own fold interval unchanged.
    """
    import pipeline as P

    yi = np.array([0.20, 0.05, 0.35, -0.10, 0.15, 0.25, 0.00])
    sei = np.array([0.10, 0.15, 0.12, 0.20, 0.11, 0.14, 0.18])

    res = P.heldout_interval_coverage(yi, sei)
    assert res['conformal'] == pytest.approx(5.0 / 7.0, abs=1e-9)
    assert res['standard'] == pytest.approx(4.0 / 7.0, abs=1e-9)
    assert res['hksj'] == pytest.approx(4.0 / 7.0, abs=1e-9)

    # k=5 subset
    res5 = P.heldout_interval_coverage(yi[:5], sei[:5])
    assert res5['conformal'] == pytest.approx(0.6, abs=1e-9)
    assert res5['standard'] == pytest.approx(0.8, abs=1e-9)
    assert res5['hksj'] == pytest.approx(0.6, abs=1e-9)

    # Non-circularity: the fold interval for the held-out study must not depend
    # on the held-out value (genuine out-of-sample LOO).
    i = 0
    conf_a = P.conformal_prediction_set(np.delete(yi, i), np.delete(sei, i), 0.05)
    yi_pert = yi.copy()
    yi_pert[i] = 99.0
    conf_b = P.conformal_prediction_set(np.delete(yi_pert, i), np.delete(sei, i), 0.05)
    assert conf_a['conformal_lo'] == pytest.approx(conf_b['conformal_lo'], abs=1e-12)
    assert conf_a['conformal_hi'] == pytest.approx(conf_b['conformal_hi'], abs=1e-12)
    assert not (conf_b['conformal_lo'] <= 99.0 <= conf_b['conformal_hi'])
