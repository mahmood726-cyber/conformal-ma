"""Known-truth validation tests for conformal vs standard prediction intervals.

These assert the measured invariants behind the repo's headline claim, using
injected ground truth (so coverage is exact, not a held-out proxy).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from truth_recovery_import import run_dist  # noqa: E402


def test_conformal_beats_standard_for_observed_new_effect():
    # k=20, normal: conformal should clearly out-cover the standard PI for y_new.
    c = run_dist("normal", mu=0.3, tau=0.25, k=20, reps=1500, alpha=0.05, seed=1)
    assert c["cov_conformal_ynew"] > c["cov_standard_ynew"]
    assert c["cov_standard_ynew"] < 0.9          # standard PI under-covers y_new
    assert c["cov_conformal_ynew"] >= 0.93       # conformal near nominal at k>=10


def test_conformal_distribution_free_across_shapes():
    # Conformal's guarantee should hold for heavy-tailed / skew / bimodal at k>=10.
    for dist in ("heavy", "skew", "bimodal"):
        c = run_dist(dist, mu=0.3, tau=0.25, k=20, reps=1500, alpha=0.05, seed=2)
        assert c["cov_conformal_ynew"] >= 0.92, f"{dist}: {c['cov_conformal_ynew']}"


def test_standard_pi_targets_theta_not_ynew():
    # Honest decomposition: the standard PI covers the new study's TRUE effect far
    # better than its OBSERVED effect -> it omits new-study sampling variance.
    c = run_dist("normal", mu=0.3, tau=0.25, k=20, reps=1500, alpha=0.05, seed=3)
    assert c["cov_standard_theta_new"] > c["cov_standard_ynew"] + 0.05


def test_deterministic():
    a = run_dist("normal", mu=0.3, tau=0.25, k=10, reps=400, alpha=0.05, seed=9)
    b = run_dist("normal", mu=0.3, tau=0.25, k=10, reps=400, alpha=0.05, seed=9)
    assert a == b
