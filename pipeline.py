"""Conformal Meta-Analysis: Held-Out Prediction Calibration for Evidence Synthesis.

Standard prediction intervals assume:
  - Normal random-effects distribution
  - Known variance components
  - Asymptotic approximation (t_{k-2})

Conformal-style prediction sets provide:
  - Held-out empirical calibration without assuming normal random effects
  - Fewer parametric assumptions than standard prediction intervals
  - Diagnostics for multimodal, skewed, heavy-tailed effect distributions
  - Held-out checks for reviews with at least k=5 studies

Method: Split conformal prediction adapted for meta-analysis.
For each review, compute leave-one-out nonconformity scores,
then construct the prediction set from the (1-alpha) quantile.

Applied to eligible Pairwise70 Cochrane reviews. Compare coverage of:
  1. Standard PI (t_{k-2}, assumes normality)
  2. Conformal-style PI (calibrated held-out coverage)
  3. HKSJ PI (Knapp-Hartung adjusted)

Key question: How often does the standard PI FAIL to cover
the effect that a new study would find?
"""

import csv
import json
import math
import os
import time
import numpy as np
import pyreadr
from pathlib import Path
from scipy import stats as sp_stats

ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = ROOT / 'data' / 'output'


def resolve_pairwise_dir() -> Path:
    """Resolve the Pairwise70 data directory without machine-specific paths."""
    candidates = []
    env_dir = os.environ.get('PAIRWISE70_DATA')
    if env_dir:
        candidates.append(Path(env_dir))
    candidates.extend([
        ROOT.parent / 'Pairwise70' / 'data',
        ROOT.parent / 'pairwise70' / 'data',
        ROOT.parent / 'Pairwise70' / 'inst' / 'extdata',
    ])

    for candidate in candidates:
        if candidate.exists() and any(candidate.glob('*.rda')):
            return candidate

    searched = ', '.join(str(path) for path in candidates)
    raise FileNotFoundError(
        'Pairwise70 .rda files were not found. Set PAIRWISE70_DATA or '
        f'place Pairwise70 next to this repo. Searched: {searched}'
    )


# ═══════════════════════════════════════════════════
# CONFORMAL PREDICTION FOR META-ANALYSIS
# ═══════════════════════════════════════════════════

def conformal_prediction_set(yi, sei, alpha=0.05):
    """Compute conformal prediction set for the next study's effect.

    Uses full conformal (leave-one-out) approach:
    1. For each study i, fit DL on the remaining k-1 studies
    2. Compute nonconformity score: |yi - theta_{-i}| / sqrt(sei^2 + tau2_{-i})
    3. The (1-alpha) quantile of scores defines the prediction set radius

    Returns: prediction set [lo, hi] for held-out calibration checks.
    """
    k = len(yi)
    if k < 4:
        return None

    scores = np.zeros(k)

    for i in range(k):
        # Leave-one-out
        yi_loo = np.delete(yi, i)
        sei_loo = np.delete(sei, i)

        # DL on k-1 studies
        wi = 1.0 / sei_loo**2
        sw = np.sum(wi)
        theta_fe = np.sum(wi * yi_loo) / sw
        Q = float(np.sum(wi * (yi_loo - theta_fe)**2))
        C = float(sw - np.sum(wi**2) / sw)
        tau2 = max(0, (Q - (k - 2)) / C) if C > 0 else 0

        ws = 1.0 / (sei_loo**2 + tau2)
        sws = np.sum(ws)
        theta_loo = float(np.sum(ws * yi_loo) / sws)

        # Nonconformity score: standardized residual
        sigma_pred = math.sqrt(sei[i]**2 + tau2)
        scores[i] = abs(yi[i] - theta_loo) / sigma_pred

    # Quantile for prediction set
    # For conformal, use ceil((1-alpha)*(k+1))/k quantile
    q_level = math.ceil((1 - alpha) * (k + 1)) / k
    q_level = min(q_level, 1.0)
    threshold = np.quantile(scores, q_level)

    # Full-data DL for center
    wi = 1.0 / sei**2
    sw = np.sum(wi)
    theta_fe = np.sum(wi * yi) / sw
    Q = float(np.sum(wi * (yi - theta_fe)**2))
    C = float(sw - np.sum(wi**2) / sw)
    tau2 = max(0, (Q - (k - 1)) / C) if C > 0 else 0
    ws = 1.0 / (sei**2 + tau2)
    sws = np.sum(ws)
    theta = float(np.sum(ws * yi) / sws)
    se_theta = float(1.0 / math.sqrt(sws))

    # Prediction set: theta +/- threshold * sigma_pred_new
    # Use median SE as proxy for new study's SE
    se_new = float(np.median(sei))
    sigma_pred_new = math.sqrt(se_new**2 + tau2)
    conformal_lo = theta - threshold * sigma_pred_new
    conformal_hi = theta + threshold * sigma_pred_new

    return {
        'theta': theta,
        'se_theta': se_theta,
        'tau2': tau2,
        'conformal_lo': conformal_lo,
        'conformal_hi': conformal_hi,
        'conformal_width': conformal_hi - conformal_lo,
        'threshold': float(threshold),
        'scores': scores.tolist(),
    }


def standard_prediction_interval(theta, se_theta, tau2, k, alpha=0.05):
    """Standard PI: theta +/- t_{k-2} * sqrt(tau2 + se^2). Assumes normality."""
    if k < 3:
        return None
    t_crit = sp_stats.t.ppf(1 - alpha / 2, k - 2)
    pi_se = math.sqrt(tau2 + se_theta**2)
    return {
        'standard_lo': theta - t_crit * pi_se,
        'standard_hi': theta + t_crit * pi_se,
        'standard_width': 2 * t_crit * pi_se,
    }


def hksj_prediction_interval(yi, sei, theta, tau2, k, alpha=0.05):
    """HKSJ-adjusted PI with the Cochrane/RevMan variance floor."""
    if k < 3:
        return None
    ws = 1.0 / (sei**2 + tau2)
    sws = np.sum(ws)
    q_hksj_raw = float(np.sum(ws * (yi - theta)**2) / (k - 1))
    q_hksj = max(1.0, q_hksj_raw)
    se_hksj = math.sqrt(q_hksj / sws)

    t_crit = sp_stats.t.ppf(1 - alpha / 2, k - 1)
    pi_se = math.sqrt(tau2 + se_hksj**2)
    return {
        'hksj_lo': theta - t_crit * pi_se,
        'hksj_hi': theta + t_crit * pi_se,
        'hksj_width': 2 * t_crit * pi_se,
        'hksj_scale_raw': q_hksj_raw,
        'hksj_scale': q_hksj,
    }


def heldout_interval_coverage(yi, sei):
    """Empirical coverage from intervals fit without the held-out study."""
    k = len(yi)
    covered = {'conformal': 0, 'standard': 0, 'hksj': 0}
    scorable = {'conformal': 0, 'standard': 0, 'hksj': 0}

    for i in range(k):
        yi_train = np.delete(yi, i)
        sei_train = np.delete(sei, i)
        k_train = len(yi_train)

        conf = conformal_prediction_set(yi_train, sei_train, alpha=0.05)
        if conf is not None:
            scorable['conformal'] += 1
            if conf['conformal_lo'] <= yi[i] <= conf['conformal_hi']:
                covered['conformal'] += 1

            std = standard_prediction_interval(
                conf['theta'], conf['se_theta'], conf['tau2'], k_train, 0.05
            )
            if std is not None:
                scorable['standard'] += 1
                if std['standard_lo'] <= yi[i] <= std['standard_hi']:
                    covered['standard'] += 1

            hksj = hksj_prediction_interval(
                yi_train, sei_train, conf['theta'], conf['tau2'], k_train, 0.05
            )
            if hksj is not None:
                scorable['hksj'] += 1
                if hksj['hksj_lo'] <= yi[i] <= hksj['hksj_hi']:
                    covered['hksj'] += 1

    return {
        name: (covered[name] / scorable[name] if scorable[name] else np.nan)
        for name in covered
    }


# ═══════════════════════════════════════════════════
# DATA LOADING
# ═══════════════════════════════════════════════════

def load_review(rda_path):
    result = pyreadr.read_r(str(rda_path))
    df = list(result.values())[0].copy()
    df.columns = df.columns.str.replace(' ', '.', regex=False)
    review_id = rda_path.stem.split('_')[0]

    import pandas as pd
    groups = []
    for (grp, num), sub in df.groupby(['Analysis.group', 'Analysis.number']):
        has_binary = (sub['Experimental.cases'].notna() & (sub['Experimental.cases'] > 0)).any()
        groups.append({'grp': grp, 'num': num, 'k': len(sub), 'binary': has_binary})
    if not groups: return None
    gdf = pd.DataFrame(groups)
    binary = gdf[gdf['binary']]
    best = binary.loc[binary['k'].idxmax()] if len(binary) > 0 else gdf.loc[gdf['k'].idxmax()]
    primary = df[(df['Analysis.group'] == best['grp']) & (df['Analysis.number'] == best['num'])]

    has_binary = (primary['Experimental.cases'].notna() & (primary['Experimental.cases'] > 0)).any()
    scale = 'ratio' if has_binary else ('ratio' if (primary['Mean'].dropna() > 0).all() else 'difference')

    if scale == 'ratio':
        v = (primary['Mean'].notna() & (primary['Mean'] > 0) & primary['CI.start'].notna() & (primary['CI.start'] > 0) & primary['CI.end'].notna() & (primary['CI.end'] > 0))
        sub = primary[v]
        if len(sub) < 4: return None
        yi = np.log(sub['Mean'].values.astype(float))
        sei = (np.log(sub['CI.end'].values.astype(float)) - np.log(sub['CI.start'].values.astype(float))) / (2 * 1.96)
    else:
        v = primary['Mean'].notna() & primary['CI.start'].notna() & primary['CI.end'].notna()
        sub = primary[v]
        if len(sub) < 4: return None
        yi = sub['Mean'].values.astype(float)
        sei = (sub['CI.end'].values.astype(float) - sub['CI.start'].values.astype(float)) / (2 * 1.96)

    ok = (sei > 0) & np.isfinite(yi) & np.isfinite(sei)
    yi, sei = yi[ok], sei[ok]
    if len(yi) < 4: return None
    return {'review_id': review_id, 'yi': yi, 'sei': sei, 'k': len(yi), 'scale': scale}


# ═══════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print("Conformal Meta-Analysis")
    print("=" * 40)

    t0 = time.time()
    pairwise_dir = resolve_pairwise_dir()
    rda_files = sorted(pairwise_dir.glob('*.rda'))
    if not rda_files:
        raise RuntimeError(f'No .rda files found in {pairwise_dir}')

    results = []
    for rda in rda_files:
        review = load_review(rda)
        if review is None or review['k'] < 5:
            continue

        yi, sei, k = review['yi'], review['sei'], review['k']

        # Conformal prediction set
        conf = conformal_prediction_set(yi, sei, alpha=0.05)
        if conf is None:
            continue

        # Standard PI
        std = standard_prediction_interval(conf['theta'], conf['se_theta'], conf['tau2'], k, 0.05)

        # HKSJ PI
        hksj = hksj_prediction_interval(yi, sei, conf['theta'], conf['tau2'], k, 0.05)

        # Held-out coverage for each method
        coverage = heldout_interval_coverage(yi, sei)
        cov_conformal = coverage['conformal']
        cov_standard = coverage['standard']
        cov_hksj = coverage['hksj']

        # Width ratios
        width_ratio_vs_standard = conf['conformal_width'] / std['standard_width'] if std and std['standard_width'] > 0 else 1
        width_ratio_vs_hksj = conf['conformal_width'] / hksj['hksj_width'] if hksj and hksj['hksj_width'] > 0 else 1

        row = {
            'review_id': review['review_id'],
            'k': k,
            'scale': review['scale'],
            'theta': round(conf['theta'], 4),
            'tau2': round(conf['tau2'], 4),
            'conformal_lo': round(conf['conformal_lo'], 4),
            'conformal_hi': round(conf['conformal_hi'], 4),
            'conformal_width': round(conf['conformal_width'], 4),
            'standard_lo': round(std['standard_lo'], 4) if std else '',
            'standard_hi': round(std['standard_hi'], 4) if std else '',
            'standard_width': round(std['standard_width'], 4) if std else '',
            'hksj_lo': round(hksj['hksj_lo'], 4) if hksj else '',
            'hksj_hi': round(hksj['hksj_hi'], 4) if hksj else '',
            'hksj_width': round(hksj['hksj_width'], 4) if hksj else '',
            'hksj_scale_raw': round(hksj['hksj_scale_raw'], 4) if hksj else '',
            'hksj_scale': round(hksj['hksj_scale'], 4) if hksj else '',
            'cov_conformal': round(cov_conformal, 3),
            'cov_standard': round(cov_standard, 3),
            'cov_hksj': round(cov_hksj, 3),
            'width_ratio_conf_std': round(width_ratio_vs_standard, 3),
            'width_ratio_conf_hksj': round(width_ratio_vs_hksj, 3),
            'conformal_wider': width_ratio_vs_standard > 1,
        }
        results.append(row)

    elapsed = time.time() - t0
    n = len(results)
    if n == 0:
        raise RuntimeError(
            f'No eligible reviews were produced from {pairwise_dir}; '
            'check the Pairwise70 schema and eligibility filters.'
        )
    print(f"  Processed: {n} reviews in {elapsed:.1f}s")

    # HEADLINE STATS
    cov_conf = np.array([r['cov_conformal'] for r in results])
    cov_std = np.array([r['cov_standard'] for r in results])
    cov_hksj = np.array([r['cov_hksj'] for r in results])
    width_ratios = np.array([r['width_ratio_conf_std'] for r in results])

    # Coverage failure rate (below nominal 95%)
    conf_undercov = np.sum(cov_conf < 0.90) / n * 100
    std_undercov = np.sum(cov_std < 0.90) / n * 100
    hksj_undercov = np.sum(cov_hksj < 0.90) / n * 100

    conformal_wider = sum(1 for r in results if r['conformal_wider'])

    print(f"\n{'='*55}")
    print("COVERAGE COMPARISON (nominal 95%)")
    print(f"{'='*55}")
    print(f"  {'Method':20s} {'Mean Cov':>10s} {'Median':>10s} {'<90% (fail)':>12s}")
    print(f"  {'Conformal':20s} {np.mean(cov_conf):>10.3f} {np.median(cov_conf):>10.3f} {conf_undercov:>10.1f}%")
    print(f"  {'Standard PI':20s} {np.mean(cov_std):>10.3f} {np.median(cov_std):>10.3f} {std_undercov:>10.1f}%")
    print(f"  {'HKSJ PI':20s} {np.mean(cov_hksj):>10.3f} {np.median(cov_hksj):>10.3f} {hksj_undercov:>10.1f}%")

    print(f"\n  Conformal wider than standard: {conformal_wider}/{n} ({100*conformal_wider/n:.1f}%)")
    print(f"  Mean width ratio (conformal/standard): {np.mean(width_ratios):.2f}")
    print(f"  Median width ratio: {np.median(width_ratios):.2f}")

    # KEY FINDING: reviews where standard PI fails but conformal succeeds
    discordant = sum(1 for r in results if r['cov_standard'] < 0.85 and r['cov_conformal'] >= 0.90)
    print(f"\n  Standard fails (<85% cov) but conformal holds (>=90%): {discordant}/{n} ({100*discordant/n:.1f}%)")

    # EXPORT
    fields = list(results[0].keys())
    with open(OUTPUT_DIR / 'conformal_results.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)

    summary = {
        'n_reviews': n,
        'coverage': {
            'conformal_mean': round(float(np.mean(cov_conf)), 3),
            'standard_mean': round(float(np.mean(cov_std)), 3),
            'hksj_mean': round(float(np.mean(cov_hksj)), 3),
            'conformal_undercov_pct': round(conf_undercov, 1),
            'standard_undercov_pct': round(std_undercov, 1),
            'hksj_undercov_pct': round(hksj_undercov, 1),
        },
        'width': {
            'conformal_wider_pct': round(100 * conformal_wider / n, 1),
            'mean_ratio': round(float(np.mean(width_ratios)), 2),
        },
        'discordant': discordant,
    }
    with open(OUTPUT_DIR / 'conformal_summary.json', 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2)

    print(f"\n  Saved to {OUTPUT_DIR}/")


if __name__ == '__main__':
    main()
