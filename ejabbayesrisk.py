#!/usr/bin/env python3
"""
Simulated Bayes-risk criticality maps for:
    p < 0.005
versus
    eJAB_01 <= 1/3
and:
    BIC BF01 <= 1/3
versus
    eJAB_01 <= 1/3

across the ten tests in the eJAB manuscript simulation section.

The script simulates raw datasets, computes test-specific p-values, computes
eJAB_01 and the BIC BF01 approximation from those p-values, estimates alpha
under theta=0 and power under theta>0, then estimates:

    risk(n, theta) = (alpha + 1 - power) / 2

and plots:

    normalized_delta_r = (risk_p - risk_eJAB) / (0.005 / 2)
    normalized_delta_r_BIC = (risk_BIC - risk_eJAB) / (0.005 / 2)

Positive normalized_delta_r means eJAB has lower estimated Bayes risk than
the named reference rule.

Recommended local run:
    python ejabbayesrisk.py --workers 8 --null-reps 20000 --alt-reps 3000

Fast smoke test:
    python ejabbayesrisk.py --workers 4 --null-reps 1000 --alt-reps 200 --fast

Notes:
- The effect-size axis is raw test-specific simulation effect size, not a
  universal Cohen's d across all tests.
- For Cox, this uses a fast univariate Cox score/log-rank-style test and uses
  number of observed events as n_eff for eJAB.
- For publication-quality simulations, increase null reps substantially because
  alpha=0.005 is a rare event.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
import math
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import integrate, stats
from scipy.special import expit
from scipy.stats import chi2


ALPHA_P = 0.005
K01 = 1 / 3
BASE_SEED = 20260625
T3_SCALE = 1 / math.sqrt(3)
HERMITE_X, HERMITE_W = np.polynomial.hermite.hermgauss(80)


TESTS = [
    {
        "key": "ttest",
        "label": "Two-sample t-test",
        "short": "t-test",
        "q": 1,
        "effect_scale": "Equal-group mean difference; approximately Cohen's d.",
    },
    {
        "key": "lm",
        "label": "Simple linear regression",
        "short": "Linear regression",
        "q": 1,
        "effect_scale": "Slope with x ~ N(0,1), error SD 1.",
    },
    {
        "key": "logistic",
        "label": "Simple binary logistic regression",
        "short": "Logistic regression",
        "q": 1,
        "effect_scale": "Log-odds slope with x ~ N(0,1), intercept 0.",
    },
    {
        "key": "cox",
        "label": "Cox proportional-hazards regression",
        "short": "Cox PH",
        "q": 1,
        "effect_scale": "Log-hazard slope with x ~ N(0,1); exponential baseline/censoring.",
    },
    {
        "key": "wilcoxon",
        "label": "Wilcoxon signed-rank test",
        "short": "Wilcoxon",
        "q": 1,
        "effect_scale": "Paired location shift with heavy-tailed t_3 errors.",
    },
    {
        "key": "mannwhitney",
        "label": "Mann-Whitney U-test",
        "short": "Mann-Whitney",
        "q": 1,
        "effect_scale": "Two-sample location shift with heavy-tailed t_3 errors.",
    },
    {
        "key": "anova",
        "label": "One-way ANOVA",
        "short": "ANOVA",
        "q": 3,
        "effect_scale": "Four groups; theta scales centered group means with SD 1.",
    },
    {
        "key": "ranova",
        "label": "One-way repeated-measures ANOVA",
        "short": "rANOVA",
        "q": 3,
        "effect_scale": "Four conditions; theta scales centered condition means; ICC fixed.",
    },
    {
        "key": "chisq",
        "label": "Chi-squared independence test",
        "short": "Chi-square",
        "q": 4,
        "effect_scale": "3x3 multinomial table; theta is diagonal log-linear association.",
    },
    {
        "key": "kruskal",
        "label": "Kruskal-Wallis H-test",
        "short": "Kruskal-Wallis",
        "q": 4,
        "effect_scale": "Five groups; theta scales centered group locations with t_3 errors.",
    },
]


def make_grids(
    fast: bool = False,
    large_n: bool = False,
    n_points: int = 70,
    theta_points: int = 60,
    n_max: float = 1e7,
):
    if fast:
        n_grid = np.array([30, 60, 120, 250, 500, 1000, 2000])
        theta_grid = np.array([0.0, 0.05, 0.10, 0.20, 0.35, 0.50, 0.80])
    elif large_n:
        n_grid = np.unique(np.round(np.geomspace(30, n_max, n_points)).astype(int))

        low_theta = np.geomspace(0.001, 0.02, min(14, theta_points))
        high_count = max(theta_points - len(low_theta), 1)
        high_theta = np.linspace(0.025, 0.9, high_count)
        theta_grid = np.unique(np.r_[0.0, low_theta, high_theta])
    else:
        n_grid = np.array([30, 50, 80, 120, 180, 270, 400, 600, 900,
                           1300, 2000, 3000, 4500, 7000, 10000])
        theta_grid = np.array([0.0, 0.02, 0.05, 0.08, 0.10, 0.15,
                               0.20, 0.30, 0.40, 0.60, 0.80])

    return n_grid, theta_grid


def ejab_w_threshold(n_eff: float, q: int, k01: float = K01) -> float:
    return (math.log(n_eff) - 2 * math.log(k01)) / (1 - n_eff ** (-1 / q))


def ejab01_from_p(p: float, n_eff: float, q: int) -> float:
    """
    eJAB_01 = sqrt(n) exp{-0.5 * (1 - n^{-1/q}) * Q_chisq_q(1-p)}
    """
    if not np.isfinite(p) or not np.isfinite(n_eff) or n_eff <= 1:
        return np.nan
    p_clip = min(max(float(p), np.finfo(float).tiny), 1.0)
    W = chi2.isf(p_clip, q)
    return math.sqrt(n_eff) * math.exp(-0.5 * (1 - n_eff ** (-1 / q)) * W)


def ejab_reject_from_p(p: float, n_eff: float, q: int, k01: float = K01) -> bool | float:
    """
    Reject H0 by eJAB when eJAB_01 <= k01.
    Equivalent critical value:
        W >= (log(n) - 2 log(k01)) / (1 - n^{-1/q})
    """
    if not np.isfinite(p) or not np.isfinite(n_eff) or n_eff <= 1:
        return np.nan

    p_clip = min(max(float(p), np.finfo(float).tiny), 1.0)
    W = chi2.isf(p_clip, q)
    threshold = ejab_w_threshold(n_eff, q, k01)
    return bool(W >= threshold)


def bic_bf01_reject_from_p(p: float, n_eff: float, q: int, k01: float = K01) -> bool | float:
    """
    Reject H0 when the BIC BF01 approximation is <= k01.

    BF01_BIC ~= n_eff^(q/2) exp(-W/2), so the rejection threshold is:
        W >= q log(n_eff) - 2 log(k01)
    """
    if not np.isfinite(p) or not np.isfinite(n_eff) or n_eff <= 1:
        return np.nan

    p_clip = min(max(float(p), np.finfo(float).tiny), 1.0)
    W = chi2.isf(p_clip, q)
    threshold = q * math.log(n_eff) - 2 * math.log(k01)
    return bool(W >= threshold)


def logistic_wald_pvalue(x: np.ndarray, y: np.ndarray,
                          max_iter: int = 40,
                          tol: float = 1e-8) -> float:
    """
    Fast two-parameter logistic regression via IRLS.
    Returns two-sided Wald p-value for slope.
    """
    n = len(y)
    X = np.column_stack((np.ones(n), x))
    beta = np.zeros(2)

    for _ in range(max_iter):
        eta = X @ beta
        mu = expit(eta)
        W = np.clip(mu * (1 - mu), 1e-9, None)
        H = (X.T * W) @ X
        score = X.T @ (y - mu)

        try:
            step = np.linalg.solve(H, score)
        except np.linalg.LinAlgError:
            return np.nan

        beta_new = beta + step

        if np.max(np.abs(step)) < tol:
            beta = beta_new
            break

        beta = beta_new

        # Rough separation guard.
        if np.any(np.abs(beta) > 35):
            return np.nan

    eta = X @ beta
    mu = expit(eta)
    W = np.clip(mu * (1 - mu), 1e-9, None)
    H = (X.T * W) @ X

    try:
        cov = np.linalg.inv(H)
    except np.linalg.LinAlgError:
        return np.nan

    if cov[1, 1] <= 0 or not np.isfinite(cov[1, 1]):
        return np.nan

    z = beta[1] / math.sqrt(cov[1, 1])
    return float(2 * stats.norm.sf(abs(z)))


def cox_score_pvalue(time: np.ndarray, event: np.ndarray, x: np.ndarray) -> tuple[float, int]:
    """
    Fast univariate Cox score-test p-value at beta=0.

    Risk set for subject i after sorting by time is i:m.
    This ignores exact tie corrections, which are negligible for continuous
    exponential event/censoring times.
    """
    order = np.argsort(time)
    e = np.asarray(event, dtype=bool)[order]
    xo = np.asarray(x, dtype=float)[order]
    m = len(xo)

    # Cumulative sums over risk sets.
    risk_n = np.arange(m, 0, -1)
    risk_x = np.cumsum(xo[::-1])[::-1]
    risk_x2 = np.cumsum((xo ** 2)[::-1])[::-1]

    mean_x = risk_x / risk_n
    var_x = risk_x2 / risk_n - mean_x ** 2

    U = np.sum((xo - mean_x)[e])
    I = np.sum(var_x[e])
    n_events = int(e.sum())

    if I <= 1e-12 or not np.isfinite(I) or n_events <= 1:
        return np.nan, n_events

    stat = U ** 2 / I
    return float(chi2.sf(stat, 1)), n_events


def rm_anova_pvalue(Y: np.ndarray) -> float:
    """
    One-way repeated-measures ANOVA p-value from sums of squares.
    Y has shape n_subjects x n_conditions.
    """
    m, k = Y.shape
    grand = Y.mean()
    cond_means = Y.mean(axis=0)
    subj_means = Y.mean(axis=1)

    ss_cond = m * np.sum((cond_means - grand) ** 2)
    ss_subj = k * np.sum((subj_means - grand) ** 2)
    ss_total = np.sum((Y - grand) ** 2)
    ss_error = ss_total - ss_cond - ss_subj

    df_cond = k - 1
    df_error = (m - 1) * (k - 1)

    if ss_error <= 0 or df_error <= 0:
        return np.nan

    F = (ss_cond / df_cond) / (ss_error / df_error)
    return float(stats.f.sf(F, df_cond, df_error))


def simulate_one(test_key: str, n: int, theta: float, rng: np.random.Generator,
                 ranova_icc: float = 0.5) -> tuple[float, int]:
    """
    Simulate one dataset and return (p_value, n_eff_for_eJAB).
    """
    if test_key == "ttest":
        n1 = n // 2
        n2 = n - n1
        x = rng.normal(0, 1, n1)
        y = rng.normal(theta, 1, n2)
        return float(stats.ttest_ind(x, y, equal_var=True).pvalue), n

    if test_key == "lm":
        x = rng.normal(0, 1, n)
        y = theta * x + rng.normal(0, 1, n)
        return float(stats.linregress(x, y).pvalue), n

    if test_key == "logistic":
        x = rng.normal(0, 1, n)
        p = expit(theta * x)
        y = rng.binomial(1, p)
        if y.sum() == 0 or y.sum() == n:
            return np.nan, n
        return logistic_wald_pvalue(x, y), n

    if test_key == "anova":
        k = 4
        sizes = np.full(k, n // k)
        sizes[: n % k] += 1

        contrast = np.arange(k) - (k - 1) / 2
        contrast = contrast / np.sqrt(np.mean(contrast ** 2))

        groups = [
            rng.normal(theta * contrast[j], 1, int(sizes[j]))
            for j in range(k)
        ]

        return float(stats.f_oneway(*groups).pvalue), n

    if test_key == "ranova":
        k = 4

        # The manuscript uses independent observations equal to
        # participants * (conditions - 1). So invert nominal n to subjects.
        m = max(3, int(math.ceil(n / (k - 1))))
        n_eff = m * (k - 1)

        contrast = np.arange(k) - (k - 1) / 2
        contrast = contrast / np.sqrt(np.mean(contrast ** 2))

        subj = rng.normal(0, math.sqrt(ranova_icc), (m, 1))
        eps = rng.normal(0, math.sqrt(1 - ranova_icc), (m, k))
        Y = subj + theta * contrast[None, :] + eps

        return rm_anova_pvalue(Y), n_eff

    if test_key == "chisq":
        # 3x3 multinomial design. Under theta=0, all cells equal.
        # Under theta>0, diagonal cells receive extra log mass.
        r = c = 3
        log_probs = np.zeros((r, c))
        np.fill_diagonal(log_probs, theta)
        probs = np.exp(log_probs)
        probs = probs / probs.sum()

        counts = rng.multinomial(n, probs.ravel()).reshape(r, c)

        try:
            p = stats.chi2_contingency(counts, correction=False).pvalue
        except Exception:
            p = np.nan

        return float(p), n

    if test_key == "cox":
        x = rng.normal(0, 1, n)

        # Exponential PH model: hazard = baseline * exp(theta * x)
        baseline = 0.1
        event_time = rng.exponential(1 / (baseline * np.exp(theta * x)))

        # Gives about 20% censoring under theta=0.
        censor_rate = 0.025
        censor_time = rng.exponential(1 / censor_rate, n)

        time = np.minimum(event_time, censor_time)
        event = event_time <= censor_time

        return cox_score_pvalue(time, event, x)

    if test_key == "wilcoxon":
        # Heavy-tailed paired differences, variance approximately 1.
        d = theta + rng.standard_t(3, n) / math.sqrt(3)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                p = stats.wilcoxon(d, alternative="two-sided", method="auto").pvalue
        except Exception:
            p = np.nan
        return float(p), n

    if test_key == "mannwhitney":
        n1 = n // 2
        n2 = n - n1
        x = rng.standard_t(3, n1) / math.sqrt(3)
        y = theta + rng.standard_t(3, n2) / math.sqrt(3)
        p = stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic").pvalue
        return float(p), n

    if test_key == "kruskal":
        k = 5
        sizes = np.full(k, n // k)
        sizes[: n % k] += 1

        contrast = np.arange(k) - (k - 1) / 2
        contrast = contrast / np.sqrt(np.mean(contrast ** 2))

        groups = [
            theta * contrast[j] + rng.standard_t(3, int(sizes[j])) / math.sqrt(3)
            for j in range(k)
        ]

        return float(stats.kruskal(*groups).pvalue), n

    raise ValueError(f"Unknown test key: {test_key}")


def normal_expectation(func) -> float:
    """
    Deterministic Gauss-Hermite expectation for X ~ N(0, 1).
    """
    x = math.sqrt(2) * HERMITE_X
    return float(np.sum(HERMITE_W * func(x)) / math.sqrt(math.pi))


def t3_pdf(x: float | np.ndarray) -> float | np.ndarray:
    return stats.t.pdf(np.asarray(x) / T3_SCALE, df=3) / T3_SCALE


def t3_cdf(x: float | np.ndarray) -> float | np.ndarray:
    return stats.t.cdf(np.asarray(x) / T3_SCALE, df=3)


@lru_cache(maxsize=None)
def logistic_ncp_per_n(theta: float) -> float:
    if theta == 0:
        return 0.0
    mean_score = normal_expectation(lambda x: x * (expit(theta * x) - 0.5))
    return 4 * mean_score ** 2


@lru_cache(maxsize=None)
def chisq_ncp_per_n(theta: float) -> float:
    if theta == 0:
        return 0.0

    log_probs = np.zeros((3, 3))
    np.fill_diagonal(log_probs, theta)
    probs = np.exp(log_probs)
    probs = probs / probs.sum()

    expected = probs.sum(axis=1, keepdims=True) @ probs.sum(axis=0, keepdims=True)
    return float(np.sum((probs - expected) ** 2 / expected))


@lru_cache(maxsize=None)
def cox_event_fraction(theta: float) -> float:
    baseline = 0.1
    censor_rate = 0.025
    return normal_expectation(
        lambda x: baseline * np.exp(theta * x) / (baseline * np.exp(theta * x) + censor_rate)
    )


@lru_cache(maxsize=None)
def mannwhitney_ncp_per_n(theta: float) -> float:
    if theta == 0:
        return 0.0

    def integrand(e):
        return t3_cdf(theta + e) * t3_pdf(e)

    area, _ = integrate.quad(
        integrand,
        -np.inf,
        np.inf,
        epsabs=1e-10,
        epsrel=1e-8,
        limit=200,
    )
    return float(3 * (area - 0.5) ** 2)


@lru_cache(maxsize=None)
def wilcoxon_ncp_per_n(theta: float) -> float:
    if theta == 0:
        return 0.0

    def abs_cdf_shifted(x):
        return t3_cdf(x - theta) - t3_cdf(-x - theta)

    def integrand(d):
        return abs_cdf_shifted(d) * t3_pdf(d - theta)

    area, _ = integrate.quad(
        integrand,
        0,
        np.inf,
        epsabs=1e-10,
        epsrel=1e-8,
        limit=200,
    )
    return float(12 * (area - 0.25) ** 2)


@lru_cache(maxsize=None)
def kruskal_ncp_per_n(theta: float) -> float:
    if theta == 0:
        return 0.0

    k = 5
    contrast = np.arange(k) - (k - 1) / 2
    contrast = contrast / np.sqrt(np.mean(contrast ** 2))
    locations = theta * contrast

    def pooled_cdf(x):
        return np.mean([t3_cdf(x - loc) for loc in locations], axis=0)

    rank_means = []
    for loc in locations:
        def integrand(x):
            return pooled_cdf(x) * t3_pdf(x - loc)

        rank_mean, _ = integrate.quad(
            integrand,
            -np.inf,
            np.inf,
            epsabs=1e-10,
            epsrel=1e-8,
            limit=200,
        )
        rank_means.append(rank_mean)

    rank_means = np.asarray(rank_means)
    return float(12 * np.mean((rank_means - 0.5) ** 2))


def asymptotic_n_eff(test_key: str, n: int, theta: float) -> float:
    if test_key == "ranova":
        return max(3, int(math.ceil(n / 3))) * 3
    if test_key == "cox":
        return max(2.0, n * cox_event_fraction(theta))
    return float(n)


def asymptotic_ncp(test_key: str, n: int, theta: float) -> float:
    if theta == 0:
        return 0.0

    if test_key == "ttest":
        return n * theta ** 2 / 4
    if test_key == "lm":
        return n * theta ** 2
    if test_key == "logistic":
        return n * logistic_ncp_per_n(theta)
    if test_key == "anova":
        return n * theta ** 2
    if test_key == "ranova":
        n_eff = asymptotic_n_eff(test_key, n, theta)
        return n_eff * (8 / 3) * theta ** 2
    if test_key == "chisq":
        return n * chisq_ncp_per_n(theta)
    if test_key == "cox":
        # Local Cox score approximation under about 80% observed events.
        return n * 0.8 * theta ** 2
    if test_key == "wilcoxon":
        return n * wilcoxon_ncp_per_n(theta)
    if test_key == "mannwhitney":
        return n * mannwhitney_ncp_per_n(theta)
    if test_key == "kruskal":
        return n * kruskal_ncp_per_n(theta)

    raise ValueError(f"Unknown test key: {test_key}")


def simulate_cell(task: dict) -> dict:
    """
    Worker function for one test x n x theta cell.
    """
    cfg = task["cfg"]
    n = int(task["n"])
    theta = float(task["theta"])
    reps = int(task["reps"])
    seed = int(task["seed"])

    rng = np.random.default_rng(seed)

    reject_p = []
    reject_ejab = []
    reject_bic = []
    p_values = []
    n_eff_values = []
    ejab_values = []

    for _ in range(reps):
        p, n_eff = simulate_one(cfg["key"], n, theta, rng)

        if not np.isfinite(p) or not np.isfinite(n_eff):
            continue

        ejab_reject = ejab_reject_from_p(p, n_eff, cfg["q"])

        if isinstance(ejab_reject, float) and np.isnan(ejab_reject):
            continue

        bic_reject = bic_bf01_reject_from_p(p, n_eff, cfg["q"])

        if isinstance(bic_reject, float) and np.isnan(bic_reject):
            continue

        ejab_value = ejab01_from_p(p, n_eff, cfg["q"])

        reject_p.append(p < ALPHA_P)
        reject_ejab.append(bool(ejab_reject))
        reject_bic.append(bool(bic_reject))
        p_values.append(p)
        n_eff_values.append(n_eff)
        ejab_values.append(ejab_value)

    valid_reps = len(reject_p)

    if valid_reps == 0:
        rejection_p = np.nan
        rejection_ejab = np.nan
        rejection_bic = np.nan
        mean_p = np.nan
        median_p = np.nan
        mean_n_eff = np.nan
        mean_ejab01 = np.nan
    else:
        rejection_p = float(np.mean(reject_p))
        rejection_ejab = float(np.mean(reject_ejab))
        rejection_bic = float(np.mean(reject_bic))
        mean_p = float(np.mean(p_values))
        median_p = float(np.median(p_values))
        mean_n_eff = float(np.mean(n_eff_values))
        mean_ejab01 = float(np.nanmean(ejab_values))

    return {
        "test": cfg["label"],
        "short": cfg["short"],
        "test_key": cfg["key"],
        "q": cfg["q"],
        "effect_scale": cfg["effect_scale"],
        "n_nominal": n,
        "theta": theta,
        "target_reps": reps,
        "valid_reps": valid_reps,
        "rejection_p": rejection_p,
        "rejection_ejab": rejection_ejab,
        "rejection_bic": rejection_bic,
        "mean_p": mean_p,
        "median_p": median_p,
        "mean_n_eff": mean_n_eff,
        "mean_ejab01": mean_ejab01,
        "method": "raw",
    }


def simulate_cell_asymptotic(task: dict) -> dict:
    """
    Analytic large-sample approximation for one test x n x theta cell.

    The p-value threshold is evaluated against a central/noncentral chi-square
    approximation to the test statistic. This avoids materializing raw datasets
    at n values such as 10^5 or 10^7.
    """
    cfg = task["cfg"]
    n = int(task["n"])
    theta = float(task["theta"])
    reps = int(task["reps"])
    q = int(cfg["q"])

    n_eff = asymptotic_n_eff(cfg["key"], n, theta)
    p_threshold = chi2.isf(ALPHA_P, q)
    ejab_threshold = ejab_w_threshold(n_eff, q)
    bic_threshold = q * math.log(n_eff) - 2 * math.log(K01)

    if theta == 0:
        rejection_p = float(chi2.sf(p_threshold, q))
        rejection_ejab = float(chi2.sf(ejab_threshold, q))
        rejection_bic = float(chi2.sf(bic_threshold, q))
    else:
        ncp = asymptotic_ncp(cfg["key"], n, theta)
        rejection_p = float(stats.ncx2.sf(p_threshold, q, ncp))
        rejection_ejab = float(stats.ncx2.sf(ejab_threshold, q, ncp))
        rejection_bic = float(stats.ncx2.sf(bic_threshold, q, ncp))

    return {
        "test": cfg["label"],
        "short": cfg["short"],
        "test_key": cfg["key"],
        "q": cfg["q"],
        "effect_scale": cfg["effect_scale"],
        "n_nominal": n,
        "theta": theta,
        "target_reps": reps,
        "valid_reps": reps,
        "rejection_p": rejection_p,
        "rejection_ejab": rejection_ejab,
        "rejection_bic": rejection_bic,
        "mean_p": np.nan,
        "median_p": np.nan,
        "mean_n_eff": n_eff,
        "mean_ejab01": np.nan,
        "method": "asymptotic",
    }


def simulate_cell_large_n_mc(task: dict) -> dict:
    """
    Large-n Monte Carlo cell using each test's large-sample statistic model.

    This simulates rejection counts from the central/noncentral chi-square
    rejection probabilities. It is intended for fine grids and very large n,
    where raw row-level dataset simulation is computationally impractical.
    """
    cfg = task["cfg"]
    n = int(task["n"])
    theta = float(task["theta"])
    reps = int(task["reps"])
    q = int(cfg["q"])
    rng = np.random.default_rng(int(task["seed"]))

    n_eff = asymptotic_n_eff(cfg["key"], n, theta)
    p_threshold = chi2.isf(ALPHA_P, q)
    ejab_threshold = ejab_w_threshold(n_eff, q)
    bic_threshold = q * math.log(n_eff) - 2 * math.log(K01)

    if theta == 0:
        prob_p = float(chi2.sf(p_threshold, q))
        prob_ejab = float(chi2.sf(ejab_threshold, q))
        prob_bic = float(chi2.sf(bic_threshold, q))
    else:
        ncp = asymptotic_ncp(cfg["key"], n, theta)
        prob_p = float(stats.ncx2.sf(p_threshold, q, ncp))
        prob_ejab = float(stats.ncx2.sf(ejab_threshold, q, ncp))
        prob_bic = float(stats.ncx2.sf(bic_threshold, q, ncp))

    return {
        "test": cfg["label"],
        "short": cfg["short"],
        "test_key": cfg["key"],
        "q": cfg["q"],
        "effect_scale": cfg["effect_scale"],
        "n_nominal": n,
        "theta": theta,
        "target_reps": reps,
        "valid_reps": reps,
        "rejection_p": float(rng.binomial(reps, prob_p) / reps),
        "rejection_ejab": float(rng.binomial(reps, prob_ejab) / reps),
        "rejection_bic": float(rng.binomial(reps, prob_bic) / reps),
        "mean_p": np.nan,
        "median_p": np.nan,
        "mean_n_eff": n_eff,
        "mean_ejab01": np.nan,
        "method": "large_n_mc",
    }


def build_tasks(n_grid: np.ndarray,
                theta_grid: np.ndarray,
                null_reps: int,
                alt_reps: int,
                selected_tests: set[str] | None = None) -> list[dict]:
    tasks = []

    for test_idx, cfg in enumerate(TESTS):
        if selected_tests is not None and cfg["key"] not in selected_tests:
            continue

        for n_idx, n in enumerate(n_grid):
            for theta_idx, theta in enumerate(theta_grid):
                reps = null_reps if theta == 0 else alt_reps
                seed = BASE_SEED + 100000 * test_idx + 1000 * n_idx + 10 * theta_idx

                tasks.append({
                    "cfg": cfg,
                    "n": int(n),
                    "theta": float(theta),
                    "reps": int(reps),
                    "seed": int(seed),
                })

    return tasks


def estimate_risks(raw: pd.DataFrame) -> pd.DataFrame:
    """
    Use theta=0 cells as alpha estimates, theta>0 cells as power estimates.
    """
    null = (
        raw[raw["theta"] == 0]
        [["test_key", "n_nominal", "rejection_p", "rejection_ejab", "rejection_bic"]]
        .rename(columns={
            "rejection_p": "alpha_p_hat",
            "rejection_ejab": "alpha_ejab_hat",
            "rejection_bic": "alpha_bic_hat",
        })
    )

    alt = raw[raw["theta"] > 0].merge(null, on=["test_key", "n_nominal"], how="left")

    alt["power_p_hat"] = alt["rejection_p"]
    alt["power_ejab_hat"] = alt["rejection_ejab"]
    alt["power_bic_hat"] = alt["rejection_bic"]

    alt["risk_p_hat"] = 0.5 * (alt["alpha_p_hat"] + 1 - alt["power_p_hat"])
    alt["risk_ejab_hat"] = 0.5 * (alt["alpha_ejab_hat"] + 1 - alt["power_ejab_hat"])
    alt["risk_bic_hat"] = 0.5 * (alt["alpha_bic_hat"] + 1 - alt["power_bic_hat"])

    alt["delta_r_hat"] = alt["risk_p_hat"] - alt["risk_ejab_hat"]
    alt["normalized_delta_r_hat"] = alt["delta_r_hat"] / (ALPHA_P / 2)
    alt["delta_r_bic_hat"] = alt["risk_bic_hat"] - alt["risk_ejab_hat"]
    alt["normalized_delta_r_bic_hat"] = alt["delta_r_bic_hat"] / (ALPHA_P / 2)

    return alt


def plot_heatmap(
    risk: pd.DataFrame,
    outdir: Path,
    value_col: str,
    reference_label: str,
    reference_short: str,
    outfile_prefix: str,
    suffix: str = "",
) -> Path:
    """
    Ten-panel heatmap of normalized risk advantage, sorted by q.
    """
    present_keys = list(risk["test_key"].drop_duplicates())
    present_tests = [cfg for cfg in TESTS if cfg["key"] in present_keys]

    n_panels = len(present_tests)
    ncols = 5 if n_panels > 5 else n_panels
    nrows = int(math.ceil(n_panels / ncols))

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(3.6 * ncols, 3.6 * nrows),
        sharex=True,
        sharey=True,
        squeeze=False,
    )

    axes_flat = axes.ravel()

    last_mesh = None
    vmin, vmax = -5, 5

    for ax, cfg in zip(axes_flat, present_tests):
        d = risk[risk["test_key"] == cfg["key"]].copy()

        pivot = (
            d.pivot(index="theta", columns="n_nominal", values=value_col)
            .sort_index()
            .sort_index(axis=1)
        )

        x = pivot.columns.values.astype(float)
        y = pivot.index.values.astype(float)
        Z = pivot.values
        Z_plot = np.clip(Z, vmin, vmax)

        # Cell edges for pcolormesh.
        if len(x) > 1:
            x_log = np.log10(x)
            x_edges = 10 ** np.r_[
                x_log[0] - (x_log[1] - x_log[0]) / 2,
                (x_log[:-1] + x_log[1:]) / 2,
                x_log[-1] + (x_log[-1] - x_log[-2]) / 2,
            ]
        else:
            x_edges = np.array([x[0] / 1.2, x[0] * 1.2])

        if len(y) > 1:
            y_edges = np.r_[
                y[0] - (y[1] - y[0]) / 2,
                (y[:-1] + y[1:]) / 2,
                y[-1] + (y[-1] - y[-2]) / 2,
            ]
        else:
            y_edges = np.array([y[0] - 0.01, y[0] + 0.01])

        last_mesh = ax.pcolormesh(
            x_edges,
            y_edges,
            Z_plot,
            shading="auto",
            vmin=vmin,
            vmax=vmax,
        )

        ax.set_xscale("log")
        ax.set_title(f"{cfg['short']} (q={cfg['q']})", fontsize=10)
        ax.tick_params(axis="both", labelsize=8)

    for ax in axes_flat[n_panels:]:
        ax.axis("off")

    for ax in axes[-1, :]:
        ax.set_xlabel("Nominal sample size n", fontsize=9)

    for row in axes:
        row[0].set_ylabel("Raw simulation effect size", fontsize=9)

    method = None
    if "method" in risk.columns and risk["method"].nunique(dropna=True) == 1:
        method = str(risk["method"].dropna().iloc[0])
    title_prefix = "Asymptotic Bayes-risk" if method == "asymptotic" else "Empirical Bayes-risk"

    fig.suptitle(
        f"{title_prefix} criticality maps: {reference_label} vs eJAB_01 <= 1/3",
        fontsize=14,
        y=1.02,
    )

    fig.subplots_adjust(
        left=0.06,
        right=0.86,
        top=0.88,
        bottom=0.08,
        wspace=0.18,
        hspace=0.28,
    )

    cax = fig.add_axes([0.895, 0.18, 0.014, 0.62])
    cbar = fig.colorbar(last_mesh, cax=cax)
    cbar.set_label(
        f"Normalized risk advantage\n({reference_short} - risk_eJAB) / 0.0025",
        fontsize=9,
        labelpad=8,
    )
    cbar.ax.tick_params(labelsize=8)

    png = outdir / f"{outfile_prefix}{suffix}.png"

    fig.savefig(png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return png


def parse_int_list(value: str) -> list[int]:
    return sorted({int(part.strip()) for part in value.split(",") if part.strip()})


def make_theoretical_grids(
    n_min: float,
    n_max: float,
    n_points: int,
    theta_min: float,
    theta_max: float,
    theta_points: int,
) -> tuple[np.ndarray, np.ndarray]:
    n_grid = np.unique(np.round(np.geomspace(n_min, n_max, n_points)).astype(int))

    low_count = min(14, theta_points)
    low_theta = np.geomspace(theta_min, min(0.02, theta_max), low_count)
    high_count = max(theta_points - low_count, 1)
    high_start = min(0.025, theta_max)
    high_theta = np.linspace(high_start, theta_max, high_count)
    theta_grid = np.unique(np.r_[low_theta, high_theta])

    return n_grid, theta_grid


def theoretical_threshold(reference: str, n, q: int):
    if reference == "berger_p005":
        return chi2.isf(ALPHA_P, q)
    if reference == "bic_bf01":
        n = np.asarray(n, dtype=float)
        return q * np.log(n) - 2 * math.log(K01)
    raise ValueError(f"Unknown theoretical reference: {reference}")


def theoretical_normalized_advantage(
    q: int,
    n_grid: np.ndarray,
    theta_grid: np.ndarray,
    reference: str,
) -> np.ndarray:
    n_mesh, theta_mesh = np.meshgrid(n_grid.astype(float), theta_grid.astype(float))
    lam = n_mesh * theta_mesh ** 2

    c_ref = theoretical_threshold(reference, n_mesh, q)
    c_ejab = (np.log(n_mesh) - 2 * math.log(K01)) / (1 - n_mesh ** (-1 / q))

    alpha_ref = chi2.sf(c_ref, q)
    alpha_ejab = chi2.sf(c_ejab, q)
    power_ref = stats.ncx2.sf(c_ref, q, lam)
    power_ejab = stats.ncx2.sf(c_ejab, q, lam)

    risk_ref = 0.5 * (alpha_ref + 1 - power_ref)
    risk_ejab = 0.5 * (alpha_ejab + 1 - power_ejab)
    return (risk_ref - risk_ejab) / (ALPHA_P / 2)


def plot_theoretical_q_panels(
    q_values: list[int],
    n_grid: np.ndarray,
    theta_grid: np.ndarray,
    reference: str,
    outdir: Path,
) -> Path:
    reference_label = "p < 0.005" if reference == "berger_p005" else "BIC BF01 <= 1/3"
    reference_short = "risk_p" if reference == "berger_p005" else "risk_BIC"
    slug = "berger_p005" if reference == "berger_p005" else "bic_bf01"
    q_slug = "_".join(str(q) for q in q_values)

    fig, axes = plt.subplots(
        1,
        len(q_values),
        figsize=(4.2 * len(q_values), 4.6),
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    axes_flat = axes.ravel()
    last_mesh = None
    vmin, vmax = -5, 5

    x = n_grid.astype(float)
    x_log = np.log10(x)
    x_edges = 10 ** np.r_[
        x_log[0] - (x_log[1] - x_log[0]) / 2,
        (x_log[:-1] + x_log[1:]) / 2,
        x_log[-1] + (x_log[-1] - x_log[-2]) / 2,
    ]

    y = theta_grid.astype(float)
    y_edges = np.r_[
        y[0] - (y[1] - y[0]) / 2,
        (y[:-1] + y[1:]) / 2,
        y[-1] + (y[-1] - y[-2]) / 2,
    ]
    y_edges[0] = max(0, y_edges[0])

    for ax, q in zip(axes_flat, q_values):
        z = theoretical_normalized_advantage(q, n_grid, theta_grid, reference)
        last_mesh = ax.pcolormesh(
            x_edges,
            y_edges,
            np.clip(z, vmin, vmax),
            shading="auto",
            vmin=vmin,
            vmax=vmax,
        )
        ax.set_xscale("log")
        ax.set_title(f"q = {q}", fontsize=11)
        ax.tick_params(axis="both", labelsize=8)
        ax.set_xlabel("Sample size n", fontsize=9)

    axes_flat[0].set_ylabel(r"Wald-scale effect size $\theta$", fontsize=9)
    fig.suptitle(
        f"Theoretical Bayes-risk criticality maps: {reference_label} vs eJAB_01 <= 1/3",
        fontsize=14,
        y=1.02,
    )
    fig.subplots_adjust(left=0.055, right=0.89, top=0.82, bottom=0.17, wspace=0.12)

    cax = fig.add_axes([0.91, 0.21, 0.012, 0.56])
    cbar = fig.colorbar(last_mesh, cax=cax)
    cbar.set_label(
        f"Normalized risk advantage\n({reference_short} - risk_eJAB) / 0.0025",
        fontsize=9,
        labelpad=8,
    )
    cbar.ax.tick_params(labelsize=8)

    outpath = outdir / f"ejab_vs_{slug}_theoretical_q_{q_slug}.png"
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return outpath


def write_theoretical_outputs(args) -> list[Path]:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    q_values = parse_int_list(args.theoretical_q_values)
    n_grid, theta_grid = make_theoretical_grids(
        n_min=args.theoretical_n_min,
        n_max=args.theoretical_n_max,
        n_points=args.theoretical_n_points,
        theta_min=args.theoretical_theta_min,
        theta_max=args.theoretical_theta_max,
        theta_points=args.theoretical_theta_points,
    )

    config_path = outdir / (
        "ejab_theoretical_q_config_"
        + "_".join(str(q) for q in q_values)
        + ".csv"
    )
    pd.DataFrame(
        {
            "q": q_values,
            "alpha_p": ALPHA_P,
            "k01": K01,
            "n_min": float(n_grid.min()),
            "n_max": float(n_grid.max()),
            "n_points": len(n_grid),
            "theta_min": float(theta_grid.min()),
            "theta_max": float(theta_grid.max()),
            "theta_points": len(theta_grid),
        }
    ).to_csv(config_path, index=False)

    outputs = [config_path]
    outputs.append(plot_theoretical_q_panels(q_values, n_grid, theta_grid, "berger_p005", outdir))
    outputs.append(plot_theoretical_q_panels(q_values, n_grid, theta_grid, "bic_bf01", outdir))
    return outputs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=str, default="bayes_risk_outputs")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--null-reps", type=int, default=20000)
    parser.add_argument("--alt-reps", type=int, default=3000)
    parser.add_argument("--fast", action="store_true")
    parser.add_argument("--n-points", type=int, default=70)
    parser.add_argument("--theta-points", type=int, default=60)
    parser.add_argument("--n-max", type=float, default=1e7)
    parser.add_argument("--theoretical", action="store_true")
    parser.add_argument("--theoretical-q-values", type=str, default="1,2,4,16,256")
    parser.add_argument("--theoretical-n-min", type=float, default=30)
    parser.add_argument("--theoretical-n-max", type=float, default=1e7)
    parser.add_argument("--theoretical-n-points", type=int, default=900)
    parser.add_argument("--theoretical-theta-min", type=float, default=0.001)
    parser.add_argument("--theoretical-theta-max", type=float, default=0.90)
    parser.add_argument("--theoretical-theta-points", type=int, default=520)
    parser.add_argument(
        "--large-n",
        action="store_true",
        help="Use a fine large-n grid and Monte Carlo rejection counts from test-specific large-sample distributions.",
    )
    parser.add_argument(
        "--asymptotic",
        action="store_true",
        help="Use analytic large-sample rejection probabilities instead of Monte Carlo or raw-data simulation.",
    )
    parser.add_argument(
        "--tests",
        type=str,
        default="all",
        help="Comma-separated test keys or 'all'. Keys: " + ",".join([t["key"] for t in TESTS]),
    )
    args = parser.parse_args()

    if args.theoretical:
        outputs = write_theoretical_outputs(args)
        print("Done.")
        for output in outputs:
            print(f"Wrote: {output}")
        return

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    n_grid, theta_grid = make_grids(
        args.fast,
        args.large_n,
        n_points=args.n_points,
        theta_points=args.theta_points,
        n_max=args.n_max,
    )
    use_asymptotic = args.asymptotic
    use_large_n_mc = args.large_n and not args.asymptotic

    selected_tests = None
    if args.tests != "all":
        selected_tests = set([x.strip() for x in args.tests.split(",") if x.strip()])
        valid = set(t["key"] for t in TESTS)
        bad = selected_tests - valid
        if bad:
            raise ValueError(f"Unknown test keys: {bad}. Valid keys: {valid}")

    tasks = build_tasks(
        n_grid=n_grid,
        theta_grid=theta_grid,
        null_reps=args.null_reps,
        alt_reps=args.alt_reps,
        selected_tests=selected_tests,
    )

    print(f"Running {len(tasks)} cells with workers={args.workers}")
    if use_asymptotic:
        method_label = "asymptotic"
    elif use_large_n_mc:
        method_label = "large_n_mc"
    else:
        method_label = "raw"
    print(f"Simulation method: {method_label}")
    print(f"Null reps per theta=0 cell: {args.null_reps}")
    print(f"Alt reps per theta>0 cell: {args.alt_reps}")
    if args.large_n:
        print(f"Large-n grid max: {int(n_grid.max())}")

    rows = []
    if use_asymptotic:
        worker = simulate_cell_asymptotic
    elif use_large_n_mc:
        worker = simulate_cell_large_n_mc
    else:
        worker = simulate_cell
    progress_interval = max(10, len(tasks) // 100)

    if args.workers == 1:
        for i, task in enumerate(tasks, start=1):
            rows.append(worker(task))
            if i % progress_interval == 0 or i == len(tasks):
                print(f"Completed {i}/{len(tasks)} cells")
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futures = [ex.submit(worker, task) for task in tasks]
            for i, fut in enumerate(as_completed(futures), start=1):
                rows.append(fut.result())
                if i % progress_interval == 0 or i == len(tasks):
                    print(f"Completed {i}/{len(tasks)} cells")

    raw = pd.DataFrame(rows)
    risk = estimate_risks(raw)

    suffix_parts = []
    if args.fast:
        suffix_parts.append("fast")
    if args.large_n:
        suffix_parts.append("large_n")
    if use_large_n_mc:
        suffix_parts.append("mc")
    if use_asymptotic:
        suffix_parts.append("asymptotic")
    suffix = "_" + "_".join(suffix_parts) if suffix_parts else ""

    raw_csv = outdir / f"ejab_simulated_10_tests_rejections{suffix}.csv"
    risk_csv = outdir / f"ejab_simulated_10_tests_risk_grid{suffix}.csv"
    config_csv = outdir / f"ejab_simulated_10_tests_config{suffix}.csv"

    raw.to_csv(raw_csv, index=False)
    risk.to_csv(risk_csv, index=False)
    pd.DataFrame(TESTS).to_csv(config_csv, index=False)

    p_plot = plot_heatmap(
        risk,
        outdir,
        value_col="normalized_delta_r_hat",
        reference_label="p < 0.005",
        reference_short="risk_p",
        outfile_prefix="ejab_vs_berger_p005_10_tests_heatmap_empirical",
        suffix=suffix,
    )
    bic_plot = plot_heatmap(
        risk,
        outdir,
        value_col="normalized_delta_r_bic_hat",
        reference_label="BIC BF01 <= 1/3",
        reference_short="risk_BIC",
        outfile_prefix="ejab_vs_bic_bf01_10_tests_heatmap_empirical",
        suffix=suffix,
    )

    print("Done.")
    print(f"Wrote: {raw_csv}")
    print(f"Wrote: {risk_csv}")
    print(f"Wrote: {config_csv}")
    print(f"Wrote: {p_plot}")
    print(f"Wrote: {bic_plot}")


if __name__ == "__main__":
    main()
