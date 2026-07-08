"""
Bayes-risk maps comparing eJAB against two reference decision rules.

Model:
    W ~ chi-square_q(lambda)
    lambda = n * theta^2

Here theta is a Wald-scale effect size. For multivariate tests, interpret
theta^2 as the noncentrality per observation after whitening.

Decision rules:
    Berger p-rule:
        reject H0 if p < 0.005
        equivalently W > qchisq(0.995, q)

    eJAB-rule:
        reject H0 if eJAB_01 <= 1/3
        eJAB_01 = sqrt(n) * exp{-0.5 * (1 - n^(-1/q)) * W}
        equivalently W >= log(9n) / (1 - n^(-1/q))

    BIC Bayes-factor approximation:
        BF01_BIC ~= n^(q/2) * exp(-W/2)
        reject H0 if BF01_BIC <= 1/3
        equivalently W >= q * log(n) + log(9)

Bayes risk:
    r(n, theta) = 0.5 * (alpha + beta)
                = 0.5 * (alpha + 1 - power)

Risk advantage:
    delta_r = r_reference - r_ejab
    positive delta_r means eJAB has lower Bayes risk.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import chi2, ncx2


ALPHA_P = 0.005
K01 = 1 / 3
Q_LABELS = {
    1: "Scalar tests: t-test, lm, glm, Cox, Wilcoxon, Mann-Whitney",
    3: "ANOVA / repeated-measures ANOVA example",
    4: "Chi-square / Kruskal-Wallis example",
}


@dataclass(frozen=True)
class Comparison:
    key: str
    file_slug: str
    title_label: str
    risk_label: str
    note: str


@dataclass(frozen=True)
class TestPanel:
    label: str
    q: int


TEST_PANELS = (
    TestPanel("t-test", 1),
    TestPanel("Linear regression", 1),
    TestPanel("Logistic regression", 1),
    TestPanel("Cox PH", 1),
    TestPanel("Wilcoxon", 1),
    TestPanel("Mann-Whitney", 1),
    TestPanel("ANOVA", 3),
    TestPanel("rANOVA", 3),
    TestPanel("Chi-square", 4),
    TestPanel("Kruskal-Wallis", 4),
)


COMPARISONS = (
    Comparison(
        key="berger_p005",
        file_slug="berger_p005",
        title_label="Berger p < 0.005",
        risk_label=r"r_{\mathrm{Berger}}",
        note=(
            "For fixed theta > 0, the unscaled Berger risk difference tends to "
            "alpha_p/2, so this normalized scale tends to 1."
        ),
    ),
    Comparison(
        key="bic_bf01",
        file_slug="bic_bf01",
        title_label="BIC BF01 <= 1/3",
        risk_label=r"r_{\mathrm{BIC}}",
        note="BIC uses BF01 ~= n^(q/2) exp(-W/2), thresholded at BF01 <= 1/3.",
    ),
)


def pvalue_threshold(q: int, alpha_p: float = ALPHA_P) -> float:
    """Critical chi-square statistic for p < alpha_p."""
    return float(chi2.ppf(1 - alpha_p, q))


def ejab_threshold(n, q: int, k01: float = K01):
    """Critical chi-square statistic for eJAB_01 <= k01."""
    n = np.asarray(n, dtype=float)
    if np.any(n <= 1):
        raise ValueError("n must be > 1 so that 1 - n^(-1/q) is positive.")
    return (np.log(n) - 2 * np.log(k01)) / (1 - n ** (-1 / q))


def bic_bf01_threshold(n, q: int, k01: float = K01):
    """Critical chi-square statistic for the BIC BF01 approximation <= k01."""
    n = np.asarray(n, dtype=float)
    if np.any(n <= 0):
        raise ValueError("n must be positive.")
    return q * np.log(n) - 2 * np.log(k01)


def reference_threshold(comparison_key: str, n, q: int, alpha_p: float, k01: float):
    """Return the reference rule's critical chi-square statistic."""
    if comparison_key == "berger_p005":
        return pvalue_threshold(q, alpha_p)
    if comparison_key == "bic_bf01":
        return bic_bf01_threshold(n, q, k01)
    raise ValueError(f"Unknown comparison: {comparison_key}")


def rule_risk(c, q: int, n, theta):
    """
    Bayes risk for a chi-square critical value c.

    c can be scalar or array-valued. n and theta are broadcast through the
    noncentrality lambda = n * theta^2.
    """
    n = np.asarray(n, dtype=float)
    theta = np.asarray(theta, dtype=float)
    lam = n * theta ** 2

    alpha = chi2.sf(c, q)
    power = ncx2.sf(c, q, lam)
    risk = 0.5 * (alpha + 1 - power)
    return risk, alpha, power


def flat_broadcast(value, shape):
    """Flatten value after broadcasting scalars to the requested grid shape."""
    return np.broadcast_to(value, shape).ravel()


def sorted_unique_q(q_values) -> tuple[int, ...]:
    """Return unique q values in ascending numeric order."""
    return tuple(sorted({int(q) for q in q_values}))


def risk_comparison_grid(
    q: int,
    comparison_key: str,
    n_min: float = 20,
    n_max: float = 1e10,
    n_points: int = 900,
    theta_min: float = 0.005,
    theta_max: float = 0.30,
    theta_points: int = 520,
    alpha_p: float = ALPHA_P,
    k01: float = K01,
) -> pd.DataFrame:
    """Create a grid of risks and eJAB risk advantages."""
    n_vals = np.unique(np.round(np.geomspace(n_min, n_max, n_points)).astype(int))
    theta_vals = np.linspace(theta_min, theta_max, theta_points)
    n_grid, theta_grid = np.meshgrid(n_vals, theta_vals)

    c_reference = reference_threshold(comparison_key, n_grid, q, alpha_p, k01)
    c_ejab = ejab_threshold(n_grid, q, k01)

    risk_reference, alpha_reference, power_reference = rule_risk(
        c_reference, q, n_grid, theta_grid
    )
    risk_ejab, alpha_ejab, power_ejab = rule_risk(c_ejab, q, n_grid, theta_grid)

    delta_r = risk_reference - risk_ejab
    normalized_delta_r = delta_r / (alpha_p / 2)

    shape = n_grid.shape
    return pd.DataFrame(
        {
            "comparison": comparison_key,
            "q": q,
            "n": n_grid.ravel(),
            "theta": theta_grid.ravel(),
            "critical_reference": flat_broadcast(c_reference, shape),
            "critical_ejab": flat_broadcast(c_ejab, shape),
            "alpha_reference": flat_broadcast(alpha_reference, shape),
            "alpha_ejab": flat_broadcast(alpha_ejab, shape),
            "power_reference": flat_broadcast(power_reference, shape),
            "power_ejab": flat_broadcast(power_ejab, shape),
            "risk_reference": flat_broadcast(risk_reference, shape),
            "risk_ejab": flat_broadcast(risk_ejab, shape),
            "delta_r": flat_broadcast(delta_r, shape),
            "normalized_delta_r": flat_broadcast(normalized_delta_r, shape),
        }
    )


def pivot_grid(df: pd.DataFrame, z_col: str):
    """Return sorted n/theta axes and a matrix for plotting."""
    n_vals = np.sort(df["n"].unique())
    theta_vals = np.sort(df["theta"].unique())
    z = df.pivot(index="theta", columns="n", values=z_col).loc[theta_vals, n_vals].values
    return n_vals, theta_vals, z


def z_label(comparison: Comparison, normalized: bool) -> str:
    """Colorbar label for a comparison."""
    if normalized:
        return rf"$({comparison.risk_label}-r_{{\mathrm{{eJAB}}}})/(0.005/2)$"
    return rf"${comparison.risk_label}-r_{{\mathrm{{eJAB}}}}$"


def plot_grid(
    df: pd.DataFrame,
    comparison: Comparison,
    title: str,
    normalized: bool = True,
    outfile: Path | None = None,
    clip: tuple[float, float] = (-5, 5),
):
    """Plot one no-contour heatmap for a single q value."""
    z_col = "normalized_delta_r" if normalized else "delta_r"
    n_vals, theta_vals, z = pivot_grid(df, z_col)
    z_plot = np.clip(z, *clip) if normalized else z
    n_grid, theta_grid = np.meshgrid(n_vals, theta_vals)

    fig, ax = plt.subplots(figsize=(8.5, 6))
    norm = TwoSlopeNorm(vmin=clip[0], vcenter=0, vmax=clip[1]) if normalized else None
    mesh = ax.pcolormesh(
        n_grid,
        theta_grid,
        z_plot,
        shading="auto",
        cmap="RdBu_r",
        norm=norm,
    )
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(z_label(comparison, normalized))

    ax.set_xscale("log")
    ax.set_xlabel("Sample size n")
    ax.set_ylabel(r"Effect size $\theta$; noncentrality $\lambda=n\theta^2$")
    ax.set_title(title)

    caption = "Positive values mean eJAB_01 <= 1/3 has lower Bayes risk."
    if normalized:
        caption += " Color clipped to [-5, 5]."
    fig.text(0.5, -0.01, caption, ha="center", fontsize=9)
    fig.tight_layout()

    if outfile:
        fig.savefig(outfile, dpi=300, bbox_inches="tight")
    return fig, ax


def plot_q_panels(
    grids: list[pd.DataFrame],
    comparison: Comparison,
    outfile: Path,
    normalized: bool = True,
    clip: tuple[float, float] = (-5, 5),
):
    """Plot q-ordered panels for one comparison, without contour lines."""
    z_col = "normalized_delta_r" if normalized else "delta_r"
    ordered_grids = sorted(grids, key=lambda grid: int(grid["q"].iloc[0]))

    fig, axes = plt.subplots(
        1,
        len(ordered_grids),
        figsize=(5.3 * len(ordered_grids), 5.3),
        sharey=True,
        constrained_layout=True,
    )
    if len(ordered_grids) == 1:
        axes = [axes]

    norm = TwoSlopeNorm(vmin=clip[0], vcenter=0, vmax=clip[1]) if normalized else None
    mesh = None
    for ax, grid in zip(axes, ordered_grids):
        q = int(grid["q"].iloc[0])
        n_vals, theta_vals, z = pivot_grid(grid, z_col)
        z_plot = np.clip(z, *clip) if normalized else z
        n_grid, theta_grid = np.meshgrid(n_vals, theta_vals)

        mesh = ax.pcolormesh(
            n_grid,
            theta_grid,
            z_plot,
            shading="auto",
            cmap="RdBu_r",
            norm=norm,
        )
        ax.set_xscale("log")
        ax.set_title(f"q = {q}")
        ax.set_xlabel("Sample size n")

    axes[0].set_ylabel(r"Effect size $\theta$")
    fig.suptitle(f"eJAB Bayes-risk advantage vs {comparison.title_label}", fontsize=15)
    fig.supxlabel(
        "Positive values mean eJAB_01 <= 1/3 has lower Bayes risk. "
        + ("Color clipped to [-5, 5]. " if normalized else "")
        + comparison.note,
        fontsize=10,
    )

    if mesh is not None:
        cbar = fig.colorbar(mesh, ax=axes, shrink=0.88, pad=0.015)
        cbar.set_label(z_label(comparison, normalized))

    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    return fig, axes


def plot_test_panels(
    grids: list[pd.DataFrame],
    comparison: Comparison,
    outfile: Path,
    clip: tuple[float, float] = (-5, 5),
):
    """Plot the ten named test panels shown in the manuscript-style figure."""
    grids_by_q = {int(grid["q"].iloc[0]): grid for grid in grids}
    missing_q = sorted({panel.q for panel in TEST_PANELS} - set(grids_by_q))
    if missing_q:
        missing = ", ".join(str(q) for q in missing_q)
        raise ValueError(f"Cannot plot test panels; missing q grid(s): {missing}")

    fig, axes = plt.subplots(
        2,
        5,
        figsize=(21, 10.5),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    contour_color = "#3b0055"
    mesh = None

    for ax, panel in zip(axes.ravel(), TEST_PANELS):
        grid = grids_by_q[panel.q]
        n_vals, theta_vals, z = pivot_grid(grid, "normalized_delta_r")
        z_plot = np.clip(z, *clip)
        n_grid, theta_grid = np.meshgrid(n_vals, theta_vals)

        mesh = ax.pcolormesh(
            n_grid,
            theta_grid,
            z_plot,
            shading="auto",
            cmap="viridis",
            vmin=clip[0],
            vmax=clip[1],
        )
        if np.nanmin(z) <= 0 <= np.nanmax(z):
            ax.contour(
                n_grid,
                theta_grid,
                z,
                levels=[0],
                colors=contour_color,
                linewidths=1.3,
            )
        if np.nanmin(z) <= 1 <= np.nanmax(z):
            ax.contour(
                n_grid,
                theta_grid,
                z,
                levels=[1],
                colors=contour_color,
                linewidths=1.1,
                linestyles="--",
            )

        ax.set_xscale("log")
        ax.set_title(f"{panel.label} (q={panel.q})", fontsize=14)

    for ax in axes[:, 0]:
        ax.set_ylabel("Raw simulation effect size", fontsize=12)
    for ax in axes[1, :]:
        ax.set_xlabel("Nominal sample size n", fontsize=12)

    title_reference = (
        "p < 0.005"
        if comparison.key == "berger_p005"
        else "BIC BF01 <= 1/3"
    )
    fig.suptitle(
        f"Asymptotic Bayes-risk criticality maps: {title_reference} vs eJAB_01 <= 1/3",
        fontsize=20,
    )
    if mesh is not None:
        cbar = fig.colorbar(mesh, ax=axes, shrink=0.72, pad=0.02)
        reference_risk = "risk_p" if comparison.key == "berger_p005" else "risk_BIC"
        cbar.set_label(
            "Normalized risk advantage\n"
            f"({reference_risk} - risk_eJAB) / 0.0025",
            fontsize=12,
        )

    fig.text(
        0.5,
        -0.015,
        "Positive values mean eJAB has lower Bayes risk under the asymptotic "
        "chi-square approximation. Solid contour: equal risk. Dashed contour: "
        "normalized advantage = 1. Color clipped to [-5, 5].",
        ha="center",
        fontsize=12,
    )
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    return fig, axes


def crossing_table(
    comparison_key: str,
    q_values=(1, 3, 4),
    theta_values=(0.01, 0.02, 0.05, 0.10, 0.20),
    n_min: float = 20,
    n_max: float = 1e10,
    n_points: int = 1200,
    alpha_p: float = ALPHA_P,
    k01: float = K01,
) -> pd.DataFrame:
    """
    Approximate smallest grid n after which eJAB risk remains lower
    for all larger n on the grid.
    """
    rows = []
    n_vals = np.unique(np.round(np.geomspace(n_min, n_max, n_points)).astype(int))

    for q in sorted_unique_q(q_values):
        n_float = n_vals.astype(float)
        for theta in theta_values:
            theta_grid = np.full_like(n_float, theta, dtype=float)

            c_reference = reference_threshold(comparison_key, n_float, q, alpha_p, k01)
            c_ejab = ejab_threshold(n_float, q, k01)

            risk_reference, _, _ = rule_risk(c_reference, q, n_float, theta_grid)
            risk_ejab, _, _ = rule_risk(c_ejab, q, n_float, theta_grid)
            delta_r = risk_reference - risk_ejab

            crossing_n = None
            for i, n in enumerate(n_vals):
                if np.all(delta_r[i:] >= 0):
                    crossing_n = int(n)
                    break

            rows.append(
                {
                    "comparison": comparison_key,
                    "q": q,
                    "theta": theta,
                    "approx_min_n_for_ejab_uniformly_lower_risk": crossing_n,
                    "note": "Grid-based; increase n_points/n_max for greater precision.",
                }
            )

    return pd.DataFrame(rows)


def asymptotic_sanity_check(
    comparison_key: str,
    q_values=(1, 3, 4),
    theta_values=(0.01, 0.05, 0.10),
    n_values=(1e6, 1e8, 1e10, 1e12),
    alpha_p: float = ALPHA_P,
    k01: float = K01,
) -> pd.DataFrame:
    """Evaluate risk differences at large n values."""
    rows = []

    for q in sorted_unique_q(q_values):
        for n in n_values:
            for theta in theta_values:
                c_reference = reference_threshold(comparison_key, n, q, alpha_p, k01)
                c_ejab = ejab_threshold(n, q, k01)

                risk_reference, alpha_reference, power_reference = rule_risk(
                    c_reference, q, n, theta
                )
                risk_ejab, alpha_ejab, power_ejab = rule_risk(c_ejab, q, n, theta)

                delta_r = risk_reference - risk_ejab
                rows.append(
                    {
                        "comparison": comparison_key,
                        "q": q,
                        "n": n,
                        "theta": theta,
                        "alpha_reference": float(alpha_reference),
                        "alpha_ejab": float(alpha_ejab),
                        "power_reference": float(power_reference),
                        "power_ejab": float(power_ejab),
                        "risk_reference": float(risk_reference),
                        "risk_ejab": float(risk_ejab),
                        "delta_r": float(delta_r),
                        "normalized_delta_r": float(delta_r / (alpha_p / 2)),
                    }
                )

    return pd.DataFrame(rows)


def parse_q_values(value: str) -> tuple[int, ...]:
    """Parse comma-separated q values."""
    return sorted_unique_q(part.strip() for part in value.split(",") if part.strip())


def build_outputs(args) -> list[Path]:
    """Generate grids, q-ordered plots, crossing tables, and checks."""
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    q_values = parse_q_values(args.q_values)
    written: list[Path] = []

    for comparison in COMPARISONS:
        grids: list[pd.DataFrame] = []

        for q in q_values:
            label = Q_LABELS.get(q, f"q = {q}")
            grid = risk_comparison_grid(
                q=q,
                comparison_key=comparison.key,
                n_min=args.n_min,
                n_max=args.n_max,
                n_points=args.n_points,
                theta_min=args.theta_min,
                theta_max=args.theta_max,
                theta_points=args.theta_points,
                alpha_p=args.alpha_p,
                k01=args.k01,
            )
            grids.append(grid)

            csv_path = outdir / f"ejab_vs_{comparison.file_slug}_risk_grid_q{q}.csv"
            grid.to_csv(csv_path, index=False)
            written.append(csv_path)

            single_plot_path = (
                outdir / f"ejab_vs_{comparison.file_slug}_normalized_risk_advantage_q{q}.png"
            )
            fig, _ = plot_grid(
                grid,
                comparison,
                title=(
                    f"eJAB Bayes-risk advantage vs {comparison.title_label}, "
                    f"q = {q}\n{label}"
                ),
                normalized=True,
                outfile=single_plot_path,
            )
            plt.close(fig)
            written.append(single_plot_path)

        panel_path = outdir / f"ejab_vs_{comparison.file_slug}_normalized_risk_advantage_by_q.png"
        fig, _ = plot_q_panels(grids, comparison, panel_path, normalized=True)
        plt.close(fig)
        written.append(panel_path)

        test_panel_path = (
            outdir / f"ejab_vs_{comparison.file_slug}_10_tests_heatmap_large_n_asymptotic.png"
        )
        fig, _ = plot_test_panels(grids, comparison, test_panel_path)
        plt.close(fig)
        written.append(test_panel_path)

        crossing_path = outdir / f"ejab_vs_{comparison.file_slug}_crossing_table.csv"
        crossing_table(
            comparison.key,
            q_values=q_values,
            n_min=args.n_min,
            n_max=args.n_max,
            n_points=max(args.n_points, 1200),
            alpha_p=args.alpha_p,
            k01=args.k01,
        ).to_csv(crossing_path, index=False)
        written.append(crossing_path)

        sanity_path = outdir / f"ejab_vs_{comparison.file_slug}_asymptotic_sanity_check.csv"
        asymptotic_sanity_check(
            comparison.key,
            q_values=q_values,
            alpha_p=args.alpha_p,
            k01=args.k01,
        ).to_csv(sanity_path, index=False)
        written.append(sanity_path)

    return written


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate eJAB Bayes-risk comparison plots ordered by q."
    )
    parser.add_argument("--outdir", default="bayes_risk_outputs")
    parser.add_argument("--q-values", default="1,3,4")
    parser.add_argument("--n-min", type=float, default=20)
    parser.add_argument("--n-max", type=float, default=1e7)
    parser.add_argument("--n-points", type=int, default=900)
    parser.add_argument("--theta-min", type=float, default=0.005)
    parser.add_argument("--theta-max", type=float, default=0.90)
    parser.add_argument("--theta-points", type=int, default=520)
    parser.add_argument("--alpha-p", type=float, default=ALPHA_P)
    parser.add_argument("--k01", type=float, default=K01)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    written = build_outputs(args)
    print("Wrote:")
    for path in written:
        print(path)


if __name__ == "__main__":
    main()
