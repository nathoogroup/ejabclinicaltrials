#!/usr/bin/env python3
"""D_n vs n under the strict rule: original p < .005 and every replication
attempt p < .005, only p and n required. Adapted from
R2Empirical/Scripts/plot_replications.py (3v4n-5/documents), plot_points() only.

D_n = -sqrt(n) * p * ln(p), computed from ln(p) for both series.

Reads:  ../Data/{EERP,RPP,RPCB,SCORE}_orig.csv, *_repl.csv, counts.csv
        (built by Scripts/build_strict_data.R)
Writes: ../Output/Dn_vs_n_points_maxn_uncensored_logy_strict005_FINALFINAL.png,
            the figure to use: log n and log D_n axes, censored p-values removed
            (reported as a bound, e.g. "p < .001"; see the censored column
            written by build_strict_data.R), largest-n point circled and
            labelled with its D_n, condition R2 stated below it
        ../Output/Dn_vs_n_points_maxn_uncensored_strict005.png: the same with a
            linear D_n axis
        ../Output/Dn_vs_n_points_maxn_uncensored_logy_quantile_strict005.png: the
            figure to use with a quantile-regression median line and
            10th-90th percentile band
        ../Output/Dn_vs_n_points_maxn_uncensored_logy_polyfit_strict005.png: the
            figure to use with a quadratic least-squares line, log10 n centred
        ../Output/Dn_vs_n_points_maxn_strict005_FINAL.png: linear D_n axis with
            the censored p-values kept
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "Data"
OUT = ROOT / "Output"
DATASETS = ["EERP", "RPP", "RPCB", "SCORE"]
COLORS = {"EERP": "#9467bd", "RPP": "#1f77b4", "RPCB": "#2ca02c", "SCORE": "#d62728"}

D_N_LABEL = "$D_n = -\\sqrt{n}\\,p\\,\\ln p$"
# condition R2 as stated in the manuscript (REVISION/REVISED.tex, Theorem 1)
R2_LABEL = ("$\\mathbf{R2}$: $D_n=-\\sqrt{n}\\,p\\ln p\\ \\overset{P}{\\to}\\ 0$ "
            "as $n\\to\\infty$ under $\\mathcal{H}_1$")


def load(name, kind, drop_censored=False):
    d = pd.read_csv(DATA / f"{name}_{kind}.csv")
    if drop_censored:
        d = d[~d["censored"]]
    return d[np.isfinite(d["logp"]) & (d["logp"] < 0) & (d["n"] > 0)]


def dn_log(n, logp):
    n = np.asarray(n, dtype=float)
    lp = np.asarray(logp, dtype=float)
    return -np.sqrt(n) * np.exp(lp) * lp


def log10_dn(n, logp):
    """log10 of D_n, on the log scale throughout: five points (RPP 46, three
    SCORE replications) have D_n below 1e-308, which underflows to 0."""
    n = np.asarray(n, dtype=float)
    lp = np.asarray(logp, dtype=float)
    return 0.5 * np.log10(n) + (lp + np.log(-lp)) / np.log(10)


def add_quantile_fit(ax, n, l10):
    """Quantile regression of log10 D_n on log10 n: median line and 10th-90th
    percentile band. Quantile regression, not least squares, so the extreme
    points (RPP 46 at 1e-1166) do not drive the fit. The linear 10th and 90th
    lines cross below n ~ 5, so the band starts at n = 5."""
    import statsmodels.api as sm
    x = np.log10(np.asarray(n, dtype=float))
    y = np.asarray(l10, dtype=float)
    X = sm.add_constant(x)
    fit = {q: sm.QuantReg(y, X).fit(q=q) for q in (0.1, 0.5, 0.9)}
    xs = np.linspace(np.log10(5), x.max(), 200)
    line = {q: f.params[0] + f.params[1] * xs for q, f in fit.items()}
    lo, hi = fit[0.5].conf_int()[1]
    ax.fill_between(10 ** xs, line[0.1], line[0.9], color="0.5", alpha=0.15,
                    lw=0, zorder=1, label="10th–90th percentile")
    ax.plot(10 ** xs, line[0.5], color="black", lw=1.5, zorder=1,
            label=f"median (slope {fit[0.5].params[1]:.1f})")
    for q in (0.1, 0.9):
        ax.plot(10 ** xs, line[q], color="0.4", lw=0.8, ls="--", zorder=1)


def add_poly_fit(ax, n, l10, degree=2):
    """Least-squares polynomial of log10 D_n on log10 n, with log10 n centred
    at its mean. Least squares, so the extreme points (RPP 46 at 1e-1166)
    pull the line below most of the data."""
    x = np.log10(np.asarray(n, dtype=float))
    y = np.asarray(l10, dtype=float)
    x0 = x.mean()
    c = np.polyfit(x - x0, y, degree)
    r2 = 1 - (y - np.polyval(c, x - x0)).var() / y.var()
    xs = np.linspace(x.min(), x.max(), 200)
    ax.plot(10 ** xs, np.polyval(c, xs - x0), color="black", lw=1.5, zorder=1,
            label=f"degree-{degree} fit ($R^2$ = {r2:.2f})")


def lighten(color, f=0.55):
    c = np.array(matplotlib.colors.to_rgb(color))
    return tuple(c + (1.0 - c) * f)


def cap_n_axis(ax, n, dn):
    """Linear n axis only: cap at the 95th percentile of n, rescale D_n to the
    points still in view. Returns an x label saying how many points fall
    beyond the cap."""
    n, dn = np.asarray(n, dtype=float), np.asarray(dn, dtype=float)
    cap = np.percentile(n, 95)
    ax.set_xlim(-0.02 * cap, 1.02 * cap)
    top = 1.05 * dn[n <= cap].max()
    ax.set_ylim(-0.02 * top, top)
    return (f"sample size $n$ (capped at 95th percentile, {cap:,.0f}; "
            f"{(n > cap).sum()} of {len(n)} points beyond)")


def plot_points(xscale="log", filename="Dn_vs_n_points_strict005.png",
                mark_max=False, r2=False, drop_censored=False, logy=False,
                quantiles=False, poly=None):
    counts = pd.read_csv(DATA / "counts.csv").set_index("project")
    fig, ax = plt.subplots(figsize=(10, 6.5))
    # logy plots log10 D_n on a linear axis labelled in powers of ten
    y = log10_dn if logy else dn_log
    all_n, all_dn = [], []
    for name in DATASETS:
        o = load(name, "orig", drop_censored)
        total = counts.loc[name, "results_p_n"]
        ax.scatter(o["n"], y(o["n"], o["logp"]), s=24, color=COLORS[name],
                   alpha=0.7, edgecolor="none", zorder=2,
                   label=f"{name} original ({len(o)}/{total})")
        r = load(name, "repl", drop_censored)
        ax.scatter(r["n"], y(r["n"], r["logp"]), s=26, marker="^",
                   facecolor=lighten(COLORS[name]), edgecolor=COLORS[name],
                   linewidth=0.5, alpha=0.85, zorder=2,
                   label=f"{name} replication ({len(r)})")
        for d in (o, r):
            all_n.extend(d["n"])
            all_dn.extend(y(d["n"], d["logp"]))

    ax.set_xscale(xscale)
    ax.set_xlabel("sample size $n$")
    ax.set_ylabel(D_N_LABEL)
    # on a linear n axis the highest points sit at the left edge, under an upper-left legend
    if quantiles:                      # before the legend, so the fit is listed
        add_quantile_fit(ax, all_n, all_dn)
    if poly:
        add_poly_fit(ax, all_n, all_dn, poly)
    # on a log D_n axis the points crowd the top and the bottom left is empty
    ax.legend(title="study", ncol=2, fontsize=9, framealpha=0.9,
              loc="lower left" if logy else
                  "upper left" if xscale == "log" else "best")
    ax.set_title("$D_n$ vs $n$, original $p\\leq0.005$ & every replication $p\\leq0.005$"
                 + (", censored $p$ removed" if drop_censored else ""))
    if xscale == "linear":
        ax.set_xlabel(cap_n_axis(ax, all_n, all_dn))
    elif logy:
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
            lambda v, _: f"$10^{{{v:.0f}}}$"))
        ax.set_ylabel(D_N_LABEL + " (log scale)")
    else:
        top = ax.get_ylim()[1]
        ax.set_ylim(-0.02 * top, top)
    if mark_max:                       # circle the largest-n point and label it
        k = int(np.argmax(all_n))
        ax.scatter([all_n[k]], [all_dn[k]], s=260, facecolor="none",
                   edgecolor="black", linewidth=1.2, zorder=3)
        man, exp = f"{(10 ** all_dn[k] if logy else all_dn[k]):.2e}".split("e")
        # (obs) marks the observed value, apart from the D_n of the R2 statement
        label = f"$D_{{\\max n}}^{{(\\mathrm{{obs}})}} = {man}\\times10^{{{int(exp)}}}$"
        ax.annotate(f"{label}\n{R2_LABEL}" if r2 else label,
                    (all_n[k], all_dn[k]),
                    # on a log D_n axis the point sits at the top: label below it
                    xytext=(-40, -260 if logy else 70), textcoords="offset points",
                    ha="right",
                    fontsize=12, zorder=3, multialignment="center", linespacing=2.0,
                    arrowprops=dict(arrowstyle="->", color="black", lw=1.0,
                                    shrinkB=9))
    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    path = OUT / filename
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print("saved:", path)


def plot_by_study(xscale="log", shared=True,
                  filename="Dn_vs_n_by_study_strict005.png"):
    """The same scatter, one panel per project (2x2). shared=True puts every
    panel on the same x and y axes; False lets each panel scale to its own
    points."""
    counts = pd.read_csv(DATA / "counts.csv").set_index("project")
    fig, axs = plt.subplots(2, 2, figsize=(14, 9.5), sharex=shared, sharey=shared)
    for ax, name in zip(axs.flat, DATASETS):
        o = load(name, "orig")
        r = load(name, "repl")
        ax.scatter(o["n"], dn_log(o["n"], o["logp"]), s=24, color=COLORS[name],
                   alpha=0.7, edgecolor="none", zorder=2,
                   label=f"original ({len(o)}/{counts.loc[name, 'results_p_n']})")
        ax.scatter(r["n"], dn_log(r["n"], r["logp"]), s=26, marker="^",
                   facecolor=lighten(COLORS[name]), edgecolor=COLORS[name],
                   linewidth=0.5, alpha=0.85, zorder=2,
                   label=f"replication ({len(r)})")
        ax.set_title(name)
        ax.set_xscale(xscale)
        if xscale == "linear":
            ax.set_xlabel(cap_n_axis(ax, pd.concat([o["n"], r["n"]]),
                                     np.concatenate([dn_log(o["n"], o["logp"]),
                                                     dn_log(r["n"], r["logp"])])))
        ax.legend(fontsize=9, framealpha=0.9,
                  loc="upper left" if xscale == "log" else "best")
    if xscale == "log":
        for ax in (axs.flat if not shared else axs[1]):
            ax.set_xlabel("sample size $n$")
    for ax in (axs.flat if not shared else axs[:, 0]):
        ax.set_ylabel(D_N_LABEL)
    if xscale == "log":
        for ax in (axs.flat if not shared else [ax]):
            top = ax.get_ylim()[1]
            ax.set_ylim(-0.02 * top, top)
    fig.suptitle("$D_n$ vs $n$ by study, original $p\\leq0.005$ & every "
                 "replication $p\\leq0.005$")
    fig.tight_layout()
    path = OUT / filename
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print("saved:", path)


if __name__ == "__main__":
    # the figure to use (2026-10-23); the variants above are kept but not drawn
    plot_points("log", "Dn_vs_n_points_maxn_uncensored_logy_strict005_FINALFINAL.png",
                mark_max=True, r2=True, drop_censored=True, logy=True)
    plot_points("log", "Dn_vs_n_points_maxn_uncensored_strict005.png",
                mark_max=True, r2=True, drop_censored=True)
    plot_points("log", "Dn_vs_n_points_maxn_uncensored_logy_quantile_strict005.png",
                mark_max=True, r2=True, drop_censored=True, logy=True, quantiles=True)
    plot_points("log", "Dn_vs_n_points_maxn_uncensored_logy_polyfit_strict005.png",
                mark_max=True, r2=True, drop_censored=True, logy=True, poly=2)
    plot_points("log", "Dn_vs_n_points_maxn_strict005_FINAL.png", mark_max=True,
                r2=True)
