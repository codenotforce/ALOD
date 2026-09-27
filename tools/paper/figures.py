"""The user's tuned September plotting parameters, without host-specific paths.

Data adapters own file formats; this module only consumes canonical run records.
Separate rc_contexts prevent E1/E2 serif settings from leaking into the E3 plot.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
from matplotlib.ticker import NullFormatter
from .selection import METHODS, KAPPAS, completed_cycles, displayed, promotions, target_prefix

STYLE_PROVENANCE = {
    "E1_E2": {
        "script": "generate_kernel_lift_figures.py",
        "sha256": "e54b2b48941673ac87935206f263a4ef276f69f012104d37d1958160f83224c1",
    },
    "E3": {
        "script": "plot_results.py",
        "sha256": "18b105625b977d3805dc0a783474f1137173b3ffe68a0a7a28926185a2eb6bc9",
        "renamed_output": "kappa_error_dof.pdf",
    },
    "sampling": "ALOD: initial and completed H-steps in E1/E2/E3. E3 cutoff uses the first E <= 0.01 among all accepted states; baselines keep every point in that prefix.",
}
E12_STYLE = {
    "font.family": "serif",
    "mathtext.fontset": "dejavuserif",
    "font.size": 10,
    "axes.labelsize": 10,
    "legend.fontsize": 8,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "savefig.dpi": 220,
}
E3_STYLE = {
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.labelsize": 10,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.dpi": 200,
}
COLORS = {
    "ALOD": "#2166ac",
    "AFEM": "#25815f",
    "UFEM": "#777777",
    "LOD": "#8c56a0",
    "eta": "#cd7725",
}
MARKERS = {"ALOD": "o", "AFEM": "^", "UFEM": "v", "LOD": "s"}
K_COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#333333")
K_MARKERS = ("o", "s", "^", "D", "v")


def error_figure(exp, runs):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.65), layout="constrained")
    markers = dict(MARKERS, ALOD="D" if exp == "e2" else "o")
    for method in METHODS:
        run = runs.get(f"{exp}-{method}")
        if run:
            rows = (
                completed_cycles(run)
                if method == "ALOD"
                else displayed(exp, method, run, runs)["rows"]
            )
            axes[0].loglog(
                [r["N_on"] for r in rows],
                [r["E"] for r in rows],
                color=COLORS[method],
                lw=1.35,
                marker=markers[method],
                ms=4,
                markevery=4 if method == "AFEM" else 1,
                label=method,
            )
    run = runs.get(f"{exp}-ALOD")
    if run:
        rows = completed_cycles(run)
        for check in promotions(run):
            row = next((r for r in rows if r["state"] >= check["state_id"]), rows[-1])
            axes[0].annotate(
                rf'$\ell:{check["ell"]}\to{check["next_ell"]}$',
                xy=(row["N_on"], row["E"]),
                xytext=(row["N_on"] * 0.5, row["E"] * 0.4),
                arrowprops=dict(arrowstyle="->", color=".35", lw=0.8),
                fontsize=9,
            )
        axes[1].loglog(
            [r["N_on"] for r in rows],
            [r["E"] for r in rows],
            color=COLORS["ALOD"],
            lw=1.45,
            marker=markers["ALOD"],
            ms=3.5,
            label="ALOD relative error",
        )
        eta = [r for r in rows if r.get("eta_relative") is not None]
        if eta:
            axes[1].loglog(
                [r["N_on"] for r in eta],
                [r["eta_relative"] for r in eta],
                color=COLORS["eta"],
                lw=1.35,
                ls="--",
                label=r"Coarse indicator: $\eta_H/\|u\|_\kappa$",
            )
    for ax in axes:
        if exp == "e2":
            ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel(r"Online unknowns $N_{\rm on}$")
        ax.set_ylabel("Relative energy error")
        ax.grid(alpha=0.2, which="major")
        ax.legend(frameon=False, fontsize=14)
    return fig


def draw_mesh(ax, geometry, *, e2=False, region=False):
    points = np.asarray(geometry["nodes"])
    tri = np.asarray(geometry["elements"], dtype=int)
    if (
        tri.ndim != 2
        or tri.shape[1] != 3
        or tri.size == 0
        or tri.min() < 0
        or tri.max() >= len(points)
    ):
        raise ValueError("invalid mesh connectivity")
    # Draw interior edges once. Drawing each triangle would darken shared edges
    # and change the user's line widths, particularly on fine reference meshes.
    edges = np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    unique, counts = np.unique(np.sort(edges, axis=1), axis=0, return_counts=True)
    width = 0.3 if len(tri) < 15000 else 0.04 if e2 else 0.1
    ax.add_collection(
        LineCollection(points[unique], color="#42556b", linewidths=width, rasterized=True)
    )
    ax.add_collection(
        LineCollection(
            points[unique[counts == 1]], color="#32485e", linewidths=0.55, rasterized=True
        )
    )
    if region:
        angle = np.linspace(0, 1.5 * np.pi, 400)
        arc = 0.6 * np.c_[np.cos(angle), np.sin(angle)]
        ax.add_patch(
            Polygon(np.vstack([[0, 0], arc, [0, 0]]), color="#efad46", alpha=0.16, zorder=2)
        )
        ax.plot(arc[:, 0], arc[:, 1], "--", color="#b7471c", lw=1.4, zorder=4)
        ax.annotate(
            "",
            xy=(-0.424264, 0.424264),
            xytext=(0, 0),
            arrowprops=dict(arrowstyle="<->", color="#b7471c"),
        )
        ax.text(
            0.20,
            -0.18,
            r"$R=0.6$",
            color="#a83411",
            ha="left",
            va="top",
            bbox=dict(fc="white", ec="none", alpha=0.88),
            zorder=5,
        )
    limits = (-1.03, 1.03) if e2 else (-0.01, 1.01)
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_aspect("equal")
    ax.set_axis_off()


def mesh_figure(exp, runs, issues):
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.65), layout="constrained")
    for ax, (method, part) in zip(
        axes, [("ALOD", "coarse"), ("ALOD", "reference"), ("AFEM", "coarse")]
    ):
        geometry = runs.get(f"{exp}-{method}", {}).get("mesh", {}).get(part)
        if geometry:
            draw_mesh(ax, geometry, e2=exp == "e2", region=exp == "e2" and method == "ALOD")
        else:
            issues.append(f"{exp}-{method}: missing {part} mesh")
            ax.text(0.5, 0.5, "Missing mesh", ha="center")
            ax.axis("off")
    return fig


def kappa_figure(runs):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), layout="constrained")
    for ax, method in zip(axes, ("ALOD", "AFEM", "LOD")):
        for k, color, marker in zip(KAPPAS, K_COLORS, K_MARKERS):
            run = runs.get(f"k{k}-{method}")
            if not run:
                continue
            rows = target_prefix(run)
            if method == "ALOD":
                rows = completed_cycles(dict(rows=rows))
            x = np.array([r["N_on"] for r in rows], float) / (k * k)
            y = np.array([r["E"] for r in rows])
            ax.loglog(
                x,
                y,
                color=color,
                marker=marker,
                markevery=max(1, len(rows) // 10),
                markersize=5.5,
                markeredgewidth=1.2,
                markerfacecolor="white",
                linewidth=1.8,
                label=rf"$k={k}$",
            )
            ax.plot(
                x[-1],
                y[-1],
                marker=marker,
                color=color,
                markersize=6.5,
                markeredgewidth=1.2,
                markerfacecolor="white",
                linestyle="none",
            )
        ax.axhline(0.01, color="#777777", linestyle=":", linewidth=0.9, label="0.01 target")
        ax.set_xlabel(r"$N_{\mathrm{on}}/k^2$", fontsize=12)
        ax.set_ylabel(method + ": relative energy error", fontsize=12)
        ax.grid(True, which="major", alpha=0.22)
        ax.grid(True, which="minor", alpha=0.07)
        ax.set_ylim(0.005, 2.2)
        ax.legend(loc="upper right", fontsize=12, framealpha=0.9)
    return fig


def save(fig, output, name, *, untitled_pdf=False):
    # Constrained layout is updated during savefig. Preserve each original
    # script's save order as well as its style to avoid subpixel layout drift.
    extensions = ("png", "pdf") if untitled_pdf else ("pdf", "png")
    for ext in extensions:
        options = {"metadata": {"Title": None}} if ext == "pdf" and untitled_pdf else {}
        fig.savefig(output / f"{name}.{ext}", bbox_inches="tight", **options)
    plt.close(fig)


def render_figures(runs, output, issues):
    # Start from Matplotlib defaults for each independent original script.
    with plt.rc_context({**matplotlib.rcParamsDefault, **E12_STYLE}):
        for exp in ("e1", "e2"):
            save(error_figure(exp, runs), output, exp + "_error_dof")
            save(mesh_figure(exp, runs, issues), output, exp + "_final_meshes")
    with plt.rc_context({**matplotlib.rcParamsDefault, **E3_STYLE}):
        save(kappa_figure(runs), output, "kappa_error_dof", untitled_pdf=True)
