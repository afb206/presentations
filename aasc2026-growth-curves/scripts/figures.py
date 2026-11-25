# Stage 2: redraw the figures for the web from the saved files in sources/ (plot only).
# Text stays as real text in the SVG (so it is searchable and takes the page font).
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyreadr
from common import *

plt.rcParams.update({
    "svg.fonttype": "none", "font.family": "sans-serif",
    "font.sans-serif": ["Public Sans", "Segoe UI", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 11, "axes.labelsize": 11, "xtick.labelsize": 10, "ytick.labelsize": 10,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 1.1,
    "legend.frameon": False,
})
INK = "#222222"
OUT = os.path.join(HERE, "figs")
os.makedirs(OUT, exist_ok=True)

def save(fig, name):
    fig.savefig(os.path.join(OUT, name), format="svg", bbox_inches="tight")
    plt.close(fig)
    print("wrote figs/" + name)

# ---------------------------------------------------------------- Figure 1: simulation, size band
def fig_sim():
    m3 = rd("mm3_fits_all.csv")
    m3 = m3[(m3.half == "test") & (m3.condition == "clean") & m3.variant.isin(["V0", "MJ"])]
    cov, wid = {}, {}
    for v in ["MJ", "V0"]:
        g = m3[m3.variant == v].groupby("design")
        cov[v] = g.cov_y.mean().reindex(DES)
        wid[v] = g.w_y.median().reindex(DES)
    P1, P2 = rd("P1_cov_y.csv"), rd("P2_w_y.csv")
    for a in ["log_nls", "natural_nls"]:
        cov[a] = P1[P1.arm == a].set_index("group").est.reindex(DES)
        wid[a] = P2[P2.arm == a].set_index("group").est.reindex(DES)
    G1, G2 = rd("G1_cov_y.csv"), rd("G2_w_y.csv")
    cov["spline_log"] = G1[G1.arm == "log_spline"].set_index("design").est.reindex(DES)
    wid["spline_log"] = G2[G2.arm == "log_spline"].set_index("design").est.reindex(DES)
    order = ["MJ", "V0", "log_nls", "natural_nls", "spline_log"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 4.1), gridspec_kw={"wspace": 0.28})
    x = list(range(4))
    for ax, d in [(a1, cov), (a2, wid)]:
        for k in order:
            hollow = k in ("spline_log", "natural_nls")
            ax.plot(x, d[k].values, color=COL[k], lw=2.0,
                    ls="--" if k == "spline_log" else ":" if k == "natural_nls" else "-",
                    marker="s" if k == "natural_nls" else "o", ms=6.5, mfc="white" if hollow else COL[k], mew=1.6,
                    label=LAB[k] + (" (other replicates)" if k == "spline_log" else ""))
        ax.set_xticks(x); ax.set_xticklabels(DES); ax.set_xlabel("Sampling design: dense to sparse")
    a1.axhline(0.95, color=INK, lw=1.2, ls=":"); a1.text(0.0, 0.93, "nominal 0.95", fontsize=9, va="top")
    a1.set_ylim(0, 1.05); a1.set_ylabel("Coverage of the true size curve")
    a2.set_yscale("log"); a2.set_ylabel("Median band width, % of K")
    a1.set_title("Coverage", loc="left", weight="bold"); a2.set_title("Band width", loc="left", weight="bold")
    h, l = a1.get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=10)
    save(fig, "fig1_sim.svg")

# ---------------------------------------------------------------- Figure 2: hold-out error by horizon
def fig_horizon():
    b = rd("F4b_H3_error_by_day_log_pairing.csv")
    G = ["Day 31", "Day 34", "Days 36-44"]
    T_ = ["V0", "MJ", "log_nls", "spline_log", "carry", "loglin"]
    fig, ax = plt.subplots(figsize=(9.0, 4.3))
    w = 0.12
    for j, t in enumerate(T_):
        d = b[b.tool == t].set_index("group").reindex(G)
        xs = [i + (j - 2.5) * w for i in range(3)]
        base = t in ("carry", "loglin")
        ax.errorbar(xs, d.est, yerr=[d.est - d.lo, d.hi - d.est], fmt="D" if base else "o", ms=6.5,
                    color=COL[t], mfc="white" if base else COL[t], mew=1.6, elinewidth=1.6, capsize=0, label=LAB[t])
    ax.set_xticks(range(3)); ax.set_xticklabels(["Day 31", "Day 34", "Days 36 to 44"])
    ax.set_xlabel("Left-out readings, after the Day 30 cut (H3)")
    ax.set_ylabel("Median absolute log error"); ax.set_ylim(0, 3.0)
    h, l = ax.get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", bbox_to_anchor=(0.5, -0.17), ncol=3, fontsize=10)
    save(fig, "fig2_horizon.svg")

# ---------------------------------------------------------------- paired GF arms (coverage of the size band by design)
def fig_paired_gf():
    P1 = rd("P1_cov_y.csv")
    des = DES + ["sparse"]
    arms = [("MJ", COL["MJ"], "o", "-"), ("log_nls", COL["log_nls"], "o", "-"), ("natural_nls", COL["natural_nls"], "o", "-"),
            ("log_published", COL["log_nls"], "s", ":"), ("natural_published", COL["natural_nls"], "s", ":"),
            ("log_improved", COL["log_nls"], "^", "--"), ("natural_improved", COL["natural_nls"], "^", "--")]
    names = {"MJ": "GP+L", "log_nls": "GF (log)", "natural_nls": "GF (nat)", "log_published": "GF (log), published settings",
             "natural_published": "GF (nat), published settings", "log_improved": "GF (log), improved (EXTRA)", "natural_improved": "GF (nat), improved (EXTRA)"}
    fig, ax = plt.subplots(figsize=(9.0, 4.4))
    for a, c, mk, ls in arms:
        d = P1[P1.arm == a].set_index("group").reindex(des)
        worked = d.lab.str.extract(r"^(\d+)/")[0].astype(float)
        ax.plot(range(len(des)), d.est, ls, color=c, lw=1.4, zorder=1)
        ax.scatter(range(len(des)), d.est, marker=mk, s=42, color=[c if w > 0 else "white" for w in worked], edgecolor=c, linewidth=1.4, zorder=2, label=names[a])
    ax.axhline(0.95, color=INK, lw=0.8, ls=(0, (2, 2)))
    ax.set_xticks(range(len(des))); ax.set_xticklabels(des)
    ax.set_xlabel("Sampling design (fewer readings to the right)"); ax.set_ylabel("Coverage of the true size curve (95% band)")
    ax.set_ylim(-0.03, 1.05)
    h, l = ax.get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize=9.5)
    save(fig, "fig5_paired_gf.svg")

# ---------------------------------------------------------------- spike: ratio of interval score, spiked over clean
def fig_spike():
    c = rd("mm3_c1_c5.csv").set_index("variant")
    lv = [("0.05", "005"), ("0.12", "012"), ("0.30", "030")]
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    w = 0.34
    for j, (v, col, lab) in enumerate([("V0", COL["V0"], "GP"), ("MJ", COL["MJ"], "GP+L")]):
        ys = [c.loc["MJ", ("c4_ratio0_spike_" + k) if v == "V0" else ("c4_ratio_spike_" + k)] for _, k in lv]
        xs = [i + (j - 0.5) * w for i in range(3)]
        ax.bar(xs, ys, width=w, color=col, label=lab)
        for x_, y_ in zip(xs, ys): ax.text(x_, y_ + 0.2, f"{y_:.1f}", ha="center", fontsize=9.5)
    ax.set_xticks(range(3)); ax.set_xticklabels([l for l, _ in lv])
    ax.set_xlabel("Spike size s (spike = s × the 90th percentile of that dataset’s readings)"); ax.set_ylabel("Interval score, spiked / clean (median)")
    ax.legend(loc="upper left")
    save(fig, "fig6_spike.svg")

# ---------------------------------------------------------------- Arabidopsis description: fit / Day 44 reading, and the GP+L logistic midpoint
def fig_arab_desc():
    d = rd("arab_mj_set_per_fit.csv")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.0, 4.0))
    bins = np.linspace(0.8, 1.6, 33)
    for arm, c, lab in [("V0", COL["V0"], "GP"), ("MJ", COL["MJ"], "GP+L")]:
        x = d[(d.arm == arm)].d44_ratio.dropna().clip(0.8, 1.6)
        a1.hist(x, bins=bins, color=c, alpha=0.75, label=f"{lab} (n = {len(x)})")
    a1.axvline(1, color=INK, lw=0.9, ls=(0, (2, 2)))
    a1.set_xlabel("Fitted size / the Day 44 reading (in-sample; values above 1.6 shown at 1.6)"); a1.set_ylabel("Plants"); a1.legend()
    m = d[d.arm == "MJ"].tm_med.dropna()
    a2.hist(m, bins=np.arange(0, 40, 1.0), color=COL["MJ"])
    a2.set_xlabel("GP+L logistic midpoint, posterior median (day)"); a2.set_ylabel("Plants (n = %d)" % len(m))
    fig.tight_layout()
    save(fig, "fig7_arab_desc.svg")

# ---------------------------------------------------------------- Figure 3: inclusion against width
def fig_calib():
    c = rd("F3_inclusion_width.csv")
    L = ["H1", "H2", "H3"]
    fig, ax = plt.subplots(figsize=(9.0, 4.4))
    for t in ["spline_natural", "natural_nls", "natural_improved", "V0", "MJ", "log_nls", "spline_log"]:
        d = c[c.tool == t].set_index("level").reindex(L)
        faint = t in ("spline_natural", "natural_nls", "natural_improved")
        ax.plot(d.wlog, 100 * d.incl, color=COL[t] if t != "spline_natural" else COL["spline_log"],
                lw=1.8, alpha=0.4 if faint else 1, ls=":" if faint else "-", marker="o", ms=6,
                mfc="white" if faint else COL[t], label=LAB[t] + (" (faint)" if faint else ""))
        ax.annotate("H3", (d.wlog["H3"], 100 * d.incl["H3"]), textcoords="offset points", xytext=(5, 3),
                    fontsize=9, color=COL[t], alpha=0.6 if faint else 1)
    ax.axhline(90, color=INK, lw=1.2, ls=":"); ax.text(0.03, 91, "nominal 90%", fontsize=9, va="bottom")
    ax.set_xlabel("Interval width: ln(upper / lower)"); ax.set_ylabel("Left-out readings inside the 90% interval (%)")
    ax.set_ylim(0, 105)
    h, l = ax.get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize=10)
    save(fig, "fig3_calibration.svg")

# ---------------------------------------------------------------- Figure 4: one plant, four tools
def fig_plant():
    cur = rd("abstract_1003_curves.csv")
    cur = cur[cur.t >= 4]                      # the Day 1 to 3 extrapolation is left out, as on the poster
    r = pyreadr.read_r(os.path.join(SRC, "arabidopsis_rosette_area_valid.rds"))
    d = list(r.values())[0]
    pl = d[d["Plant ID"] == "1135FA1003"]
    tools = ["log_nls", "spline_log", "V0", "MJ"]
    fig, axs = plt.subplots(2, 4, figsize=(10.4, 5.2), sharex=True, gridspec_kw={"hspace": 0.12, "wspace": 0.22})
    ymax = float(pl.area_mm2.max()) * 1.6
    for j, t in enumerate(tools):
        c = cur[cur.tool == t].sort_values("t")
        a = axs[0, j]
        a.fill_between(c.t, c.y_lo, c.y_hi, color=COL[t], alpha=0.25, lw=0)
        a.plot(c.t, c.y, color=COL[t], lw=2)
        a.plot(pl.DayFloat, pl.area_mm2, "o", ms=3.5, color=INK, alpha=0.8, zorder=3)
        a.set_ylim(0, ymax); a.set_title(LAB[t], color=COL[t] if t != "V0" else "#2b8f98", weight="bold", fontsize=11)
        b = axs[1, j]
        b.fill_between(c.t, c.agr_lo, c.agr_hi, color=COL[t], alpha=0.25, lw=0)
        b.plot(c.t, c.agr, color=COL[t], lw=2)
        b.set_ylim(*[v for v in (-0.5 * 1.0, None)] if False else (None, None))
        b.set_xlabel("Day")
        if j == 0:
            a.set_ylabel("Rosette area (mm²)"); b.set_ylabel("Growth rate (mm²/day)")
        else:
            a.tick_params(labelleft=True);
    # common AGR limits from the medians (bands may run off the panel; stated in the caption)
    amax = float(cur.agr.max()) * 1.6
    for j in range(4):
        axs[1, j].set_ylim(0, amax)
    save(fig, "fig4_plant1003.svg")

if __name__ == "__main__":
    fig_sim(); fig_horizon(); fig_calib(); fig_plant(); fig_paired_gf(); fig_spike(); fig_arab_desc()
