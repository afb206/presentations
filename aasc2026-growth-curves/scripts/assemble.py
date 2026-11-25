# Stage 3: fill template.html -> index.html.
#  * numbers in the text are {{tokens}} read from saved files (registered with put(), written to source_list.csv)
#  * tables and figures replace the slot <div class="slot"> that carries their token; slots with no data stay visible
import html, re, os
import numpy as np
import pandas as pd
import pyreadr
from common import *

def f2(x): return f"{x:.2f}"
def f3(x): return f"{x:.3f}"
def pc(x, d=0): return f"{100 * x:.{d}f}%"
def esc(s): return html.escape(str(s))

# ======================================================================================= load
F1a, F3, F4b = rd("F1a_error_each_tool.csv"), rd("F3_inclusion_width.csv"), rd("F4b_H3_error_by_day_log_pairing.csv")
G1, G2, G4 = rd("G1_cov_y.csv"), rd("G2_w_y.csv"), rd("G4_fail.csv")
P1, P2 = rd("P1_cov_y.csv"), rd("P2_w_y.csv")
M3 = rd("mm3_fits_all.csv")
M3 = M3[(M3.half == "test") & (M3.condition == "clean") & M3.variant.isin(["V0", "MJ"])]

# ======================================================================================= tokens (same rules as the poster)
for v, name in [("MJ", "mj"), ("V0", "v0")]:
    d = M3[(M3.variant == v) & (M3.design == "thin4")]
    put(f"c1_cov_{name}", f2(d.cov_y.mean()), "mm3_fits_all.csv", f"mean cov_y, variant {v}, clean, test half, design thin4 (n={len(d)})")
    put(f"c1_w_{name}", f"{d.w_y.median():.0f}", "mm3_fits_all.csv", f"median w_y (% of K), variant {v}, clean, test half, thin4")
r = one(P1, group="thin4", arm="log_nls")
put("c1_cov_ln", f2(r.est), "P1_cov_y.csv", "est, group thin4, arm log_nls (paired, on GP+L's datasets)")
put("c1_lab_ln", r.lab.split(",")[0].replace("/", " of "), "P1_cov_y.csv", "lab (first part: fits that worked / datasets), thin4, log_nls")
put("c1_wp_ln", f"{one(P2, group='thin4', arm='log_nls').est:.0f}", "P2_w_y.csv", "est (median band width, % of K), thin4, log_nls (paired)")
r = one(G4, design="thin4", arm="log_spline")
put("c1_fail_sp", r.lab.replace(" failed", "").replace("/", " of "), "G4_fail.csv", "lab, thin4, log_spline (older grid, other replicates)")

def e(level, tool): return one(F1a, level=level, tool=tool)
for lv in ["H1", "H2"]:
    put(f"c2_mj_{lv}", f3(e(lv, "MJ").est), "F1a_error_each_tool.csv", f"est, {lv}, MJ")
    put(f"c2_carry_{lv}", f3(e(lv, "carry").est), "F1a_error_each_tool.csv", f"est, {lv}, carry")
def b(group, tool): return one(F4b, group=group, tool=tool)
for t_, n in [("log_nls", "ln"), ("MJ", "mj"), ("V0", "v0"), ("spline_log", "sp")]:
    put(f"c2_{n}_d36", f2(b("Days 36-44", t_).est), "F4b_H3_error_by_day_log_pairing.csv", f"est, Days 36-44, {t_}")
for t_, n in [("spline_log", "sp"), ("log_nls", "ln")]:
    r = b("Days 36-44", t_)
    put(f"c2_f_{n}", f"{int(r.n_fail)} of {int(r.n)}", "F4b_H3_error_by_day_log_pairing.csv", f"n_fail of n, Days 36-44, {t_}")
def c(level, tool): return one(F3, level=level, tool=tool)
for t_, n in [("MJ", "mj"), ("log_nls", "ln"), ("V0", "v0"), ("spline_log", "sp"), ("natural_nls", "nn"), ("spline_natural", "sn")]:
    put(f"c3_inc_{n}", pc(c("H3", t_).incl), "F3_inclusion_width.csv", f"incl, H3, {t_}")
r = c("H3", "spline_natural")
put("c3_fail_sn", f"{int(r.n_fail)} of {int(r.n)}", "F3_inclusion_width.csv", "n_fail of n, H3, spline_natural")

# ======================================================================================= tables
def chip(t_):
    cls = CHIP.get(t_)
    return f'<span class="chip {cls}">{LAB[t_]}</span>' if cls else esc(LAB[t_])

def parse_lab(lab):                      # "49/50, 2% fail" -> (worked, n)
    a = lab.split(",")[0].split("/")
    return int(a[0]), int(a[1])

def cell(cov, width, worked, n):
    fail = n - worked
    if worked == 0:
        return f'<td class="num"><b>{cov:.2f}</b><span class="sub cellfail">no fit worked (0 of {n})</span></td>'
    fcls = ' cellfail' if fail else ''
    return f'<td class="num"><b>{cov:.2f}</b><span class="sub">{width:.0f}% of K<br><span class="{fcls.strip()}">{fail} of {n} failed</span></span></td>'

def tbl_sim():
    DD = DES + ["sparse"]
    rows = []
    for v in ["MJ", "V0"]:
        cs = []
        for dsg in DD:
            d = M3[(M3.variant == v) & (M3.design == dsg)]
            cs.append(cell(d.cov_y.mean(), d.w_y.median(), int(d.ok.sum()), len(d)))
        rows.append((chip(v), "mm3_fits_all", cs))
    for a in ["log_nls", "natural_nls"]:
        cs = []
        for dsg in DD:
            p1, p2 = one(P1, group=dsg, arm=a), one(P2, group=dsg, arm=a)
            w, n = parse_lab(p1.lab)
            cs.append(cell(p1.est, p2.est if w else 0, w, n))
        rows.append((chip(a), "P1 and P2 (paired)", cs))
    cs = []
    for dsg in DD:
        g1, g2 = one(G1, design=dsg, arm="log_spline"), one(G2, design=dsg, arm="log_spline")
        w, n = parse_lab(g1.lab)
        cs.append(cell(g1.est, g2.est if w else 0, w, n))
    rows.append((chip("spline_log"), "G1 and G2 (other replicates)", cs))
    nread = {dsg: int(M3[M3.design == dsg].n_read.median()) for dsg in DD}
    head = "".join(f'<th class="num">{dsg}<span class="sub">{nread[dsg]} readings</span></th>' for dsg in DD)
    body = "".join(f"<tr><td>{a}<span class='sub'>{s}</span></td>{''.join(cs)}</tr>" for a, s, cs in rows)
    return ('<div class="tblwrap"><table><caption style="text-align:left;font-weight:600;padding:.3rem 0">Coverage of the 95% band for the true size curve; beneath it the median band width (% of the plateau size K) and the number of datasets with no fit</caption>'
            f'<tr><th>Arm</th>{head}</tr>{body}</table></div>'
            '<p class="figcap">Coverage counts a failed fit as a miss. Width is among fits that worked. GP+L and GP: 50 datasets per thinned design and 90 at dense, 10 at sparse (the <em>n</em> is in each cell: failed of datasets). '
            'GF: refit on GP+L’s datasets, so paired (50 per design, 10 at sparse). SPL (log): the older grid, other replicates, unpaired (500 per design, 100 at sparse). Differences between cells smaller than about twice the Monte Carlo error should not be read. Simulated truths are very smooth.</p>')

def tbl_holdout():
    order = ["MJ", "V0", "log_nls", "natural_nls", "spline_log", "spline_natural", "carry", "loglin", "natural_improved"]
    names = {"H1": "H1: last reading left out", "H2": "H2: last two readings left out", "H3": "H3: everything after Day 30 left out"}
    out = []
    for lv in ["H1", "H2", "H3"]:
        rows = ""
        for t_ in order:
            f1 = F1a[(F1a.level == lv) & (F1a.tool == t_)]
            f3_ = F3[(F3.level == lv) & (F3.tool == t_)]
            assert len(f1) == 1
            r = f1.iloc[0]
            if len(f3_):
                assert int(f3_.iloc[0].n) == int(r.n) and int(f3_.iloc[0].n_fail) == int(r.n_fail)
                inc, wid = pc(f3_.iloc[0].incl), f2(f3_.iloc[0].wlog)
            else:
                inc = wid = "–"
            cls = ' class="extra"' if t_ == "natural_improved" else (' class="rule"' if t_ == "carry" else "")
            fcell = f'{int(r.n_fail)} of {int(r.n)}' + (' <span class="cellfail">●</span>' if r.n_fail else "")
            rows += (f'<tr{cls}><td>{chip(t_)}</td><td class="num">{r.est:.3f}<span class="sub">{r.lo:.3f} to {r.hi:.3f}</span></td>'
                     f'<td class="num">{fcell}<span class="sub">{int(r.n_plants)} plants scored</span></td><td class="num">{inc}</td><td class="num">{wid}</td></tr>')
        out.append(f'<h3 style="margin-top:1.2rem">{names[lv]}</h3><div class="tblwrap"><table><tr><th>Arm</th><th class="num">Median absolute log error</th>'
                   '<th class="num">Readings with no fit</th><th class="num">Inside 90% interval</th><th class="num">Width ln(upper / lower)</th></tr>' + rows + '</table></div>')
    out.append('<p class="figcap">Failed readings and the plants scored are different counts and are kept apart: “no fit” counts left-out readings; the line beneath counts plants that were scored. '
               'The pale rows are baselines (no interval) and the EXTRA arm, which is reported only and is not part of the main comparison. '
               'The two numbers under each error are the <em>lo</em> and <em>hi</em> columns of F1a_error_each_tool.csv, as drawn in the poster’s figure. '
               '<span class="verify">file does not say what interval these are; confirm before publishing</span></p>')
    return "".join(out)

def tbl_horizon():
    G = ["Day 31", "Day 34", "Days 36-44"]
    T_ = ["V0", "MJ", "log_nls", "spline_log", "carry", "loglin"]
    head = "".join(f'<th class="num">{g.replace("-", " to ")}</th>' for g in G)
    body = ""
    for t_ in T_:
        cs = ""
        for g in G:
            r = b(g, t_)
            cs += f'<td class="num">{r.est:.2f}<span class="sub">{int(r.n_fail)} of {int(r.n)} no fit</span></td>'
        body += f"<tr><td>{chip(t_)}</td>{cs}</tr>"
    return (f'<div class="tblwrap"><table><caption style="text-align:left;font-weight:600;padding:.3rem 0">H3 error by day after the cut, median absolute log error (per-plant median, paired readings)</caption>'
            f'<tr><th>Arm</th>{head}</tr>{body}</table></div>')

def fig_block(img, alt, caption):
    return f'<figure class="fig"><img src="figs/{img}" alt="{esc(alt)}"><figcaption>{caption}</figcaption></figure>'

def failnote_f3():
    L = ["H1", "H2", "H3"]; parts = []
    for t_ in ["spline_log", "spline_natural", "natural_nls", "natural_improved"]:
        d = F3[F3.tool == t_].set_index("level").reindex(L)
        parts.append(f"{LAB[t_]} " + "/".join(str(int(x)) for x in d.n_fail))
    n = F3[F3.tool == "MJ"].set_index("level").reindex(L).n
    return "; ".join(parts) + " of " + "/".join(str(int(x)) for x in n) + " left-out readings"

def failnote_f4b():
    G = ["Day 31", "Day 34", "Days 36-44"]; parts = []
    for t_ in ["spline_log", "log_nls"]:
        d = F4b[F4b.tool == t_].set_index("group").reindex(G)
        parts.append(f"{LAB[t_]} " + "/".join(str(int(x)) for x in d.n_fail))
    n = F4b[F4b.tool == "MJ"].set_index("group").reindex(G).n
    return "; ".join(parts) + " of " + "/".join(str(int(x)) for x in n)

HOLD = f'<div class="holdout"><b>Label.</b> {esc(HOLDOUT_LABEL)}</div>'

def blk_fig_sim():
    G4s = G4[G4.arm == "log_spline"].set_index("design").reindex(DES)
    sp = ", ".join(f"{d} {int(round(100 * r.est))}%" for d, r in G4s.iterrows())
    return fig_block("fig1_sim.svg", "Coverage and band width of the size band by sampling design, five arms",
        "<b>Figure 1.</b> Simulation only. Coverage (left) and median band width (right) of the 95% band for the true size curve, by sampling design. "
        "GP+L and GP: <code>mm3_fits_all.csv</code> (clean, test half, same datasets). GF: refit on GP+L’s datasets (paired). "
        f"SPL (log), dashed and hollow: the older grid, other replicates, unpaired; failed fits there: {sp}. Failed fits count as misses in coverage. "
        "Simulated truths are very smooth.")

def blk_fig_horizon():
    return (fig_block("fig2_horizon.svg", "Median absolute log error by day after the cut",
        "<b>Figure 2.</b> " + esc(HOLDOUT_LABEL) + " Bars are the lo and hi columns of the file. Failed readings (Day 31/34/36 to 44, H3): " + esc(failnote_f4b()) + ".")
        + tbl_horizon()
        + fig_block("fig3_calibration.svg", "Share of left-out readings inside the 90% interval against interval width",
        "<b>Figure 3.</b> " + esc(HOLDOUT_LABEL) + " Points run H1, H2, H3 (longer gap). Failed fits count as misses. Failed readings (H1/H2/H3): " + esc(failnote_f3()) + ". Plain natural-scale arms and the EXTRA arm are drawn faint."))

def blk_fig_plant():
    return fig_block("fig4_plant1003.svg", "One plant, four tools: size and growth rate",
        "<b>Figure 4.</b> One plant (Arabidopsis C24, full data; plant 1135FA1003). Illustration of how the tools differ, not a comparison of accuracy. "
        "Curves are tool medians for whole days 4 to 44 with their bands (the Day 1 to 3 extrapolation is left out); dots are the plant’s readings. "
        "Growth-rate bands that run off the top of a panel are cut off by the axis. Source: <code>abstract_1003_curves.csv</code>, <code>arabidopsis_rosette_area_valid.rds</code>.")

def blk_arab_valid():
    d = list(pyreadr.read_r(os.path.join(SRC, "arabidopsis_rosette_area_valid.rds")).values())[0]
    g = d.groupby(["group", "Plant ID"]).size().groupby("group")
    rows = ""
    for grp in ["moving", "stationary"]:
        s = g.get_group(grp)
        rows += (f'<tr><td>{grp}</td><td class="num">{len(s)}</td><td class="num">{int(s.sum())}</td>'
                 f'<td class="num">{int(s.min())} to {int(s.max())} (median {s.median():.0f})</td></tr>')
    tot = g.size().sum(), d.shape[0]
    rows += f'<tr class="rule"><td><b>all</b></td><td class="num"><b>{tot[0]}</b></td><td class="num"><b>{tot[1]}</b></td><td></td></tr>'
    put("arab_n_moving", str(int(g.size()["moving"])), "arabidopsis_rosette_area_valid.rds", "number of distinct Plant ID, group moving")
    return ('<div class="tblwrap"><table><caption style="text-align:left;font-weight:600;padding:.3rem 0">Plants and readings in the valid set (<code>arabidopsis_rosette_area_valid.rds</code>)</caption>'
            '<tr><th>Group</th><th class="num">Plants</th><th class="num">Readings</th><th class="num">Readings per plant</th></tr>' + rows + '</table></div>')

def blk_references():
    txt = open(os.path.join(SRC, "poster.qmd")).read()
    seg = txt[txt.index("{.references}") + len("{.references}"):txt.index("{.funding}")]
    items = [s.strip() for s in seg.replace("\\fontsize{14pt}{17pt}\\selectfont", "").replace("::::", "").replace(":::", "").split("·") if s.strip()]
    li = ""
    for it in items:
        it = esc(it)
        it = re.sub(r"doi:(\S+)", r'<a href="https://doi.org/\1">doi:\1</a>', it)
        li += f"<li>{it}</li>"
    return f'<ol class="tight">{li}</ol>'


def tbl_session_info():
    rows = [("R", "4.5.2 (2025-10-31)", "sessionInfo_hal9000.txt; Andrew’s own session (pasted in chat, 25 Nov 2026, not a saved file) agrees"),
            ("asreml", "4.2.0.382", "asreml_version_TODO.txt (quotes a saved log, line 119) and Andrew’s pasted packageVersion(\"asreml\")"),
            ("rstan", "2.32.7", "sessionInfo_hal9000.txt"), ("StanHeaders", "2.39.1", "sessionInfo_hal9000.txt"),
            ("brms", "2.23.0", "sessionInfo_hal9000.txt (listed as loaded; header lines of the file)"),
            ("posterior", "1.7.0", "sessionInfo_hal9000.txt"), ("nlme", "3.1.169", "sessionInfo_hal9000.txt"), ("MASS", "7.3.65", "sessionInfo_hal9000.txt")]
    body = "".join(f"<tr><td class='m'>{esc(a)}</td><td>{esc(b_)}</td><td>{esc(c_)}</td></tr>" for a, b_, c_ in rows)
    return ('<div class="tblwrap"><table><caption>Software versions</caption><tr><th>Package</th><th>Version</th><th>Where it is recorded</th></tr>' + body + '</table></div>'
            '<p class="note">The R session on the analysis server (hal9000) has no ASReml; the spline fits ran under Andrew’s own login. Platform of the server session: Ubuntu 24.04.3 LTS, BLAS/LAPACK 3.12.0. '
            'Whether every fit ran under exactly these library versions is not recorded in the saved files.</p>')

def tbl_step3():
    a = rd("mm3_c1_c5.csv").set_index("variant"); b = rd("mm3_c6_c7.csv")
    mj, ms = a.loc["MJ"], a.loc["M2S"]
    bm = b[(b.variant == "MJ") & (b.scope == "all curves")].iloc[0]; bs = b[(b.variant == "M2S") & (b.scope == "all curves")].iloc[0]
    f1 = lambda x: f"{x:.1f}"; s3 = lambda x: f"{x:+.3f}"
    R = [("Paired change in coverage vs GP (size / growth rate AGR / RGR), dense designs", "mean change over 90 paired datasets (9 dense scenarios)",
          f"{s3(mj.c1_dcov_y)} / {s3(mj.c1_dcov_agr)} / {s3(mj.c1_dcov_rgr)}", f"{s3(ms.c1_dcov_y)} / {s3(ms.c1_dcov_agr)} / {s3(ms.c1_dcov_rgr)}"),
         ("Median paired change in relative interval score (size / AGR / RGR)", "negative = lower score than GP",
          f"{mj.c2_dis_med_y:.0f} / {mj.c2_dis_med_agr:.0f} / {mj.c2_dis_med_rgr:.0f}", f"{ms.c2_dis_med_y:.0f} / {ms.c2_dis_med_agr:.0f} / {ms.c2_dis_med_rgr:.0f}"),
         ("Fit failure rate on clean test data", "", pc(mj.c2_fail_rate_v, 1), pc(ms.c2_fail_rate_v, 1)),
         ("Change in the share of fits flagged by the diagnostics rule vs GP", "", s3(mj.c3_dflag), s3(ms.c3_dflag)),
         ("Spike: interval score spiked / clean, s = 0.05 / 0.12 / 0.30", "median over datasets",
          f"{f1(mj.c4_ratio_spike_005)} / {f1(mj.c4_ratio_spike_012)} / {f1(mj.c4_ratio_spike_030)}", f"{f1(ms.c4_ratio_spike_005)} / {f1(ms.c4_ratio_spike_012)} / {f1(ms.c4_ratio_spike_030)}"),
         ("Spike: same ratio for GP on the same datasets", "", f"{f1(mj.c4_ratio0_spike_005)} / {f1(mj.c4_ratio0_spike_012)} / {f1(mj.c4_ratio0_spike_030)}", f"{f1(ms.c4_ratio0_spike_005)} / {f1(ms.c4_ratio0_spike_012)} / {f1(ms.c4_ratio0_spike_030)}"),
         ("Two-phase curves: coverage of size, variant / GP; failure rate", "", f"{mj.c5_cov_v:.3f} / {mj.c5_cov_0:.3f}; {pc(mj.c5_fail_v)}", f"{ms.c5_cov_v:.3f} / {ms.c5_cov_0:.3f}; {pc(ms.c5_fail_v)}"),
         ("End-of-series AGR error, thin designs (median, % of peak AGR), variant / GP", "", f"{bm.c6_thin_v:.1f} / {bm.c6_thin_0:.1f}", f"{bs.c6_thin_v:.1f} / {bs.c6_thin_0:.1f}"),
         ("Mid-gap size error, thin5 + thin4 (median, % of K), variant / GP", "", f"{bm.c7_v:.1f} / {bm.c7_0:.1f}", f"{bs.c7_v:.1f} / {bs.c7_0:.1f}"),
         ("Band width, variant / GP: end window; mid-gap", "", f"{bm.c6_bw_ratio_all:.2f}; {bm.c7_bw_ratio_thin54:.2f}", f"{bs.c6_bw_ratio_all:.2f}; {bs.c7_bw_ratio_thin54:.2f}")]
    body = "".join(f"<tr><td>{esc(a_)}<br><small>{esc(n_)}</small></td><td class='num'>{esc(x)}</td><td class='num'>{esc(y)}</td></tr>" for a_, n_, x, y in R)
    return ('<div class="tblwrap"><table><caption>Step 3 quantities, simulation only (x50t run, clean and spiked test data)</caption>'
            '<tr><th>Quantity</th><th>GP+L (joint)</th><th>GP+L two-stage (not one of the five arms)</th></tr>' + body + '</table></div>'
            '<p class="note">Read from <code>mm3_c1_c5.csv</code> and <code>mm3_c6_c7.csv</code> (all curves). The files also carry a verdict column for each rule, declared before the runs; those words are not reproduced here. '
            'The spike ratio is larger for GP+L than for GP, but the Step 3 report (<code>REPORT_A2_mm_step3_results_2026-11-22.md</code>) notes that GP+L’s own clean score is about half GP’s, so the same or a smaller spiked score gives a larger ratio: a ratio of scores is not a direct measure of the spiked score. Smooth simulated truths; spikes hit only the earliest reading.</p>')

def _lead(x):                           # "0.0887 (0.0353 to 0.204)" -> 0.0887 ; "1.04e+03 (..)" -> 1040
    return float(str(x).split(" ")[0])

def tbl_arab_flags():
    d = rd("arab_mj_set_per_fit.csv")
    rows = []
    for arm, lab in [("V0", "GP (stored fits, adapt_delta 0.8)"), ("MJ", "GP+L (adapt_delta 0.99)")]:
        g = d[d.arm == arm]
        r44 = g.d44_ratio.dropna()
        rows.append((lab, len(g), int(g.ok.sum()), f"{int(g.flagged.sum())} ({pc(g.flagged.mean(), 1)})", int((g.divergent > 0).sum()), int((g.treedepth_hits > 0).sum()),
                     f"{g.max_rhat.median():.3f} ({g.max_rhat.max():.2f})", f"{g.min_ess_bulk.median():.0f} ({g.min_ess_bulk.min():.1f})",
                     f"{r44.median():.3f}; {pc((r44 > 1).mean(), 1)} of {len(r44)}"))
    body = "".join("<tr>" + "".join(f"<td>{esc(x)}</td>" for x in r) + "</tr>" for r in rows)
    return ('<div class="tblwrap"><table><caption>Arabidopsis C24, moving plants: diagnostics of the two Gaussian-process fits, description only</caption>'
            '<tr><th>Arm</th><th>Plants</th><th>Fits that worked</th><th>Flagged fits</th><th>Fits with any divergence</th><th>Fits with tree-depth hits</th>'
            '<th>Max R-hat: median (largest)</th><th>Min bulk ESS: median (smallest)</th><th>Fitted size / Day 44 reading: median; share above 1</th></tr>' + body + '</table></div>'
            '<p class="note">From <code>arab_mj_set_per_fit.csv</code>. Real data with no known truth, so nothing here is accuracy. Flagged fits are counted as worked; failures would count as misses (there were none). '
            'Fit-against-reading comparisons are in-sample.</p>')

def tbl_two_mode():
    d = rd("arab_mj_diag_8chain_vs_original.csv"); c = rd("arab_mj_diag_per_chain.csv")
    c["early"] = c.tm_med < 12
    body = ""
    for _, r in d.iterrows():
        cc = c[(c.plant_id == r.plant_id) & (c.fit.str.startswith("original" if r.fit.startswith("original") else "8 chains"))]
        body += (f"<tr><td class='m'>{esc(r.plant_id)}</td><td>{esc(r.fit)}</td><td>{int((~cc.early).sum())} / {int(cc.early.sum())}</td><td>{'yes' if r.flagged else 'no'}</td>"
                 f"<td class='num'>{r.max_rhat:.2f}</td><td class='num'>{r.min_ess_bulk:.0f}</td><td class='num'>{int(r.divergent)}</td><td class='num'>{int(r.treedepth_hits)}</td>"
                 f"<td class='num'>{r.fit_grid_d44:,.0f}</td><td class='num'>{r.band_width_d44:,.0f}</td><td class='num'>{r.day44_reading:,.0f}</td></tr>")
    return ('<div class="tblwrap"><table><caption>Four plants whose GP+L chains settle in two places (diagnostic, description only)</caption>'
            '<tr><th>Plant</th><th>Fit</th><th>Chains: main growth / early step</th><th>Flagged</th><th>Max R-hat</th><th>Min bulk ESS</th><th>Divergences</th><th>Tree-depth hits</th>'
            '<th>Fit, Day 44 (mm²)</th><th>Band width, Day 44 (mm²)</th><th>Day 44 reading (mm²)</th></tr>' + body + '</table></div>'
            '<p class="note">From <code>arab_mj_diag_8chain_vs_original.csv</code> and <code>arab_mj_diag_per_chain.csv</code>. A chain is counted as “early step” when its logistic midpoint is before Day 12; '
            'the other chains put the midpoint between Day 21 and Day 31. More chains did not settle the split. The fitted median curve at Day 44 barely changes between the original and the 8-chain fit; the band and the noise estimate do change when the chains split. '
            'Why the chains split is a hypothesis that has not been tested: all four plants have a zero reading early on and a steep Day 6 to Day 8 step, so a logistic on that early rise may fit about as well as one on the main rise.</p>')

def tbl_improved():
    rows = []
    for L in ["H1", "H2", "H3"]:
        a = rd(f"EXTRA_holdout_{L}_five_tools_summaryA_failures_as_misses.csv"); a = a[a.group == "all hidden readings"].iloc[0]
        fc = rd(f"EXTRA_holdout_{L}_failure_counts.csv").set_index("tool")
        for t_, key in [("natural_nls", "natural_nls"), ("natural_improved", "natural_improved"), ("log_nls", "log_nls"), ("MJ", "MJ")]:
            rows.append((L, LAB[t_] if t_ != "natural_improved" else "GF (nat), improved (EXTRA)", a[f"inclusion_{key}"], f"{_lead(a[f'abs_log_err_{key}_all5']):.3f}", f"{_lead(a[f'width_mm2_{key}_all5']):,.0f}",
                         f"{int(fc.loc[t_, 'plants_with_a_failure'])} plants ({int(fc.loc[t_, 'readings_failed'])} readings)"))
    body = "".join("<tr>" + "".join(f"<td>{esc(x)}</td>" for x in r) + "</tr>" for r in rows)
    return ('<div class="tblwrap"><table><caption>Extra arm <code>natural_improved</code> (GF on the natural scale with a variance that can grow with size), all left-out readings</caption>'
            '<tr><th>Level</th><th>Arm</th><th>Inside the 90% interval (failures as misses)</th><th>Median absolute log error</th><th>Median interval width (mm²)</th><th>Fits that failed</th></tr>' + body + '</table></div>'
            '<p class="note">' + esc(HOLDOUT_LABEL) + ' Post hoc arm, proposed after the H1 to H3 results; not part of the declared design. Error and width are over the readings where all five arms worked, so declared-arm figures here can differ slightly from the main table (for example H3 GF (log) error). '
            'Read from <code>EXTRA_holdout_H*_five_tools_summaryA_failures_as_misses.csv</code> and <code>EXTRA_holdout_H*_failure_counts.csv</code>.</p>'
            '<p><b>What it showed.</b> Letting the noise grow with size widened the natural-scale intervals and raised the share inside them, but it stayed far below 90%. The log-scale GF arm reaches a much higher share with wider intervals. '
            'The constant residual sd is therefore at most part of the reason the plain natural-scale arm’s intervals are narrow. Other possible causes were listed in the report and were not tested.</p>'
            '<p><b>Corrections carried from the report’s erratum (23 Nov 2026).</b> (1) The statement that the log-scale intervals are about 8 times wider than natural_improved’s holds only for H1; from the median widths the ratio is about 8 for H1, about 10 for H2 and about 4 for H3. '
            '(2) Most of natural_improved’s failures (4 of 7 plants in H1, 7 of 10 in H2, 39 of 69 in H3) come from the frozen fallback rule, where a chosen plain-nls fit counts as a failure under this arm’s declared rule, so the failure rate reflects that rule as much as the variance function. '
            'The source report is <code>REPORT_A3_EXTRA_natural_improved_H1_H2_H3_2026-11-23.md</code>.</p>')

def _need(fname, *strings):             # hand-keyed numbers must appear in the named file
    txt = open(os.path.join(SRC, fname), encoding="utf-8").read()
    for x in strings:
        assert x in txt, (fname, x)

def _tbl(caption, head, rows, note=""):
    th = "".join(f"<th>{esc(h)}</th>" for h in head)
    body = "".join("<tr>" + "".join(f"<td>{x}</td>" for x in r) + "</tr>" for r in rows)
    return f'<div class="tblwrap"><table><caption>{caption}</caption><tr>{th}</tr>{body}</table></div>' + (f'<p class="note">{note}</p>' if note else "")

def tbl_arab_cleaning():
    _need("README_arabidopsis_status.txt", "- 1 row", "51 plants", "Day 29", "Day 41", "8 plants", "485 raw Plant IDs", "425 plants", "384", "100")
    rows = [("Raw sheet, plant rows", "485 plant IDs", "one leaked label row (<code>1135FA</code>) removed first (1 row)"),
            ("Dead or failed germination", "51 plants removed", "peak area below 100 mm²"),
            ("Bad capture dates", "Day 29 and Day 41 removed for all plants", "within-plant crash rates well above neighbouring days"),
            ("Stunted plants", "8 plants removed", "peak area below 500 mm², at a natural gap in the sorted peak areas"),
            ("After cleaning", "425 plants", "see the table above for moving and stationary counts, read from the saved cleaned file")]
    return _tbl("Arabidopsis cleaning steps, in order", ["Step", "Count", "Rule"], rows,
                "From <code>README_arabidopsis_status.txt</code> (a status note, not a count table: no table of plants per stage was saved). The moving and stationary groups are reconstructed from the sampling pattern, not recorded; "
                "the note gives 384 moving and 100 stationary before cleaning and says this matches the paper. Days 38 and 40 are absent and Day 43 has far fewer plants, so pooled summaries on those days reflect sample size.")

def tbl_wheat_desc():
    _need("README_data_provenance.md", "15 trait files", "14 design files", "142,646 rows", "137,762 rows", "FPWW028", "5 of 7 TLS dates", "roughly 9-12x")
    _need("README_height_data_status.txt", "lodging NOT yet filtered")
    rows = [("Source", "FIP 1.0, ETH Zürich; plant height from the CroPyDB trait export (not the image-bundled copy)"),
            ("Files used", "15 trait files and 14 design files, joined by plot"),
            ("Rows", "142,646 in the input file; 137,762 in the cleaned file"),
            ("Exclusions", "lot FPWW028 (small 2019 lot) removed; laser-scan heights dropped on days with a drone measurement in lot FPWW022_lot3 (5 of 7 dates); six FPWW012 plants removed by a provisional rule"),
            ("Canopy cover", "not used (senescence artefact)"),
            ("Replication", "check plots are measured roughly 9 to 12 times as often as other genotypes, from earlier exploration; no count table was saved"),
            ("Lodging", "not yet filtered; detection is not finalised, so the wheat data are described only")]
    return _tbl("Wheat height data, as handled", ["Item", "What the notes record"], rows,
                "From <code>README_data_provenance.md</code> (wheat section) and <code>README_height_data_status.txt</code>. The status note carries a placeholder date. Nothing is fitted to the wheat data on this page.")

def tbl_maize_desc():
    _need("README_data_provenance.md", "1,752 accessions", "5 timepoints (37, 44, 60, 73, 115 DAP)", "8,760 rows", "138,329 rows", "NSL437884", "79 daily values")
    rows = [("Source", "Ames panel NDVI from drone flights (Wang et al., Dryad); published P-spline parameters and daily predictions are kept as a benchmark"),
            ("Accessions and time points", "1,752 accessions; 5 time points (37, 44, 60, 73 and 115 days after planting)"),
            ("Rows", "8,760 NDVI rows; 1,752 parameter rows; 138,329 daily-prediction rows (79 days per genotype)"),
            ("One odd accession", "NSL437884 has the text “FALSE” where its predictions should be: no benchmark predictions, but its NDVI readings and parameters are kept"),
            ("What NDVI is", "a bounded greenness index, not a size, so growth rates from it are not true AGR or RGR")]
    return _tbl("Maize NDVI panel, as handled", ["Item", "What the notes record"], rows,
                "From <code>README_data_provenance.md</code> (maize section). No table of genotypes and time points was saved. The prediction rows are consistent with 1,751 genotypes × 79 days (138,329), which is my arithmetic from those two numbers, not a saved count. Described only.")

def tbl_limits():
    rows = [("The hold-out is one genotype and one dataset, mostly near the plateau; a trivial baseline is as close as GP+L there", "Both are stated; H3 is the test before the plateau"),
            ("GP+L’s plateau prior is centred on the largest remaining reading, so it runs low after a cut well before the plateau (H3)", "Documented limit. Not tuned on these data; any variant is declared first and tested in simulation"),
            ("GP here is the frozen simulation-settings Gaussian process with no mean function", "Labelled as such on every hold-out result"),
            ("Spline and natural-scale growth-function fits fail on some plants and readings", "Failure counts sit beside every error; failures count as misses for coverage; failures are not patched"),
            ("GP+L chains settle in two places for some plants (four known; others not examined)", "The fitted curve is stable but the band and noise estimate depend on the mode"),
            ("Simulated truths are very smooth; the simulated logistic is logistic + 5 mm²; spikes hit only the earliest reading", "Stated with every simulation result"),
            ("A spike-cleaning rule (SC1) works only for blatant spikes on dense and thin7 designs", "Tested in simulation only; not used on real data"),
            ("Wheat lodging detection is not finalised", "The wheat look needs a saved list of cut points; not run"),
            ("Unexplained: <code>log_improved</code> equals <code>log_nls</code>; the growth-function “published” arms equal <code>log_nls</code> in about half of fits in an earlier grid; GP+L flags more fits under proportional noise; spline fits hit singularities", "Reported as unexplained"),
            ("The saved software record is the server’s <code>sessionInfo</code>; ASReml ran under Andrew’s own login", "Whether every fit used exactly these library versions is not recorded")]
    return _tbl("Limits and open items", ["Limit or open item", "How it is handled here"], [(a, b) for a, b in rows],
                "From section 7 of the project state file <code>CONTROLLER_STATE_2026-11-23_rev6a.md</code> (a project document, not copied here), with the last row added from the export.")

def scripts_links():
    L = [l.split() for l in open(os.path.join(SRC, "script_md5_methods.txt"))]
    return ("the method files used are listed here by name and the first eight characters of their md5; they are not published in this repository: " +
            "; ".join(f"<code>{esc(n)}</code> {m[:8]}…" for m, n in L) +
            '. Full md5s are in <a href="sources/script_md5_methods.txt"><code>sources/script_md5_methods.txt</code></a>')

# ======================================================================================= slots
FILLED = {
    "fig_tradeoff": (fig_block("fig1_sim.svg", "Coverage and band width by sampling design, five arms",
                       "<b>Simulation.</b> Coverage of the true size curve (left) and median band width (right) as readings become fewer. Smooth simulated truths; failed fits count as misses. Details in section 6.")
                     + fig_block("fig3_calibration.svg", "Share of left-out readings inside the 90% interval against interval width",
                       "<b>Left-out readings.</b> Share of left-out readings inside the 90% interval, against interval width, for H1 to H3. " + esc(HOLDOUT_LABEL) + " Details in section 8.")),
    "fig_cov_by_design": blk_fig_sim(),
    "tbl_sim_cov_width": tbl_sim(),
    "fig_holdout_horizon": blk_fig_horizon(),
    "tbl_holdout": tbl_holdout(),
    "tbl_session_info": tbl_session_info(),
    "tbl_step3": tbl_step3(),
    "tbl_arab_flags": tbl_arab_flags(),
    "fig_arab_desc": fig_block("fig7_arab_desc.svg", "Fitted size over the Day 44 reading, and the GP+L logistic midpoint, moving plants",
        "<b>Fit against the Day 44 reading (left) and GP+L logistic midpoint (right).</b> Arabidopsis C24, moving plants, from <code>arab_mj_set_per_fit.csv</code>. Description only, in-sample: no truth is known, so this is not accuracy. The few plants at an early midpoint are discussed below."),
    "fig_two_mode": tbl_two_mode(),
    "tbl_improved": tbl_improved(),
    "tbl_arab_cleaning": tbl_arab_cleaning(),
    "tbl_wheat_desc": tbl_wheat_desc(),
    "tbl_maize_desc": tbl_maize_desc(),
    "tbl_limits": tbl_limits(),
    "fig_paired_gf": fig_block("fig5_paired_gf.svg", "Coverage of the true size curve for GP+L and the growth-function arms, by sampling design",
        "<b>Growth-function arms refit on GP+L’s datasets.</b> Coverage of the true size curve by design, from <code>P1_cov_y.csv</code> (paired refit). Open markers: no fit worked in that design (0 of the datasets). Dashed line: 0.95. Simulation only: smooth truths, logistic = logistic + 5 mm², failed fits count as misses."),
    "fig_spike": fig_block("fig6_spike.svg", "Interval score ratio, spiked over clean, for GP and GP+L, by spike size",
        "<b>Early spike.</b> Median ratio of the interval score on spiked data to the same dataset clean (<code>mm3_c1_c5.csv</code>, columns <code>c4_ratio</code> and <code>c4_ratio0</code>; 240 datasets each). The spike replaces the earliest reading only. Simulation only; smooth simulated truths."),
}
def replace_slot(htm, token, content):
    key = "{{" + token + "}}"
    pos = htm.find(key)
    if pos < 0:
        print("slot token not found:", token); return htm
    start = htm.rfind('<div class="slot">', 0, pos)
    i, depth = start, 0
    for m in re.finditer(r"<div\b|</div>", htm[start:]):
        depth += 1 if m.group() == "<div" else -1
        if depth == 0:
            end = start + m.end(); break
    return htm[:start] + content + htm[end:]

htm = open(os.path.join(HERE, "template.html")).read()
for tok, content in FILLED.items():
    htm = replace_slot(htm, tok, content)
# the Home trade-off block should not carry a duplicated Figure 1 caption style; fine as is.
htm = htm.replace("{{fig_plant}}", blk_fig_plant()).replace("{{tbl_arab_valid}}", blk_arab_valid())
htm = htm.replace("{{references_list}}", blk_references()).replace("{{scripts_links}}", scripts_links())
htm = htm.replace("{{poster_pdf_url}}", "assets/poster_A0_DRAFT21.pdf")
htm = htm.replace("Poster (PDF)</a>", "Poster (PDF, draft 21)</a>")

# ---- the number-to-file section (built last, so it lists every token actually used)
need = set(re.findall(r"\{\{(\w+)\}\}", htm))
missing = sorted(t_ for t_ in need if t_ not in T and not re.match(r"^(tbl|fig|scripts)_|sim_|holdout_", t_) and t_ not in
                 {"references_list", "scripts_links", "tbl_source_list"})
print("tokens left unfilled in the text (not slots):", missing)
used = sorted(t_ for t_ in need if t_ in T)
for k in used:
    htm = htm.replace("{{" + k + "}}", T[k][0])
rows = [(k, T[k][0], T[k][1], md5_of(T[k][1]), T[k][2]) for k in sorted(T) if ("{{" + k + "}}") in open(os.path.join(HERE, "template.html")).read() or k in ("arab_n_moving",)]
pd.DataFrame(rows, columns=["token", "value", "file", "md5", "read_as"]).to_csv(os.path.join(HERE, "source_list.csv"), index=False)

def src_tables():
    r1 = "".join(f'<tr><td class="m">{esc(a)}</td><td>{esc(b_)}</td><td class="m">{esc(c_)}</td><td class="m">{d_}</td><td>{esc(e_)}</td></tr>' for a, b_, c_, d_, e_ in
                 [(r[0], r[1], r[2], r[3], r[4]) for r in rows])
    figs = [("Figure 1, Home", "scripts/figures.py: fig_sim", "mm3_fits_all.csv, P1_cov_y.csv, P2_w_y.csv, G1_cov_y.csv, G2_w_y.csv"),
            ("Figure 2", "scripts/figures.py: fig_horizon", "F4b_H3_error_by_day_log_pairing.csv"),
            ("Figure 3, Home", "scripts/figures.py: fig_calib", "F3_inclusion_width.csv"),
            ("Figure 4", "scripts/figures.py: fig_plant", "abstract_1003_curves.csv, arabidopsis_rosette_area_valid.rds"),
            ("Table, simulation", "scripts/assemble.py: tbl_sim", "mm3_fits_all.csv (mean cov_y, median w_y, count of ok), P1_cov_y.csv, P2_w_y.csv, G1_cov_y.csv, G2_w_y.csv (read as stored; lab gives fits that worked)"),
            ("Tables, left-out readings", "scripts/assemble.py: tbl_holdout, tbl_horizon", "F1a_error_each_tool.csv, F3_inclusion_width.csv (n and n_fail checked equal), F4b_H3_error_by_day_log_pairing.csv"),
            ("Table, Arabidopsis plants", "scripts/assemble.py: blk_arab_valid", "arabidopsis_rosette_area_valid.rds (counts of distinct Plant ID and rows by group)"),
            ("Step 3 table, spike figure", "scripts/assemble.py: tbl_step3; scripts/figures.py: fig_spike", "mm3_c1_c5.csv, mm3_c6_c7.csv"),
            ("Paired growth-function figure", "scripts/figures.py: fig_paired_gf", "P1_cov_y.csv"),
            ("Arabidopsis tables and figure", "scripts/assemble.py: tbl_arab_flags, tbl_two_mode; scripts/figures.py: fig_arab_desc", "arab_mj_set_per_fit.csv, arab_mj_diag_8chain_vs_original.csv, arab_mj_diag_per_chain.csv"),
            ("Extra arm table", "scripts/assemble.py: tbl_improved", "EXTRA_holdout_H1-H3 summaryA and failure_counts files; REPORT_A3_EXTRA_natural_improved_H1_H2_H3_2026-11-23.md (erratum text)"),
            ("Data, limits, software tables", "scripts/assemble.py: tbl_arab_cleaning, tbl_wheat_desc, tbl_maize_desc, tbl_limits, tbl_session_info", "README_data_provenance.md, README_arabidopsis_status.txt, README_height_data_status.txt, sessionInfo_hal9000.txt, asreml_version_TODO.txt (numbers are hand-keyed and checked to appear in the file)"),
            ("References", "scripts/assemble.py: blk_references", "poster.qmd (references block)")]
    r2 = ""
    for a, s, fl in figs:
        mds = "<br>".join(f"{esc(x.strip().split(' ')[0])} {md5_of(x.strip().split(' ')[0])}" for x in fl.split(",") if x.strip().split(" ")[0] in set(USED) | {"poster.qmd"} and True)
        r2 += f"<tr><td>{esc(a)}</td><td class='m'>{esc(s)}</td><td>{esc(fl)}</td><td class='m'>{mds}</td></tr>"
    r3 = "".join(f'<tr><td class="m">{esc(l.split()[1])}</td><td class="m">{l.split()[0]}</td></tr>' for l in open(os.path.join(SRC, "MD5SUMS.txt")))
    return (f'<h3>Numbers printed in the text</h3><div class="tblwrap"><table class="src"><tr><th>Token</th><th>Value</th><th>File</th><th>md5</th><th>How it was read</th></tr>{r1}</table></div>'
            f'<h3>Figures and tables</h3><div class="tblwrap"><table class="src"><tr><th>Item</th><th>Script</th><th>Inputs and what was taken</th><th>Input md5</th></tr>{r2}</table></div>'
            f'<h3>All input files</h3><div class="tblwrap"><table class="src"><tr><th>File (in <code>sources/</code>)</th><th>md5</th></tr>{r3}</table></div>'
            '<p>The same lists are in <a href="source_list.csv"><code>source_list.csv</code></a> and <a href="sources/MD5SUMS.txt"><code>sources/MD5SUMS.txt</code></a>. '
            'Every file above is a copy of a file the poster used, and each md5 was checked against the poster’s own list when it was copied.</p>')
htm = replace_slot(htm, "tbl_source_list", src_tables())
open(os.path.join(HERE, "index.html"), "w").write(htm)
left = sorted(set(re.findall(r"\{\{(\w+)\}\}", htm)))
print("written index.html; slots still empty:", left)
