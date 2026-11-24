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
    seg = txt[txt.index("{.references}"):txt.index("{.funding}")]
    items = [s.strip() for s in seg.replace("\\fontsize{14pt}{17pt}\\selectfont", "").replace("::::", "").replace(":::", "").split("·") if s.strip()]
    li = ""
    for it in items:
        it = esc(it)
        it = re.sub(r"doi:(\S+)", r'<a href="https://doi.org/\1">doi:\1</a>', it)
        li += f"<li>{it}</li>"
    return f'<ol class="tight">{li}</ol>'

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
htm = htm.replace("{{references_list}}", blk_references())
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
