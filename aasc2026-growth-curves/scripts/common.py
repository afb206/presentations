# Shared helpers for the report build. Reads saved files in sources/, computes nothing beyond the
# summaries named where they are used. Every number that reaches the page is registered with put()
# (token, value as printed, file, md5, how it was read) and written to source_list.csv.
import hashlib, os, shutil, sys
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(HERE, "sources")
POSTER_SOURCES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "sources")  # pass the poster sources/ folder to refresh the copies

# the files this build uses (all are copies of files the poster used; see the poster's MD5SUMS.txt)
USED = ["F1a_error_each_tool.csv", "F3_inclusion_width.csv", "F4b_H3_error_by_day_log_pairing.csv",
        "G1_cov_y.csv", "G2_w_y.csv", "G4_fail.csv", "P1_cov_y.csv", "P2_w_y.csv",
        "mm3_fits_all.csv", "abstract_1003_curves.csv", "abstract_1003_sources.txt",
        "arabidopsis_rosette_area_valid.rds"]

def md5(path):
    return hashlib.md5(open(path, "rb").read()).hexdigest()

def md5_of(fname):
    return md5(os.path.join(SRC, fname))

def rd(fname):
    return pd.read_csv(os.path.join(SRC, fname), low_memory=False)

def one(df, **q):
    for k, v in q.items():
        df = df[df[k] == v]
    assert len(df) == 1, (q, len(df))
    return df.iloc[0]

# ---- the registry of numbers printed in the text ---------------------------------------------
T = {}
def put(tok, val, file, how):
    T[tok] = (val, file, how)

# ---- arm names, colours (palette A, same as scripts/figstyle.py of the poster) ---------------
COL = {"MJ": "#00808B", "V0": "#6CC3CB", "log_nls": "#414042", "natural_nls": "#9A9A9A",
       "spline_log": "#EC8525", "spline_natural": "#F4B67A", "carry": "#9A9A9A", "loglin": "#B9B9B9",
       "natural_improved": "#9A9A9A"}
LAB = {"MJ": "GP+L", "V0": "GP", "log_nls": "GF (log)", "natural_nls": "GF (nat)",
       "spline_log": "SPL (log)", "spline_natural": "SPL (nat)", "carry": "Repeat last reading",
       "loglin": "Log-linear baseline", "natural_improved": "GF (nat), improved (EXTRA)"}
CHIP = {"MJ": "c-gpl", "V0": "c-gp", "log_nls": "c-gflog", "natural_nls": "c-gfnat",
        "spline_log": "c-spl", "spline_natural": "c-spl"}
DES = ["dense", "thin7", "thin5", "thin4"]

HOLDOUT_LABEL = ("One genotype (Arabidopsis C24, moving plants). GP = frozen, simulation settings, no mean "
                 "function. Scored on left-out readings only; readings from one plant are not independent.")

if __name__ == "__main__":
    # Stage 1: copy the used files from the poster's sources/ and check each md5 against its MD5SUMS.txt
    os.makedirs(SRC, exist_ok=True)
    listed = {}
    for line in open(os.path.join(POSTER_SOURCES, "MD5SUMS.txt")):
        parts = line.split()
        if len(parts) == 2:
            listed[parts[1]] = parts[0]
    bad = 0
    for f in USED:
        shutil.copy2(os.path.join(POSTER_SOURCES, f), os.path.join(SRC, f))
        got = md5_of(f)
        want = listed.get(f)
        status = "match" if got == want else ("NOT LISTED" if want is None else "MISMATCH")
        bad += status != "match"
        print(f"{status:9s} {got}  {f}")
    with open(os.path.join(SRC, "MD5SUMS.txt"), "w") as fh:
        for f in USED:
            fh.write(f"{md5_of(f)}  {f}\n")
    print("copied", len(USED), "files; problems:", bad)
