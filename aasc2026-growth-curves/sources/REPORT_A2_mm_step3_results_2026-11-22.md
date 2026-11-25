# REPORT A2: Mechanistic-mean GP, Step 3 results (x50t, all conditions)    (22 Nov 2026, Claude Code as hal9000)

**SIMULATION ONLY.**
- Simulated truths are very smooth (dense design 14 readings; wheat 43).
- **The simulated "logistic" truth is a 3-parameter logistic plus 5 mm²**, so the logistic mean is not exact early in the season.
- **M-J ran at adapt_delta 0.99 against V0 at its own frozen setting (0.8);** M-2S as built at 0.8.
- These are labelled variants. **A pass would have earned a labelled variant beside the frozen arm; nothing changes on the poster (X15; Monday 23 Nov decision).**

**STOP for Controller review.**

## Verdict

| Variant | c1 coverage | c2 interval score | c3 flagged share | c4 spike degradation | c5 two-phase | c6 end AGR | c7 mid-gap | Non-logistic | **Go** |
|---|---|---|---|---|---|---|---|---|---|
| M-2S (0.8) | PASS | PASS | **FAIL** | **FAIL** (all 3 s) | PASS | PASS | PASS | PASS | **NO GO** |
| M-J (0.99) | PASS | PASS | PASS (narrow) | **FAIL** (all 3 s) | PASS | PASS | PASS | PASS | **NO GO** |

**Under the declared go rule, neither variant goes forward.** M-J's only failure is c4. A descriptive extra (section 4) shows what c4 is measuring there.

## Runs, stop conditions, md5s

- **No stop condition hit.**
- **Gate:** Andrew, 20 Nov 23:41. **PASSED**: 8 of 8 refits identical to the stored ad99 fits (log `gp_mm_step3_gate_andrewb_20261120_234113.log`, 9db5aca4…).
- **Clean stage:** Andrew, 20 Nov 23:44–23:53. 1,260 fits, 6 chunks (log `gp_mm_step3_clean_andrewb_20261120_234401.log`, c2204da8…).
- **Rest stage:** Andrew, 21 Nov 00:01–00:27. 4,020 fits, 22 chunks (log `gp_mm_step3_rest_andrewb_20261121_000106.log`, ab3f7273…).
- **Fitting cost:** about 13.3 core-hours over both stages.
- **Code:** both stages used the gate file and code md5 `07_gp-mm-step3.R` f1a286d6… (approval line normalised; the file on disk with the flag TRUE is 4a00b538…) and `07_gp-mm-step3-fns.R` 2e124d62….
- **Frozen files, Step 2 files and the 10 stored S1/S2 V0 chunks:** md5s identical at start and end of every run.
- **Summary:** Andrew, 22 Nov 19:35, `07_gp-mm-step3-summary.R` v0.2a (9b1f7952…). It read 28 Step 3 chunks made with one code version, plus the stored V0 (log `gp_mm_step3_summary_x50t_andrewb_20261122_193535.log`, b99525d3…).
- **Outputs:** in this folder (md5s at the end).
- **Before the results:** the Step 3 scripts, summary v0.2 (c6 decided rule, c4 check, known-answer tests) and the v0.2a bug fix were all reported.

**Scenario sets (X29 rule by label; test half, common scenarios):**
- c1 = dense 1, 4, 6–11, 14;
- c2, c3 = the 24 non-sparse scenarios;
- c4 = the 24 non-sparse scenarios, for each spike level;
- c5 = two_phase 4, 9, 14, 24, 39, 54.

Sparse 66 is in no go-rule set.

## 1. c1–c5 (B2's `gv2_go_rule`, unchanged)

| | M-2S | M-J |
|---|---|---|
| c1 paired Δ coverage, size / AGR / RGR (SE) | +0.003 (0.004) / −0.001 (0.005) / +0.003 (0.005): PASS | +0.004 (0.003) / −0.002 (0.004) / +0.002 (0.004): PASS |
| c2 median paired Δ relative interval score, size / AGR / RGR | −45 / −210 / −84: PASS | −40 / −169 / −74: PASS |
| c2 failure rate (clean test) | 0.4% (V0 0%) | 0% |
| c3 paired Δ flagged share (SE) | **+0.297 (0.036): FAIL** | +0.050 (0.026): PASS; the limit is 2 × 0.026 = 0.052 |
| c5 two_phase size coverage (variant / V0), failure rate | 0.948 / 0.967, 0%: PASS | 0.971 / 0.967, 0%: PASS |

## 2. c4: spike degradation, per variant against V0 (= the go-rule c4), with B2's all-arms version beside

**Median over paired datasets of is_y(spiked) / is_y(same clean dataset).** PASS requires the variant's median to be smaller than V0's at all three s.

| Pair | s | Datasets (dropped) | Variant ratio (bootstrap SE) | V0 ratio (SE) | Verdict |
|---|---|---|---|---|---|
| M-2S vs V0 | 0.05 | 205 (35) | 5.85 (0.19) | 3.25 (0.11) | FAIL |
| | 0.12 | 157 (83) | 8.68 (0.38) | 5.55 (0.23) | FAIL |
| | 0.30 | 123 (117) | 12.00 (0.79) | 9.40 (0.39) | FAIL |
| **M-J vs V0** | 0.05 | 240 (0) | **4.42** (0.19) | 3.05 (0.13) | FAIL |
| | 0.12 | 240 (0) | **6.34** (0.38) | 4.49 (0.24) | FAIL |
| | 0.30 | 240 (0) | **9.44** (0.39) | 6.98 (0.27) | FAIL |

- **M-2S's dropped datasets are its own prefit failures;** M-J dropped none.
- **B2's all-arms version** (a dataset is dropped if any arm failed: 205 / 157 / 123 datasets) gives M-J 4.60 / 7.41 / 10.29 against V0 3.25 / 5.55 / 9.40. The verdict is the same. That version is shown for information only, as decided (`mm3_c4_pairwise_and_allarms.csv`).

## 3. c6, c7 and the non-logistic requirement (clean test fits, both arms worked)

| | M-2S | M-J | V0 |
|---|---|---|---|
| c6 thin pooled median (% of peak AGR), all curves | 9.6 | 10.2 | 69.7 |
| c6, each thin design not worse | thin7, thin5, thin4 all yes | all yes | n/a |
| c6 dense: d = median(variant) − median(V0) (bootstrap SE, B = 2000, seed 20261121) | −7.4 (1.2) | −6.8 (1.1) | n/a |
| **c6** | **PASS** | **PASS** | |
| c7 thin5 + thin4 pooled median (% of K) | 3.8 | 4.0 | 16.2 / 15.9 |
| **c7** | **PASS** | **PASS** | |
| Non-logistic only: c6 thin, c7 | 14.9 vs 72.4; 4.7 vs 15.0: PASS, PASS | 14.5 vs 72.3; 3.1 vs 14.2: PASS, PASS | |
| Band width, variant / V0: c6 window, c7 gap (thin5 + thin4) | 0.30, 0.23 | 0.42, 0.49 | |

**The c6 and c7 gains come with narrower bands, not wider ones.** Both variants' dense medians are also below V0's (M-2S 5.7, M-J 6.3, V0 13.1), so the open question about reading the dense part of c6 doesn't change this result.

**By design**, c6 median (variant / V0):

| Design | M-2S | M-J |
|---|---|---|
| dense | 5.7 / 13.1 | 6.3 / 13.1 |
| thin7 | 7.5 / 29.2 | 8.4 / 29.2 |
| thin5 | 8.8 / 72.5 | 9.8 / 72.5 |
| thin4 | 14.9 / 91.1 | 15.3 / 91.1 |

c7 by design is in `mm3_c6_c7_by_design.csv`. On thin7 and dense, c7 is about level with V0 (M-J 2.0 against 1.9 on thin7; 2.8 against 2.3 on dense); on thin5 and thin4 it is much lower.

## 4. Descriptive extra on c4 (NOT part of the go rule; declared after the results were seen)

Log `gp_mm_step3_extra_hal9000_20261122_194010.log` (a6b346a3…). These are **absolute** relative interval scores (size, % of the true range): test half, non-sparse, pairs where both worked, medians.

| Condition | M-J | V0 | Median paired Δ (M-J − V0) | M-J lower in | M-2S | V0 | Δ | M-2S lower in |
|---|---|---|---|---|---|---|---|---|
| clean | 43.6 | 83.6 | −39.8 | 82% | 30.4 | 82.6 | −45.4 | 84% |
| spike_005 | 147.8 | 185.9 | −13.1 | 73% | 114.9 | 136.3 | −2.8 | 59% |
| spike_012 | 199.3 | 245.2 | −21.1 | 82% | 142.5 | 146.6 | −0.3 | 50% |
| spike_030 | 280.8 | 347.2 | −35.2 | 85% | 180.6 | 179.6 | −0.3 | 52% |

**Description:**
- **On spiked data, M-J's absolute interval score is lower than V0's in most pairs.**
- **c4 divides by each dataset's own clean score.** M-J's clean score is about half V0's, so the same or a smaller absolute spiked score gives a larger ratio.
- **The verdict stands as declared (FAIL).** Whether c4's ratio form is the right test for a variant that is much tighter on clean data is a question for the Controller, not something I've changed.
- **M-2S's spiked comparison is on its surviving datasets only** (many thin spiked datasets fail the prefit), so its rows are not like for like with M-J's.

## 5. Counts, failures, flags (all conditions; `mm3_counts.csv`)

**Failures** (each counts as a miss):
- **M-J:** none in 2,120 fits.
- **V0:** none.
- **M-2S:** 507 of 2,120. By design and condition, across both halves (noise conditions are test only):

| Condition | Dense | Sparse | thin4 | thin5 | thin7 |
|---|---|---|---|---|---|
| clean | 0 | 20 | 3 | 0 | 0 |
| noise_const | 5 | 0 | 0 | 0 | 0 |
| spike_005 | 0 | 0 | 52 | 16 | 4 |
| spike_012 | 0 | 0 | 83 | 48 | 35 |
| spike_030 | 0 | 0 | 97 | 80 | 64 |

All are prefit failures: "too few readings" for sparse, otherwise mostly "prefit at bound". The 5 noise_const failures are dense datasets with 11–12 of 14 readings left after clamping to zero, with the prefit at the lK or tm bound.

**Flagged share of worked fits, test half:**

| Arm | clean | noise_const | noise_prop | spike 0.05 / 0.12 / 0.30 |
|---|---|---|---|---|
| V0 | 10.4% | 4.4% | 0% | 5.0 / 5.8 / 6.2% |
| M-2S | **36.8%** | 10.6% | **54.4%** | 7.3 / 3.8 / 5.7% |
| M-J | 12.0% | 2.2% | **18.9%** | 0.4 / 0.4 / 0% |

M-J flags more than V0 under noise_prop (17 of 90 against 0). That condition isn't part of c3, which uses clean data, but it is worth knowing.

## 6. Per-scenario table (`mm3_per_scenario.csv`; test half, clean; 25 scenarios × 3 arms)

**Columns:** n, worked, failed, flagged, coverage of size / AGR / RGR (failures count as misses), median relative interval score (size, AGR), median band widths, and medians of c6, c6 width, c7 and c7 width. Highlights, as description:

- **Thin designs:** V0's size coverage is 1.000 on almost every thin scenario, with very wide bands. For example, its median relative interval score on scenario 53 is 452, against 68 for M-2S and 117 for M-J. The variants' narrower bands cost a little coverage in places: M-2S drops to 0.86–0.92 on 22, 39, 52, 53 and 55; M-J's lowest is 0.926 (scenario 22).
- **M-2S flags** are concentrated on dense low-noise scenarios (1 and 4: 7 of 10 each) and on thin4 (51–54: 7–9 of 10).
- **M-J flags** are at most 2 of 10 on every scenario except dense high-noise 11 (5 of 10).
- **Sparse 66:** M-2S fails 10 of 10 (by rule); V0 is flagged on 8 of 10; M-J works on 10 of 10 with none flagged.

## Files (this folder)

| File | md5 |
|---|---|
| `mm3_c1_c5.csv` | c4f3f25fddd9ac5c78cb52dc195e522d |
| `mm3_c4_pairwise_and_allarms.csv` | a2b35fcd04f342b05a3a041b76968eef |
| `mm3_c6_c7.csv` | 885d061acae9152f9a602269670b0a53 |
| `mm3_c6_c7_by_design.csv` | 233cf2bfddb9909e4ecb64e19c92429d |
| `mm3_per_scenario.csv` | 766eb4b34770fd044219115bd6d80921 |
| `mm3_counts.csv` | 68a7955703d4560dcf32bb65e14d87d8 |
| `mm3_interval_score_by_day.csv` | 6fc1bba0fe6d818165a53a6fe08b63cb |
| `mm3_fits_all.csv` | e9b60f59293c29e14df9eb33e169bc0f |
| `mm3_figures.pdf` | ea2b05f7cfd6feb1024dc709d2ad5a9e |

Step 3 chunks: `data/data-sim/m6_gp_mm/x50t/` (28 files; md5s in the summary log).

## Limits to carry forward

- **Simulated truths are very smooth;** dense has 14 readings against wheat's 43.
- **The logistic truth includes y0 = 5.**
- **One simulation screen,** 10 test replicates per scenario.
- **M-J ran at 0.99, V0 at 0.8.**
- **This addresses simulation findings only** (edges and thin data). It does nothing for wheat's blank early flights.
- **X15 stands:** no real-data GP statement.
- **The early mid-gap results** for these designs were already seen in the Step 1 diagnostic, so Step 3 isn't independent of that motivation. The judged replicates (121–130) are new data, though.

## Not checked

- **The c6 / c7 per-fit values** against an independent hand calculation (only the known-answer tests on synthetic inputs).
- **The bootstrap SE** against an independent implementation.
- **Any figure in `mm3_figures.pdf`** was not opened for this report.
- **Why M-J flags under noise_prop** was not examined.
