# REPORT (EXTRA): natural_improved (varPower) beside V0, M-J, natural_nls and log_nls, levels H1, H2, H3
23 Nov 2026. Session A3, Claude Code as hal9000. Runs by Andrew on hal9000 (R 4.5.2, OMP_NUM_THREADS=1, 4 workers).

EXTRA ARM UNDER X66, PROPOSED AFTER THE H1-H3 RESULTS (POST-HOC): a test of A3's hypothesis that the plain natural_nls arm's narrow
intervals come from its constant residual sd on the natural scale. One form, run once per level. NOT part of the declared design; the
declared four-tool tables and reports are unchanged.
REAL DATA: Arabidopsis C24, moving plants. HOLD-OUT CHECK UNDER X60: predictive error and interval inclusion for the hidden readings
only; everything else is description (D157). X15 withdrawn for this check (X61, provisional).

## Provenance
- Declaration: `DECLARATION_natural_improved_2026-11-23.md` (ac00db0f680bfbab65c70b8fa054ebaf), approved 02:00 (X66).
- Code: `code/08_paine-holdout.R` (c55acd8bf9ad83216b22852e4a02c077 with both approval lines FALSE; the gate normalises those lines)
  and `code/08_paine-holdout-fns.R` (c1119f5a69939a8c81b60ca20580b2e0). The gate `gate_hal9000_20261123_015214.rds` passed, including
  check (5): the stored declared runs use exactly the extra arm's hidden readings (1,974 of 1,974 plant-arm-levels).
- Runs (each about 60 core-s and 16 s wall; frozen md5s identical at start and end):

  | Level | Folder | Fits rds | Merge rds |
  |---|---|---|---|
  | H1 | `extraH1_andrewb_20261123_015525/` | 9bf0ffcb... | 767aade8... |
  | H2 | `extraH2_andrewb_20261123_015639/` | 1d72d9d0... | 00d71c75... |
  | H3 | `extraH3_andrewb_20261123_015702/` | c536ec55... | 1c3fb2a7... |

- natural_nls and log_nls are read from the stored declared runs, not refitted (md5s checked). V0 and M-J are A2's stored predictions.
- Error, width and direction are over readings where ALL FIVE tools worked. That set is a little smaller than the four-tool set, so
  the declared arms' medians here can differ slightly from the declared tables (for example H3 log_nls error 0.446 here, 0.482 in the
  declared table). The declared tables remain the reference for the declared arms.

## Failures
| Level | natural_improved plants failed (readings) | Reasons | natural_nls, for reference |
|---|---|---|---|
| H1 | 7 (7) | 4 chosen model was a fallback plain-nls fit; 3 no model could be fitted | 3 |
| H2 | 10 (20) | 7 fallback chosen; 3 no model | 5 |
| H3 | 69 (301) | 39 fallback chosen; 28 no model; 2 plants (1135FA1024, 1135FA1362) with predicted median <= 0 on some readings | 32 |

## Headline, summary A: failures as misses
| Level | Tool | Inclusion (90%) | abs log error, median | Width (mm^2), median | Median above reading |
|---|---|---|---|---|---|
| H1 (327 readings; 320 all five worked) | natural_nls | 36.1% | 0.089 | 571 | 74.1% |
| | **natural_improved** | **55.7%** | 0.088 | 1,040 | 72.5% |
| | log_nls | 92.7% | 0.152 | 8,280 | 71.2% |
| | M-J | 99.7% | 0.062 | 11,100 | 51.6% |
| | V0 | 96.6% | 1.17 | 60,800 | 96.6% |
| H2 (658; 634) | natural_nls | 34.3% | 0.098 | 492 | 76.0% |
| | **natural_improved** | **49.5%** | 0.097 | 885 | 72.9% |
| | log_nls | 92.7% | 0.191 | 8,740 | 60.4% |
| | M-J | 100.0% | 0.087 | 13,900 | 44.5% |
| | V0 | 93.0% | 1.48 | 124,000 | 99.7% |
| H3 (1,460; 1,147) | natural_nls | 43.4% | 0.258 | 1,470 | 38.5% |
| | **natural_improved** | **41.9%** | 0.235 | 1,840 | 41.3% |
| | log_nls | 88.1% | 0.446 | 7,780 | 30.9% |
| | M-J | 99.0% | 0.563 | 11,500 | 0.8% |
| | V0 | 85.2% | 1.25 | 135,000 | 98.3% |

- **Summary B** (only readings where all five worked): natural_improved inclusion is H1 56.9%, H2 50.8%, H3 53.4%; natural_nls on the
  same readings is 36.6%, 35.3% and 48.7%.
- **Without 1135FA1003:** no figure changes beyond the third significant digit, in either summary.
- **Paired natural_improved minus natural_nls** (same readings):

  | Level | Readings | Error difference | Width ratio |
  |---|---|---|---|
  | H1 | 320 | +0.002 | 1.71 (1.51 to 1.98) |
  | H2 | 634 | +0.000 | 1.59 (1.30 to 1.99) |
  | H3 | 1,157 | -0.000 | 1.17 (0.93 to 1.36) |

  The central predictions are essentially the same; the intervals are wider.
- **Fitted power delta among worked fits** (description; 1 would be noise proportional to size): median 0.50 (IQR 0.36 to 0.62) in H1,
  0.40 (0.24 to 0.55) in H2, 0.16 (0.04 to 0.37) in H3.
- **5% quantile <= 0:** natural_improved has 11 such readings in H1, 46 in H2 and 213 in H3.

## What this says about the hypothesis (plain statement; post-hoc)
- Adding a variance function that lets the noise grow with size WIDENED the natural arm's intervals (by about 1.6 to 1.7 times in H1
  and H2) and RAISED inclusion from about 35% to about 50% to 57%. In H3 it did neither by much, and it failed for more plants.
- Inclusion stays FAR BELOW the nominal 90%. log_nls (ln size, constant sd) reaches about 88% to 93% with intervals about 8 times
  wider.
- So the constant residual sd on the natural scale is at most PART of the cause of the plain arm's narrow intervals. It does not
  explain most of the shortfall.
- The fitted powers (about 0.2 to 0.5) are well below the value of 1 that proportional noise would give. A model fitted to readings
  up to the cut, where most readings are small, may not learn how large the noise is at the sizes being predicted.
- Further possible causes, ALL HYPOTHESES, NOT TESTED (no further variant is proposed here; any would need its own declaration):
  (i) parameter uncertainty understated by the linearised vcov when extrapolating beyond the last reading; (ii) the asymptotic
  curve shapes miss the late growth, so the error at the hidden day is a systematic offset that noise terms cannot cover (the
  natural arms' medians are ABOVE the reading for about 72% to 76% of H1 and H2 readings); (iii) day-to-day deviations of a plant
  are correlated, which an independent-noise model ignores.
- No second variant will be run after these results (declaration).

## Files (per extra folder)
`EXTRA_holdout_<L>_five_tools_headline.csv`, `..._summaryA_failures_as_misses.csv`, `..._summaryB_all_worked.csv`,
`EXTRA_holdout_<L>_failure_counts.csv`, `EXTRA_holdout_<L>_five_tools_long.csv`, `paine_improved_<L>_predictions.csv`,
`..._fits.csv` (with delta), `..._aic_tables.csv`, `..._fits.rds`, `EXTRA_paine_improved_<L>_merge.rds`, `EXTRA_LABEL.txt`, `log.txt`.

## ERRATUM, 23 Nov 2026 (Controller 02:10). The text above is left as it was; these corrections apply. Nothing was rerun.
1. "with intervals about 8 times wider" (log_nls against natural_improved) holds only for H1. From the headline table (ratio of median
   widths, summary A, readings where all five worked):
   - H1: 8,280 / 1,040 = about 8 times;
   - H2: 8,740 / 885 = about 10 times;
   - H3: 7,780 / 1,840 = about 4 times.
2. Most of natural_improved's failures (4 of 7 plants in H1, 7 of 10 in H2, 39 of 69 in H3) come from the frozen fallback rule, not
   from the variance function itself. When gnls fails for a model, `paine_fit_all` fits plain nls for that model instead, and that fit
   still competes in the AIC. Under this arm's declared failure rule, a chosen fallback (no variance function) counts as a failure. The
   arm's failure rate therefore reflects that rule as much as the variance function.
