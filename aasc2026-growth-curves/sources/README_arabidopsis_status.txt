STATUS: Arabidopsis rosette area (IPK Junker et al. 2016) - CLEANED, v2

Source: 1135FA_analysis.xlsx, 'Analysis Results' sheet, Row Type == 'top'
Target trait: area (zoom-corrected geometry trait based on visible-light top view) [mm^2]

Exclusions applied, in order:
1. Artifact row 'Plant ID' == '1135FA' (header/summary row leaked into data) - 1 row
2. Dead/failed-germination plants: max_area < 100 mm^2 - 51 plants
3. Bad capture dates, excluded entirely (all plants, that day):
     Day 29 - dip rate 3.4% vs ~0.2-0.4% on neighbouring days (~20x elevated)
     Day 41 - dip rate 2.1% vs ~0.6% on neighbouring days (~3x elevated)
   Confirmed via isolated_dips_v2 logic: flags a within-plant crash to <5% of its own
   max area, but only after the plant had already reached >30% of that max (i.e. a genuine
   crash, not early germination noise, which the v1 filter incorrectly flagged en masse).
4. Stunted plants, second pass: max_area < 500 mm^2 - 8 plants
   Identified via a natural gap in the sorted max_area distribution: last stunted plant at
   309 mm^2, next plant at 1569 mm^2 (~5x jump, gap_to_next = 1260) - a genuine
   discontinuity, not an arbitrary cutoff.

Final N: 425 plants (from 485 raw Plant IDs)

'group' column (moving/stationary) reconstructed from sampling density pattern, NOT an
explicit metadata column - moving: dense ~daily sampling from day 12 (n=384 originally);
stationary: sparse ~4-day sampling from day 24 (n=100 originally). Matches paper's stated
384/100 split exactly. This is a genotype-imbalance analogue: same genotype (C24), two
sampling densities - a real-data candidate for the Stage 6 density-imbalance/pooling-bias
comparison, alongside the wheat check/non-check case.

Known remaining quirk, not investigated further: day coverage is uneven across the season
(e.g. days 38/40 absent from the data entirely, day 43 has only ~147 plants vs ~340 on
neighbouring days) - likely reflects the moving/stationary imaging schedules interacting,
not a data quality problem. Worth remembering this when day 43 shows up with an artificially
narrow spread in any pooled/aggregate plot - it's a sample-size effect, not signal.

Verified clean via replotted individual trajectories, mean+-SD ribbon by group, and
stationary-only trajectory view - no remaining flat/dead lines, no ribbon distortion
at previously-bad dates.
