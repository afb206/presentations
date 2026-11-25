STATUS as of [today's date]:
- height_data_clean: FPWW028 excluded (small/problematic 2019 lot), lodging NOT yet filtered
- 2016-06-27 11:00:00 (and possibly others) look like bad capture dates, not real lodging - unresolved
- lodging_candidates_v4 is the latest attempt but likely still contaminated by bad-date artefacts
- Next step: fix date_level_check using corrected diff_next_rel logic, re-derive bad_dates,
  then re-run lodging_candidates with clean dates excluded
- lodged_by_reference (QMER-comparison method) is comparatively more trustworthy so far
