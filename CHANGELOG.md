# Changelog

## 2026-10-06

- Exp 5 (alarm policy) re-read as counts. The 24 machines were re-run with
  the same seeds (the published table reproduces exactly) and every series
  is now in `results/exp4_alarm_series.json`; the result file had kept one.
  315 / 112.5 / 82.5 / 22.5 / 2.5 false alarms a month are 126 / 45 / 33 / 9 /
  1 events in 1 728 held-out samples. `results/alarm_counts.json` adds
  intervals from resampling machines, the three policies the table left out,
  and the worst machine's warning time (23.5 h against the 27.5 h median).
- "None observed" sized by cutting each quiet day again 2 000 times: three in
  a row expects 1.6 a month at the 95th percentile, 0.13 at the 99th, 0.07 at
  the 99.9th. Correction: finding 6 quoted 29 h of warning, which is the
  95th-percentile policy and was not in the table; it now quotes 27.5 h.
- A 99th and a 99.9th percentile of 72 samples sit between the same two
  samples and are exceeded 1.4-2.7 % of the time for any distribution.
- Sensitivity: with lag-1 correlation 0.6 in the baseline, three in a row at
  the 99th percentile gives 13 false alarms a month instead of 0.15.
- `tests/test_alarm_counts.py` (5 tests). The bearing maker's name removed
  from identifiers and comments in `src/` and `tests/` (geometry unchanged);
  three stray `._*` files removed.

## 2026-09-21

- §3 band selection re-read with n attached: every rate is 90 records
  (3 defect types × 30). Correction: on white background the kurtogram is
  *ahead* of the fixed band at −18 dB (52 vs 42 of 90, p = 0.14, inside the
  noise), not "at or below everywhere"; with an interferer it is
  significantly behind at every SNR from 0 to −12 dB (8–10 records, p =
  0.003–0.03) and ahead by 5 records (p ≈ 0.4) at −15/−18. Oracle − fixed
  at −18 dB: +45 records, p ≈ 10⁻¹³. Rule 4 reworded. New
  `results/band_selection_intervals.json`, `src/band_selection_intervals.py`.
- `tests/test_readme_numbers.py` (5 tests) pins §1, §3, §5, §6 tables to
  `results/`.
