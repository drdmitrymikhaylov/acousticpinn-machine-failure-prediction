# Changelog

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
