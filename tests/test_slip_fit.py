"""Pin the Exp 2 subsection "What 'verified to 3 %' rested on" (added
2026-11-18) to results/slip_fit_check.json and results/exp1_slip.json.

The stored fit is recomputed here from the stored amplitudes with the same
rule (log a against k^2, a > 0.15), so the table's k_max column, the
harmonic-range table and the per-harmonic reading are all re-derived; the
36-record and long-record numbers are read from the JSON.  numpy + pytest.
"""
import json
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
SLIPS = ["0.0025", "0.005", "0.01", "0.02", "0.04"]
LABEL = {"0.0025": "0.25%", "0.005": "0.5%", "0.01": "1%", "0.02": "2%", "0.04": "4%"}


def load(name):
    with open(os.path.join(RES, name)) as f:
        return json.load(f)


def readme():
    with open(os.path.join(ROOT, "README.md"), encoding="utf-8") as f:
        return f.read()


def flat(t):
    return " ".join(t.split())


C = load("slip_fit_check.json")
E = load("exp1_slip.json")
TEXT = readme()
FLAT = flat(TEXT)
K = np.arange(1, 31)


def sigma_fit(a, k, cut=0.15):
    use = a > cut
    slope = np.polyfit(k[use] ** 2, np.log(a[use]), 1)[0]
    return float(np.sqrt(max(-slope, 1e-12) / (2 * np.pi ** 2)))


def row1(key):
    r = C["rows"][key]
    return "| %s | %.1f | %.1f | %+.1f %% | %.0f %% | %.3f | %.3f [%.3f, %.3f] | %.1f %% | %d / 36 |" % (
        LABEL[key], r["k_max_theory"], r["k_max_measured"], 100 * (r["k_max_ratio"] - 1),
        100 * r["decay_covered"], r["observed_a_at_k30"], r["k_max_ratio_36_seeds"], *r["k_max_ratio_36_seeds_ci95"],
        100 * r["per_seed_sigma_sd"] / r["slip"], r["per_seed_within_3pct"])


def row2(key):
    r = C["rows"][key]
    return "| %s | %s |" % (LABEL[key], " | ".join(
        "%+.1f %%" % (100 * (r["slip"] / r["sigma_by_kmax"][k] - 1)) for k in ("8", "12", "16", "20", "30")))


def test_stored_fit_recomputed_and_reproduced():
    """The table's k_max column from the stored amplitudes; the JSON agrees; the
    re-simulation reproduced the stored means."""
    for key in SLIPS:
        s = float(key)
        a = np.array(E["iid"][key]["amplitude_shape_corrected"])
        sig = sigma_fit(a, K)
        assert abs(sig - E["iid_fit"][key]["sigma_fitted"]) < 1e-12
        r = C["rows"][key]
        assert r["refit_matches_stored"] and abs(r["sigma_refit"] - sig) < 1e-12
        assert abs(r["k_max_measured"] - 0.2251 / sig) < 1e-9
        assert abs(r["k_max_theory"] - 0.2251 / s) < 1e-12
        assert abs(r["k_max_ratio"] - s / sig) < 1e-12
        assert abs(r["decay_covered"] - min(1.0, 2 * np.pi ** 2 * min(30, 0.2251 / s) ** 2 * s * s)) < 1e-12
        assert r["observed_a_at_k30"] == a[29]
        assert r["six_seed_mean_reproduced_max_abs"] < 1e-12
    assert C["reference_reproduced_max_abs"] < 1e-12
    assert "re-simulation to 2 × 10⁻¹⁶" in FLAT and C["reference_reproduced_max_abs"] < 3e-16
    # the deviations finding 2 now quotes
    dev = {key: round(100 * (C["rows"][key]["k_max_ratio"] - 1), 1) for key in SLIPS}
    assert dev == {"0.0025": -6.4, "0.005": -1.8, "0.01": 0.8, "0.02": 1.0, "0.04": 2.5}, dev
    assert "84 against 90 at 0.25 % slip is −6.4 %, and 5.8 against 5.6 at 4 % is +2.5 %" in FLAT
    assert "within 1–2.5% for slips of 0.5–4%; 6% low at 0.25% slip" in FLAT
    assert "verified to 3% over a 16× range of slip |" not in TEXT   # the old cell is gone


def test_main_table_rows_and_coverage():
    for key in SLIPS:
        assert row1(key) in TEXT, row1(key)
    r = C["rows"]["0.0025"]
    assert round(r["predicted_a_at_k30"], 3) == 0.895 and round(r["observed_a_at_k30"], 3) == 0.897
    assert "fallen to 0.895 of the first, and the stored mean reads 0.897" in FLAT
    assert round(100 * r["decay_covered"]) == 11 and round(r["k_max_theory"]) == 90
    assert r["k_max_ratio_36_seeds_ci95"][1] < 1.0 and r["k_max_ratio_36_seeds_ci95"][1] < 0.97
    assert C["n_seed_many"] == 36 and C["n_seed_stored"] == 6
    # 90 harmonics of f_defect against a 4 kHz band
    assert round(90 * C["f_defect"] / 1000, 1) == 9.7 and C["band_hz"] == [6000.0, 10000.0]
    assert int(4000 // C["f_defect"]) == 37
    assert "past the ~37 that a 4 kHz demodulation band can supply" in FLAT
    # single-record spread 2.4-4.2 %; 36-record spreads
    sds = [100 * C["rows"][k]["per_seed_sigma_sd"] / C["rows"][k]["slip"] for k in SLIPS]
    assert (round(min(sds), 1), round(max(sds), 1)) == (2.4, 4.2)
    assert "recovers the slip to 2.4–4.2 % (sd)" in FLAT
    # 1-2.5 % where measurable
    devs = [100 * abs(C["rows"][k]["k_max_ratio"] - 1) for k in SLIPS[1:]]
    assert 0.5 < min(devs) < 1.0 and round(max(devs), 1) == 2.5
    assert "confirmed to 1–2.5 % where it can be measured, 0.5–4 % slip" in FLAT


def test_long_record():
    lr = C["long_record"]["0.0025"]
    assert lr["duration_s"] == 32.0 and lr["n_seed"] == 12 and C["duration_s"] == 8.0
    assert "%.3f [%.3f, %.3f]" % (lr["k_max_ratio"], *lr["k_max_ratio_ci95"]) == "0.926 [0.921, 0.931]"
    assert lr["k_max_ratio"] < C["rows"]["0.0025"]["k_max_ratio_36_seeds"]
    assert "moves it to **0.926 [0.921, 0.931]**, further from 1" in FLAT
    sd8 = 100 * C["rows"]["0.0025"]["per_seed_sigma_sd"] / 0.0025
    sd32 = 100 * lr["per_seed_sigma_sd"] / 0.0025
    assert (round(sd8, 1), round(sd32, 1)) == (2.4, 1.1)
    assert "spread halves (2.4 % → 1.1 %)" in FLAT


def test_harmonic_range_table_recomputed():
    for key in SLIPS:
        r = C["rows"][key]
        a = np.array(E["iid"][key]["amplitude_shape_corrected"])
        for Kmax in (8, 12, 16, 20, 30):
            assert abs(sigma_fit(a[:Kmax], K[:Kmax]) - r["sigma_by_kmax"][str(Kmax)]) < 1e-12
        assert row2(key) in TEXT, row2(key)
    r = C["rows"]["0.0025"]
    vals = [round(100 * (0.0025 / r["sigma_by_kmax"][k] - 1), 1) for k in ("8", "16", "20")]
    assert vals == [-3.6, 3.8, -26.7], vals
    assert "runs from −3.6 % to +3.8 % to −26.7 % as harmonics are added" in FLAT


def test_per_harmonic_reading_and_leave_one_out():
    a = np.array(E["iid"]["0.0025"]["amplitude_shape_corrected"])
    ref = np.array(E["reference_amplitude"])
    assert [round(v, 3) for v in ref[18:20]] == [0.078, 0.032]
    assert [round(v, 3) for v in a[17:20]] == [0.999, 0.774, 0.972]
    assert "reference falls from 0.078 at k = 19 to 0.032 at k = 20" in FLAT
    assert "at k = 19 is 0.774 against 0.999 at k = 18 and 0.972 at k = 20" in FLAT
    with np.errstate(divide="ignore", invalid="ignore"):
        sk = np.sqrt(np.clip(-np.log(a), 0, None) / (2 * np.pi ** 2 * K ** 2)) / 0.0025
    stored = np.array(C["rows"]["0.0025"]["sigma_per_harmonic"]) / 0.0025
    assert np.allclose(sk, stored, atol=1e-12)
    assert (round(sk[2:17].min(), 2), round(sk[2:17].max(), 2)) == (0.88, 1.29)
    assert round(sk[18], 1) == 2.4 and round(sk[17], 2) == 0.19
    assert "harmonics 3 to 17 put the slip at 0.88 to 1.29 of the truth, k = 19 at 2.4 and k = 18 at 0.19" in FLAT
    # leave-one-out: < 0.5 % from 1 % slip up, recomputed
    for key in ("0.01", "0.02", "0.04"):
        a = np.array(E["iid"][key]["amplitude_shape_corrected"])
        sig = sigma_fit(a, K)
        use = a > 0.15
        worst = 0.0
        for j in np.flatnonzero(use):
            m = use.copy()
            m[j] = False
            sl = np.polyfit(K[m] ** 2, np.log(a[m]), 1)[0]
            worst = max(worst, abs(np.sqrt(-sl / (2 * np.pi ** 2)) / sig - 1))
        assert worst < 0.005 and abs(worst - C["rows"][key]["leave_one_harmonic_out_max_abs"]) < 1e-12
    assert "no single harmonic moves the fit by more than 0.5 % (leave-one-out)" in FLAT


def test_two_sentences_checked():
    # uncorrected fit: all 30 harmonics 0.015-0.019; with the 0.15 cut 0.032 ... 0.070
    allk, cut = [], []
    for key in SLIPS:
        raw = np.array(E["iid"][key]["amplitude"])
        allk.append(sigma_fit(raw, K, cut=0.0))
        cut.append(sigma_fit(raw, K))
        assert abs(cut[-1] - C["rows"][key]["sigma_uncorrected"]) < 1e-12
    assert (round(min(allk), 3), round(max(allk), 3)) == (0.015, 0.019)
    assert ["%.3f" % v for v in cut] == ["0.032", "0.032", "0.037", "0.049", "0.070"]
    assert "it returns 0.015–0.019 at every slip, as stated" in FLAT
    assert "returns **0.032, 0.032, 0.037, 0.049, 0.070** at 0.25 % to 4 %" in FLAT
    # the pedestal: predicted vs observed a(30) at 2 % and 4 %, and the 0.05 cut
    for key, pred, obs in (("0.02", "0.001", "0.101"), ("0.04", "0.000", "0.105")):
        r = C["rows"][key]
        assert "%.3f" % r["predicted_a_at_k30"] == pred and "%.3f" % r["observed_a_at_k30"] == obs
    off = [round(100 * (C["rows"][k]["slip"] / C["rows"][k]["sigma_by_cut"]["0.05"] - 1)) for k in ("0.02", "0.04")]
    assert off == [67, 338], off
    assert "puts k_max off by **+67 % and +338 %**" in FLAT
    assert "`results/slip_fit_check.json`" in TEXT
