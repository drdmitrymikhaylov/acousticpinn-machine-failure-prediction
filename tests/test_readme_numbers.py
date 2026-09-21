"""Pin the numbers quoted on the project page to results/*.json.

Run: python -m pytest tests/test_readme_numbers.py   (stdlib only)
"""
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
DEFECTS = ("BPFO", "BPFI", "BALL")


def load(name):
    with open(os.path.join(RES, name)) as f:
        return json.load(f)


def pooled(d, snr, strategy):
    return sum(d["cells"][f"{snr}|{strategy}|{df}"]["recall_at_5pct_fa"] for df in DEFECTS) / len(DEFECTS)


def close(a, b, tol):
    assert abs(a - b) <= tol, (a, b, tol)


# ---------------------------------------------------------------- §1

def test_section1_raw_vs_envelope_table():
    d = load("exp2_detection.json")
    assert d["config"]["n_per_cell"] == 30
    for snr, raw, fixed in ((-6, 0.56, 1.00), (-12, 0.24, 1.00), (-15, 0.07, 0.99)):
        close(pooled(d, snr, "raw"), raw, 0.005)
        close(pooled(d, snr, "fixed"), fixed, 0.005)


# ---------------------------------------------------------------- §3

def test_section3_minus18_rates_and_kurtogram_error():
    d = load("exp2_detection.json")
    close(pooled(d, -18, "oracle"), 0.97, 0.005)
    close(pooled(d, -18, "fixed"), 0.47, 0.005)
    close(pooled(d, -18, "kurtogram"), 0.58, 0.005)
    close(d["kurtogram_centre_error_hz"]["median"], 340, 0.5)
    close(d["kurtogram_centre_error_hz"]["p90"], 2600, 20)
    i = load("exp2_interferer.json")
    close(i["kurtogram_centre_error_hz"]["p90"], 6600, 50)


def test_section3_counts_intervals_and_tests_match_cells():
    b = load("band_selection_intervals.json")
    for name, fn in (("white", "exp2_detection.json"), ("interferer", "exp2_interferer.json")):
        d = load(fn)
        r = b[name]
        assert r["n_pooled"] == 90
        for snr in d["config"]["snrs"]:
            row = r["snr"][str(snr)]
            for st in ("raw", "fixed", "kurtogram", "oracle"):
                k = round(pooled(d, snr, st) * 90)
                assert row[st]["k"] == k, (name, snr, st)
                lo, hi = row[st]["wilson95"]
                assert lo - 1e-9 <= k / 90 <= hi + 1e-9
            assert row["kurtogram_minus_fixed"]["diff_records"] == row["kurtogram"]["k"] - row["fixed"]["k"]
    w = b["white"]["snr"]
    assert (w["-18"]["fixed"]["k"], w["-18"]["kurtogram"]["k"], w["-18"]["oracle"]["k"]) == (42, 52, 87)
    assert w["-18"]["kurtogram_minus_fixed"]["diff_records"] == 10
    close(w["-18"]["kurtogram_minus_fixed"]["p_two_sided"], 0.14, 0.005)
    assert w["-18"]["oracle_minus_fixed"]["p_two_sided"] < 1e-12
    # white background: kurtogram at or below fixed down to -15 dB, ahead only at -18
    for snr in ("0", "-3", "-6", "-9", "-12", "-15"):
        assert w[snr]["kurtogram"]["k"] <= w[snr]["fixed"]["k"]
    assert w["-18"]["kurtogram"]["k"] > w["-18"]["fixed"]["k"]
    # per-defect split at -18 dB white
    d = load("exp2_detection.json")
    fixed = [round(d["cells"][f"-18|fixed|{df}"]["recall_at_5pct_fa"] * 30) for df in DEFECTS]
    kurt = [round(d["cells"][f"-18|kurtogram|{df}"]["recall_at_5pct_fa"] * 30) for df in DEFECTS]
    assert fixed == [12, 17, 13] and kurt == [15, 18, 19]
    i = b["interferer"]["snr"]
    for snr in ("0", "-3", "-6", "-9", "-12"):
        t = i[snr]["kurtogram_minus_fixed"]
        assert -10 <= t["diff_records"] <= -8 and t["p_two_sided"] < 0.05, snr
    for snr in ("-15", "-18"):
        assert i[snr]["kurtogram_minus_fixed"]["diff_records"] == 5
        assert i[snr]["kurtogram_minus_fixed"]["p_two_sided"] > 0.3
    assert (i["-18"]["fixed"]["k"], i["-18"]["kurtogram"]["k"], i["-18"]["oracle"]["k"]) == (37, 42, 84)
    assert i["-15"]["oracle_minus_fixed"]["diff_records"] == 19


# ---------------------------------------------------------------- §5

def test_section5_alarm_policy_table():
    a = load("exp4_alarm.json")["summary"]
    cfg = a["config"]
    assert cfg["n_machines"] == 24 and cfg["baseline_n"] == 144 and cfg["growth_n"] == 288
    assert cfg["interval_min"] == 10.0
    pol = a["policies"]
    table = {"0.95|1": (315, 46.7), "0.99|1": (112.5, 44.3), "0.999|1": (82.5, 41.5),
             "0.95|2": (22.5, 34.0), "0.99|2": (2.5, 30.8), "0.99|3": (0.0, 27.5)}
    for key, (fa, warn) in table.items():
        close(pol[key]["false_alarms_per_month"], fa, 0.6)
        close(pol[key]["warning_before_end_h"], warn, 0.06)
    # "none observed" in 1728 held-out opportunities: 24 machines x 72 samples
    assert cfg["n_machines"] * cfg["baseline_n"] // 2 == 1728
    close(pol["0.99|3"]["median_snr_at_alarm_db"], -18, 0.5)


# ---------------------------------------------------------------- §6

def test_section6_resonance_error_table():
    e = load("exp5_resonance_pinn.json")
    assert e["config"]["n_machines"] == 12 and e["config"]["n_win"] == 16
    want = {"-6": (130, 12, 17, 56), "-9": (253, 15, 25, 492),
            "-12": (189, 178, 77, 3092), "-15": (343, 847, 98, 2836)}
    for snr, (kurt, ds, joint, pinn) in want.items():
        r = e["resonance_error"][snr]
        close(r["kurtogram"]["median_abs_err_hz"], kurt, 0.5)
        close(r["damped_sine"]["median_abs_err_hz"], ds, 0.5)
        close(r["joint"]["median_abs_err_hz"], joint, 0.5)
        close(r["pinn"]["median_abs_err_hz"], pinn, 0.5)
    c = e["cells"]["-15"]
    close(c["damped_sine"]["auc"], 0.90, 0.005)
    close(c["joint"]["auc"], 0.91, 0.005)
    close(c["oracle"]["auc"], 1.00, 0.005)
    close(e["resonance_error"]["-15"]["zeta_rel_err"]["damped_sine"], 5.5, 0.1)
