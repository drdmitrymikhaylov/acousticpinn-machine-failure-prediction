"""Exp 5 as counts: recompute the alarm table from the 24 published series and
pin the second reading on the page to results/alarm_counts.json.

Run: python -m pytest tests/test_alarm_counts.py   (numpy only)
"""
import json
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
ORDER = ["0.95|1", "0.99|1", "0.999|1", "0.95|2", "0.99|2", "0.999|2", "0.95|3", "0.99|3", "0.999|3"]
LABEL = {0.95: "95th pct", 0.99: "99th pct", 0.999: "99.9th pct"}


def load(name):
    with open(os.path.join(RES, name)) as f:
        return json.load(f)


def readme():
    with open(os.path.join(ROOT, "README.md"), encoding="utf-8") as f:
        return " ".join(f.read().split())


def close(a, b, tol=1e-9):
    assert abs(a - b) <= tol, (a, b, tol)


def num(x):
    """315.0 -> '315', 112.5 -> '112.5', 0.13 -> '0.13'."""
    if x >= 10 or abs(x - round(x)) < 1e-9:
        s = "%.1f" % x
        return s[:-2] if s.endswith(".0") else s
    return "%.1f" % x if x >= 1 else "%.2f" % x


def count_row(p):
    lo, hi = p["ci95_machines_resampled"]
    rate = ("none observed" if p["false_alarm_events"] == 0
            else "%s [%s, %s]" % (num(p["false_alarms_per_month"]), num(lo), num(hi)))
    return "| %s | %d | %d | %d | %s | %s | %.1f h | %.1f h |" % (
        LABEL[p["quantile"]], p["consecutive"], p["false_alarm_events"],
        p["machines_with_a_false_alarm"], rate, num(p["reshuffled_split"]["mean_per_month"]),
        p["warning_h"]["median"], p["warning_h"]["min"])


def memory_row(rho, row):
    single = row["0.99|1"]["false_alarms_per_month"]
    triple = row["0.99|3"]["false_alarms_per_month"]
    ratio = single / triple
    return "| %s | %s | %s | %s | × %s |" % (
        rho, num(single), num(triple), num(row["0.95|3"]["false_alarms_per_month"]),
        ("%.0f" % ratio) if ratio >= 10 else ("%.1f" % ratio))


def hits(exceed, m):
    n = exceed.shape[-1]
    h = exceed[..., : n - m + 1].copy()
    for j in range(1, m):
        h &= exceed[..., j: n - m + 1 + j]
    return h


def first_alarm(series, thr, m):
    run = 0
    for i, v in enumerate(series):
        run = run + 1 if v > thr else 0
        if run >= m:
            return i
    return None


def test_counts_recomputed_from_the_24_series():
    s = load("exp4_alarm_series.json")
    c = load("alarm_counts.json")
    pub = load("exp4_alarm.json")["summary"]["policies"]
    base = np.array([m["baseline"] for m in s["machines"]])
    grow = np.array([m["growth"] for m in s["machines"]])
    assert base.shape == (24, 144) and grow.shape == (24, 288)
    # machine 0 is the example series the original result file kept
    ex = load("exp4_alarm.json")["example"]
    assert np.allclose(base[0], ex["baseline"], atol=1e-6) and np.allclose(grow[0], ex["growth"], atol=1e-6)
    for key in ORDER:
        p = c["policies"][key]
        q, m = p["quantile"], p["consecutive"]
        thr = np.quantile(base[:, :72], q, axis=1)
        per_machine = hits(base[:, 72:] > thr[:, None], m).sum(axis=1)
        assert [int(v) for v in per_machine] == p["events_per_machine"]
        assert int(per_machine.sum()) == p["false_alarm_events"]
        assert int((per_machine > 0).sum()) == p["machines_with_a_false_alarm"]
        close(p["false_alarms_per_month"], per_machine.sum() / 1728 * 4320)
        close(p["false_alarms_per_month"], pub[key]["false_alarms_per_month"])
        warn = 48.0 - np.array([first_alarm(grow[k], thr[k], m) for k in range(24)]) / 6.0
        close(p["warning_h"]["median"], float(np.median(warn)))
        close(p["warning_h"]["median"], pub[key]["warning_before_end_h"])
        close(p["warning_h"]["min"], float(warn.min()))


def test_count_table_on_the_page():
    c = load("alarm_counts.json")["policies"]
    text = readme()
    for key in ORDER:
        assert count_row(c[key]) in text, count_row(c[key])
    assert [c[k]["false_alarm_events"] for k in ORDER] == [126, 45, 33, 9, 1, 1, 0, 0, 0]
    assert "one event" in text and "126 events" in text


def test_none_observed_covers_three_different_rates():
    c = load("alarm_counts.json")["policies"]
    text = readme()
    r95, r99, r999 = (c[k]["reshuffled_split"] for k in ("0.95|3", "0.99|3", "0.999|3"))
    assert num(r95["mean_per_month"]) == "1.6" and num(r99["mean_per_month"]) == "0.13"
    assert num(r999["mean_per_month"]) == "0.07"
    assert round(100 * r95["share_of_splits_with_no_event"]) == 57
    assert round(100 * r99["share_of_splits_with_no_event"]) == 95
    for phrase in ("expects 1.6 false alarms a month", "57 % of the cuts",
                   "expects 0.13 a month", "95 % of the cuts",
                   "one false alarm in about eight months"):
        assert phrase in text, phrase
    assert 7.0 < 1 / r99["mean_per_month"] < 8.5
    # finding 6 now quotes the row that is in the table
    close(c["0.99|3"]["warning_h"]["median"], 27.5)
    close(c["0.95|3"]["warning_h"]["median"], 29.1667, 1e-3)
    assert "keeping 27.5 h of the 48 h warning" in text
    assert "leaving 27.5 of the 48 hours" in text
    close(c["0.99|3"]["warning_h"]["min"], 23.5)
    assert "still has 23.5 h against a median of 27.5 h" in text


def test_a_percentile_of_72_samples():
    q = load("alarm_counts.json")["quantile_of_72_samples"]
    assert q["0.99"]["between_order_statistics"] == [71, 72] == q["0.999"]["between_order_statistics"]
    assert round(100 * q["0.99"]["interpolation_weight_on_upper"]) == 29
    assert round(100 * q["0.999"]["interpolation_weight_on_upper"]) == 93
    lo, hi = q["0.999"]["expected_rate_bounds"]
    close(lo, 1 / 73); close(hi, 2 / 73)
    assert (round(lo * 4320), round(hi * 4320)) == (59, 118)
    assert round(100 * q["0.99"]["observed_rate"], 1) == 2.6 and round(100 * q["0.999"]["observed_rate"], 1) == 1.9
    assert round(100 * q["0.99"]["reshuffled_mean_rate"], 1) == 2.2
    assert round(100 * q["0.999"]["reshuffled_mean_rate"], 1) == 1.4
    text = readme()
    for phrase in ("29 % and 93 % of the way", "between 1.4 % and 2.7 % of the time",
                   "59 to 118 false alarms a month", "2.6 % and 1.9 %", "2.2 % and 1.4 %"):
        assert phrase in text, phrase


def test_memory_table_and_its_anchor():
    c = load("alarm_counts.json")
    text = readme()
    for rho in ("0.0", "0.3", "0.6", "0.9"):
        assert memory_row(rho, c["memory"][rho]) in text, memory_row(rho, c["memory"][rho])
    # at zero correlation the synthetic fleet must agree with the reshuffled bearing fleet
    for key in ORDER:
        a = c["memory"]["0.0"][key]["false_alarms_per_month"]
        b = c["policies"][key]["reshuffled_split"]["mean_per_month"]
        assert abs(a - b) <= max(0.06 * b, 0.1), (key, a, b)
    # the simulated baseline has no memory
    s = load("exp4_alarm_series.json")
    base = np.array([m["baseline"] for m in s["machines"]])
    d = base - base.mean(axis=1, keepdims=True)
    lag1 = (d[:, 1:] * d[:, :-1]).sum(axis=1) / (d ** 2).sum(axis=1)
    close(c["baseline_memory"]["lag1_autocorrelation_mean"], float(lag1.mean()), 1e-9)
    assert abs(lag1.mean()) < 0.02 and np.abs(lag1).max() < 0.2
    assert "lag-1 autocorrelation of the 24 baselines averages −0.007" in text
    m6 = c["memory"]["0.6"]["0.99|3"]["false_alarms_per_month"]
    assert round(m6) == 13 and "gives 13 a month" in text
