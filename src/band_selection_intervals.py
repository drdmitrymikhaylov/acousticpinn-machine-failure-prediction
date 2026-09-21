"""Section 3 read with its sample size attached.

exp2_detection.json and exp2_interferer.json hold, per SNR x strategy x
defect type, the detection rate at 5 % false alarms over 30 records.  The
page quotes the pooled rate over the three defect types (90 records) as a
bare number.  This script attaches:

  - the pooled count k / 90 and its Wilson 95 % interval,
  - a two-proportion z-test for kurtogram vs fixed band and oracle vs fixed
    band at each SNR (records are independent draws, so the pooled test is
    the honest one; the cells share the same machines across strategies,
    but per-record outcomes are not stored, so no paired test is possible),
  - the per-defect breakdown, to see whether one defect type carries a gap.

Writes results/band_selection_intervals.json.
"""
from __future__ import annotations

import json
import math
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
RES = ROOT / "results"
STRATEGIES = ("raw", "fixed", "kurtogram", "oracle")
DEFECTS = ("BPFO", "BPFI", "BALL")


def wilson(k, n, z=1.96):
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return centre - half, centre + half


def two_prop_z(k1, n1, k2, n2):
    p = (k1 + k2) / (n1 + n2)
    se = math.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
    if se == 0:
        return 0.0, 1.0
    z = (k1 / n1 - k2 / n2) / se
    pval = math.erfc(abs(z) / math.sqrt(2))  # two-sided
    return z, pval


def read(path):
    d = json.load(open(path))
    n = d["config"]["n_per_cell"]
    out = {"n_per_defect": n, "n_pooled": n * len(DEFECTS), "snr": {}}
    for snr in d["config"]["snrs"]:
        row = {}
        ks = {}
        for st in STRATEGIES:
            per = {df: d["cells"][f"{snr}|{st}|{df}"]["recall_at_5pct_fa"] for df in DEFECTS}
            k = int(round(sum(per.values()) * n))
            ks[st] = k
            lo, hi = wilson(k, n * len(DEFECTS))
            row[st] = {"per_defect": per, "k": k, "rate": k / (n * len(DEFECTS)),
                       "wilson95": [lo, hi]}
        N = n * len(DEFECTS)
        for a, b in (("kurtogram", "fixed"), ("oracle", "fixed"), ("oracle", "kurtogram")):
            z, p = two_prop_z(ks[a], N, ks[b], N)
            row[f"{a}_minus_{b}"] = {"diff_records": ks[a] - ks[b], "diff_rate": (ks[a] - ks[b]) / N,
                                     "z": z, "p_two_sided": p}
        out["snr"][str(snr)] = row
    out["kurtogram_centre_error_hz"] = d["kurtogram_centre_error_hz"]
    return out


def main():
    out = {"white": read(RES / "exp2_detection.json"),
           "interferer": read(RES / "exp2_interferer.json")}
    json.dump(out, open(RES / "band_selection_intervals.json", "w"), indent=1)
    for name, r in out.items():
        print("==", name)
        for snr, row in r["snr"].items():
            f, k, o = row["fixed"], row["kurtogram"], row["oracle"]
            kf, of = row["kurtogram_minus_fixed"], row["oracle_minus_fixed"]
            print(f"{snr:>4} dB  fixed {f['k']:2d}/90 [{f['wilson95'][0]:.2f},{f['wilson95'][1]:.2f}]  "
                  f"kurt {k['k']:2d}/90 [{k['wilson95'][0]:.2f},{k['wilson95'][1]:.2f}]  "
                  f"oracle {o['k']:2d}/90  kurt-fixed {kf['diff_records']:+3d} z={kf['z']:+.2f} p={kf['p_two_sided']:.3f}  "
                  f"oracle-fixed {of['diff_records']:+3d} p={of['p_two_sided']:.1e}")


if __name__ == "__main__":
    main()
