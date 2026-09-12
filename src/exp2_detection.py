"""Does adaptive band selection earn its keep?

A deployed diagnostic module has to choose a demodulation band before it can
look for a defect frequency.  Four strategies are compared over a population
of machines whose structural resonance is *not* known in advance:

  raw        no demodulation at all -- look for the defect frequency directly
             in the spectrum of the vibration signal.  This is the naive
             configuration and the one to beat.
  fixed      one wide band (2-20 kHz) configured once for the whole fleet.
  kurtogram  the band is chosen per record by maximising the kurtosis of the
             band-limited envelope.
  oracle     the band is centred on the true resonance of that record.  Not
             available in the field; it bounds what band selection can buy.

Every record draws its own shaft speed, resonance frequency, damping and slip,
so a strategy cannot win by being tuned to one machine.  Healthy records use
the identical background, so the only difference between the classes is the
defect itself.

Reported: detection rate at a threshold fixed to give 5% false alarms on the
healthy population, and ROC AUC, as a function of signal-to-noise ratio.
"""

import json
import sys
from multiprocessing import Pool

sys.path.insert(0, 'src')

import numpy as np

from bearing_sim import Resonance, band_limited_noise, simulate, simulate_healthy
from envelope import envelope_spectrum, harmonic_statistic, kurtogram
from kinematics import SKF_6205

FS = 48000.0
DURATION = 1.0
N_HARM = 6
HALFWIDTH_HZ = 2.0
SNRS = [0, -3, -6, -9, -12, -15, -18]
N_PER_CELL = 30
FAULTS = ["BPFO", "BPFI", "BALL"]
STRATEGIES = ["raw", "fixed", "kurtogram", "oracle"]
FIXED_BAND = (2000.0, 20000.0)


def draw_machine(rng):
    """A machine: shaft speed, structural resonance, damping, slip."""
    fr = rng.uniform(20.0, 60.0)
    return dict(
        fr=fr,
        f_n=rng.uniform(4000.0, 14000.0),
        zeta=rng.uniform(0.010, 0.050),
        slip=rng.uniform(0.005, 0.020),
        slip_model="walk" if rng.random() < 0.5 else "iid",
    )


def defect_frequency(fault, fr):
    o = SKF_6205.orders()
    if fault == "BPFO":
        return o["BPFO"] * fr, 0.0, None
    if fault == "BPFI":
        return o["BPFI"] * fr, 0.8, fr
    if fault == "BALL":
        return o["BALL_IMPACT"] * fr, 0.6, o["FTF"] * fr
    raise ValueError(fault)


def statistic(x, f_defect, strategy, f_n):
    """Detection statistic in dB under one band-selection strategy."""
    bw_min = N_HARM * f_defect
    if strategy == "raw":
        # no demodulation: the defect line is sought in the plain spectrum
        w = np.hanning(len(x))
        amp = np.abs(np.fft.rfft(x * w)) / (w.sum() / 2.0)
        freqs = np.fft.rfftfreq(len(x), 1.0 / FS)
        s, _, _ = harmonic_statistic(freqs, amp, f_defect, N_HARM, HALFWIDTH_HZ)
        return s, np.nan
    if strategy == "fixed":
        band = FIXED_BAND
    elif strategy == "oracle":
        half = max(0.75 * bw_min, 1500.0)
        band = (max(f_n - half, 200.0), min(f_n + half, FS / 2 - 1))
    elif strategy == "kurtogram":
        _, best = kurtogram(x, FS, n_levels=6, f_min=500.0, bw_min=bw_min)
        band = (best[1] - best[2] / 2, best[1] + best[2] / 2)
    else:
        raise ValueError(strategy)
    f, E = envelope_spectrum(x, FS, band=band)
    s, _, _ = harmonic_statistic(f, E, f_defect, N_HARM, HALFWIDTH_HZ)
    return s, 0.5 * (band[0] + band[1])


INTERFERER = None          # set by main() to (band, power_db) or left None


def add_interferer(x, rng):
    if INTERFERER is None:
        return x
    band, power_db = INTERFERER
    p = 10 ** (power_db / 10.0)          # relative to the unit impulse power
    return x + band_limited_noise(len(x), FS, band, p, rng)


def one_record(args):
    snr, fault, idx = args
    rng = np.random.default_rng(hash((snr, fault, idx)) % (2 ** 32))
    m = draw_machine(rng)
    seed = int(rng.integers(1 << 30))
    res = Resonance(f_n=m["f_n"], zeta=m["zeta"])

    if fault == "HEALTHY":
        # a healthy record still has to be tested against *some* hypothesis:
        # the defect frequency the operator is monitoring for.
        f_defect, _, _ = defect_frequency("BPFO", m["fr"])
        x, _ = simulate_healthy(m["fr"], FS, DURATION, snr, seed=seed)
        x = add_interferer(x, np.random.default_rng(seed + 7))
    else:
        f_defect, mod_depth, mod_hz = defect_frequency(fault, m["fr"])
        x, _ = simulate(f_defect, m["fr"], fs=FS, duration=DURATION,
                        snr_db=snr, slip=m["slip"], slip_model=m["slip_model"],
                        resonance=res, load_modulation=mod_depth,
                        modulation_hz=mod_hz, seed=seed)
        x = add_interferer(x, np.random.default_rng(seed + 7))

    out = {"snr": snr, "fault": fault, "f_n": m["f_n"], "fr": m["fr"]}
    for st in STRATEGIES:
        s, fc = statistic(x, f_defect, st, m["f_n"])
        out[st] = s
        if st == "kurtogram":
            out["kurtogram_fc"] = fc
    return out


def auc(pos, neg):
    pos, neg = np.asarray(pos), np.asarray(neg)
    allv = np.concatenate([pos, neg])
    r = np.argsort(np.argsort(allv)) + 1.0
    rp = r[: len(pos)].sum()
    return float((rp - len(pos) * (len(pos) + 1) / 2.0) / (len(pos) * len(neg)))


def main(interferer=None, out_path="results/exp2_detection.json"):
    global INTERFERER
    INTERFERER = interferer
    jobs = []
    for snr in SNRS:
        for fault in FAULTS + ["HEALTHY"]:
            n = N_PER_CELL * (len(FAULTS) if fault == "HEALTHY" else 1)
            jobs += [(snr, fault, i) for i in range(n)]
    with Pool(2) as p:
        recs = p.map(one_record, jobs, chunksize=8)

    out = {"config": dict(fs=FS, duration=DURATION, n_harm=N_HARM,
                          halfwidth_hz=HALFWIDTH_HZ, n_per_cell=N_PER_CELL,
                          fixed_band=FIXED_BAND, snrs=SNRS,
                          interferer=interferer),
           "cells": {}}
    for snr in SNRS:
        healthy = [r for r in recs if r["snr"] == snr and r["fault"] == "HEALTHY"]
        for st in STRATEGIES:
            neg = np.array([r[st] for r in healthy])
            thr = float(np.quantile(neg, 0.95))     # 5% false alarms
            for fault in FAULTS:
                pos = np.array([r[st] for r in recs
                                if r["snr"] == snr and r["fault"] == fault])
                out["cells"][f"{snr}|{st}|{fault}"] = {
                    "threshold_db": thr,
                    "recall_at_5pct_fa": float((pos > thr).mean()),
                    "auc": auc(pos, neg),
                    "median_stat_db": float(np.median(pos)),
                    "median_healthy_db": float(np.median(neg)),
                }
    # how well the kurtogram finds the resonance
    err = [abs(r["kurtogram_fc"] - r["f_n"]) for r in recs
           if r["fault"] != "HEALTHY" and np.isfinite(r.get("kurtogram_fc", np.nan))]
    out["kurtogram_centre_error_hz"] = {
        "median": float(np.median(err)), "p90": float(np.quantile(err, 0.9))}
    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=1)

    hdr = "  ".join(f"{s:>10s}" for s in STRATEGIES)
    print(f"recall at 5% false alarms\n{'SNR':>5} {'fault':>6}  {hdr}")
    for snr in SNRS:
        for fault in FAULTS:
            row = "  ".join(
                f"{out['cells'][f'{snr}|{st}|{fault}']['recall_at_5pct_fa']:10.2f}"
                for st in STRATEGIES)
            print(f"{snr:5d} {fault:>6}  {row}")
    print("\nkurtogram centre error vs true resonance: "
          f"median {out['kurtogram_centre_error_hz']['median']:.0f} Hz, "
          f"p90 {out['kurtogram_centre_error_hz']['p90']:.0f} Hz")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "clean"
    if which == "clean":
        main()
    else:
        # a strong band-limited interferer over the lower half of the fixed
        # band: cavitation or gear mesh from a neighbouring machine
        print("with band-limited interferer 2-8 kHz at +10 dB")
        main(interferer=((2000.0, 8000.0), 10.0),
             out_path="results/exp2_interferer.json")
