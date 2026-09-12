"""Identifying the resonance with the equation of motion, not a kurtogram.

Every envelope method has to choose a demodulation band, and exp2 showed
that the choice is most of the difference between strategies: the kurtogram
buys about 6 dB over a fixed band, the oracle (the true resonance) a few dB
more.  The oracle is not available in the field because nobody measures the
structural resonance of every bearing housing.  But the resonance is in the
signal: each impact rings the structure, and the ringing between impacts
obeys the free-decay equation of a damped oscillator,

    x'' + 2 zeta w_n x' + w_n^2 x = 0 .

This file identifies (w_n, zeta) from the record itself with a physics-
informed network: short windows after the strongest impacts are fitted by a
network x(t) whose residual against the free-decay equation is part of the
loss, with w_n and zeta as trainable parameters shared across all windows of
the record.  The band for demodulation is then centred on the identified
resonance, and the strategy is run through the exp2 detection benchmark
against raw / fixed / kurtogram / oracle.

Two baselines that use no network are kept, because the question is whether
the network buys anything:

  damped-sine    nonlinear least squares of A e^{-zeta w t} sin(w_d t + phi)
                 on the same windows (the classical estimate)
  kurtogram      the band chosen by envelope kurtosis (exp2's adaptive choice)

Records are simulated as in exp2: every machine draws its own shaft speed,
resonance, damping and slip.
"""

import json
import sys
import time

sys.path.insert(0, 'src')

import numpy as np
import torch
from scipy.optimize import least_squares
from scipy.signal import butter, filtfilt

from bearing_sim import Resonance, simulate, simulate_healthy
from envelope import analytic_band, envelope_spectrum, harmonic_statistic, kurtogram
from exp2_detection import (DURATION, FS, HALFWIDTH_HZ, N_HARM, defect_frequency,
                            draw_machine)

torch.set_default_dtype(torch.float64)
SNRS = [-6, -9, -12, -15]
N_MACHINES = 12
WIN_S = 2.5e-3          # window after each impact
N_WIN = 16              # impacts used per record
SKIP_S = 0.15e-3        # first part of the window is the impact itself, not free decay
W_REF = 2 * np.pi * 8000.0


HP_HZ = 1000.0          # the shaft harmonics (+12 dB, tens of Hz) would otherwise
                        # appear inside a 2.5 ms window as an offset and a slope


def impact_windows(x, fs, n_win=N_WIN, win_s=WIN_S):
    """Windows after the strongest peaks of the wide-band envelope, cut from
    the record high-passed above HP_HZ so that only the ringing is in them."""
    b, a = butter(4, HP_HZ / (fs / 2), "high")
    x = filtfilt(b, a, x)
    z = analytic_band(x, fs, 2000.0, fs / 2 - 1)
    env = np.abs(z)
    L = int(win_s * fs)
    guard = int(0.02 * fs)          # keep impacts at least 20 ms apart
    order = np.argsort(env)[::-1]
    starts, taken = [], np.zeros(len(x), bool)
    for i in order:
        if i + L >= len(x) or taken[max(i - guard, 0): i + guard].any():
            continue
        starts.append(i); taken[max(i - guard, 0): i + guard] = True
        if len(starts) >= n_win:
            break
    return np.array([x[s: s + L] for s in starts])


def damped_sine_fit(W, fs, f0=8000.0):
    """Nonlinear least squares per window, median over windows."""
    t = np.arange(W.shape[1]) / fs
    est = []
    for w in W:
        w = w / (np.abs(w).max() + 1e-12)
        def resid(p):
            A, lam, fd, ph = p
            return A * np.exp(-lam * t) * np.sin(2 * np.pi * fd * t + ph) - w
        best = None
        for f_start in (f0 * 0.6, f0, f0 * 1.5):
            s = least_squares(resid, [1.0, 1000.0, f_start, 0.0],
                              bounds=([0, 0, 500, -np.pi], [10, 5e4, fs / 2, np.pi]))
            if best is None or s.cost < best.cost:
                best = s
        A, lam, fd, ph = best.x
        w_n = np.sqrt((2 * np.pi * fd) ** 2 + lam ** 2)
        est.append((w_n / (2 * np.pi), lam / w_n))
    est = np.array(est)
    return float(np.median(est[:, 0])), float(np.median(est[:, 1]))


def joint_fit(W, fs, f_lo=3000.0, f_hi=15000.0):
    """One (w_n, zeta) for all windows, amplitudes and phases per window solved
    linearly (variable projection): coarse grid, then local refinement.  This
    is the exact-physics estimator the network is measured against."""
    t = np.arange(W.shape[1]) / fs
    Wn = W / (np.abs(W).max(axis=1, keepdims=True) + 1e-12)
    keep = t > SKIP_S
    Y = Wn[:, keep].T                    # (L', K)
    tk = t[keep]

    def cost(p):
        f_n, zeta = p
        wn = 2 * np.pi * f_n
        wd = wn * np.sqrt(max(1 - zeta ** 2, 1e-9))
        env = np.exp(-zeta * wn * tk)
        A = np.stack([env * np.cos(wd * tk), env * np.sin(wd * tk)], 1)
        coef, *_ = np.linalg.lstsq(A, Y, rcond=None)
        return float(np.sum((A @ coef - Y) ** 2))

    fs_grid = np.arange(f_lo, f_hi + 1, 50.0)
    z_grid = np.geomspace(0.005, 0.1, 12)
    best = min(((cost((f, z)), f, z) for f in fs_grid for z in z_grid))
    r = least_squares(lambda p: np.sqrt(cost(p)), [best[1], best[2]],
                      bounds=([f_lo, 0.002], [f_hi, 0.3]))
    return float(r.x[0]), float(r.x[1])


class Ringing(torch.nn.Module):
    """x_k(t) for K windows from one network with a window embedding; Fourier
    features so that a 2.5 ms window with ~20 cycles is representable.  With
    sin(pi n t/T) the highest feature is n/(2T) Hz: n_freq = 80 reaches
    16 kHz, above any resonance the benchmark draws (4-14 kHz) and below the
    24 kHz Nyquist frequency of the record.  Both limits matter.  The first
    version stopped at 9.6 kHz and could not represent the upper half; a
    version that reached 25.6 kHz could oscillate *between* the samples and
    satisfy the equation of motion at the samples for any w_n at all."""

    def __init__(self, n_win, n_freq=80):
        super().__init__()
        self.emb = torch.nn.Embedding(n_win, 16)
        self.k = torch.nn.Parameter(torch.arange(1, n_freq + 1, dtype=torch.float64) * np.pi,
                                    requires_grad=False)
        self.net = torch.nn.Sequential(torch.nn.Linear(2 * n_freq + 16, 96), torch.nn.Tanh(),
                                       torch.nn.Linear(96, 96), torch.nn.Tanh(),
                                       torch.nn.Linear(96, 1))
        self.log_wn = torch.nn.Parameter(torch.tensor(np.log(2 * np.pi * 8000.0)))
        self.logit_zeta = torch.nn.Parameter(torch.tensor(np.log(0.03 / 0.97)))

    def forward(self, t_s, idx):
        ff = torch.cat([torch.sin(t_s * self.k), torch.cos(t_s * self.k)], 1)
        return self.net(torch.cat([ff, self.emb(idx)], 1))

    @property
    def wn(self): return torch.exp(self.log_wn)

    @property
    def zeta(self): return torch.sigmoid(self.logit_zeta)


def pinn_fit(W, fs, steps=800, seed=0, lam=50.0):
    torch.manual_seed(seed)
    K, L = W.shape
    scale = np.abs(W).max(axis=1, keepdims=True) + 1e-12
    Wn = torch.tensor(W / scale)
    t = torch.arange(L, dtype=torch.float64) / fs
    T = float(t[-1])
    t_s = (t / T).reshape(-1, 1)
    m = Ringing(K)
    idx_all = torch.arange(K).repeat_interleave(L)
    t_all = t_s.repeat(K, 1)
    y_all = Wn.reshape(-1, 1)
    # the equation of motion is enforced at collocation points drawn between
    # the samples, not at the samples: a network checked only where the data
    # are can bend between them and satisfy any w_n (see Ringing).
    g = torch.Generator().manual_seed(seed)
    n_c = 4 * L
    t_c = (SKIP_S / T + (1 - SKIP_S / T) * torch.rand(K * n_c, 1, generator=g,
                                                       dtype=torch.float64))
    t_c.requires_grad_(True)
    idx_c = torch.arange(K).repeat_interleave(n_c)

    def losses():
        x = m(t_all, idx_all)
        l_data = torch.mean((x - y_all) ** 2)
        xc = m(t_c, idx_c)
        dx = torch.autograd.grad(xc, t_c, torch.ones_like(xc), create_graph=True)[0] / T
        ddx = torch.autograd.grad(dx, t_c, torch.ones_like(dx), create_graph=True)[0] / T
        wn = m.wn
        res = ddx + 2 * m.zeta * wn * dx + wn ** 2 * xc
        # normalised by a FIXED frequency scale.  Dividing by the trainable
        # wn**4 made "small residual" cheapest at large wn and biased every
        # estimate toward 12 kHz (first run: median error 3 kHz).
        l_phys = torch.mean(res ** 2) / W_REF ** 4
        return l_data, l_phys

    opt = torch.optim.Adam(m.parameters(), lr=3e-3)
    for _ in range(steps):
        opt.zero_grad()
        a, b = losses()
        (a + lam * b).backward()
        opt.step()
    lb = torch.optim.LBFGS(m.parameters(), max_iter=100, line_search_fn="strong_wolfe")

    def closure():
        lb.zero_grad()
        a, b = losses()
        l = a + lam * b
        l.backward()
        return l
    lb.step(closure)
    a, b = losses()
    with torch.no_grad():
        return float(m.wn / (2 * np.pi)), float(m.zeta), float(a), float(b)


def band_from_resonance(f_n, f_defect):
    half = max(0.75 * N_HARM * f_defect, 1500.0)
    return (max(f_n - half, 200.0), min(f_n + half, FS / 2 - 1))


def stat_in_band(x, f_defect, band):
    f, E = envelope_spectrum(x, FS, band=band)
    s, _, _ = harmonic_statistic(f, E, f_defect, N_HARM, HALFWIDTH_HZ)
    return s


def auc(pos, neg):
    pos, neg = np.asarray(pos), np.asarray(neg)
    allv = np.concatenate([pos, neg])
    r = np.argsort(np.argsort(allv)) + 1.0
    return float((r[: len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def main():
    t0 = time.time()
    out = {"config": {"snrs": SNRS, "n_machines": N_MACHINES, "win_s": WIN_S, "n_win": N_WIN},
           "cells": {}, "resonance_error": {}}
    for snr in SNRS:
        stats = {k: {"pos": [], "neg": []} for k in ("fixed", "kurtogram", "oracle",
                                                       "damped_sine", "joint", "pinn")}
        errs = {"kurtogram": [], "damped_sine": [], "joint": [], "pinn": []}
        zerr = {"damped_sine": [], "joint": [], "pinn": []}
        for i in range(N_MACHINES):
            rng = np.random.default_rng(1000 * (snr + 100) + i)
            mch = draw_machine(rng)
            seed = int(rng.integers(1 << 30))
            f_defect, _, _ = defect_frequency("BPFO", mch["fr"])
            bw_min = N_HARM * f_defect
            for label in ("pos", "neg"):
                if label == "pos":
                    x, _ = simulate(f_defect, mch["fr"], fs=FS, duration=DURATION, snr_db=snr,
                                    slip=mch["slip"], slip_model=mch["slip_model"],
                                    resonance=Resonance(f_n=mch["f_n"], zeta=mch["zeta"]),
                                    seed=seed)
                else:
                    x, _ = simulate_healthy(mch["fr"], FS, DURATION, snr, seed=seed)
                W = impact_windows(x, FS)
                fn_ds, z_ds = damped_sine_fit(W, FS)
                fn_j, z_j = joint_fit(W, FS)
                fn_pinn, z_pinn, ld, lp = pinn_fit(W, FS, seed=i)
                if label == "pos":
                    shown = (fn_ds, z_ds, fn_j, z_j, fn_pinn, z_pinn, ld, lp)
                _, best = kurtogram(x, FS, n_levels=6, f_min=500.0, bw_min=bw_min)
                bands = {"fixed": (2000.0, 20000.0),
                         "kurtogram": (best[1] - best[2] / 2, best[1] + best[2] / 2),
                         "oracle": band_from_resonance(mch["f_n"], f_defect),
                         "damped_sine": band_from_resonance(fn_ds, f_defect),
                         "joint": band_from_resonance(fn_j, f_defect),
                         "pinn": band_from_resonance(fn_pinn, f_defect)}
                for k, b in bands.items():
                    stats[k][label].append(stat_in_band(x, f_defect, b))
                if label == "pos":
                    errs["kurtogram"].append(abs(best[1] - mch["f_n"]))
                    errs["damped_sine"].append(abs(fn_ds - mch["f_n"]))
                    errs["joint"].append(abs(fn_j - mch["f_n"]))
                    errs["pinn"].append(abs(fn_pinn - mch["f_n"]))
                    zerr["damped_sine"].append(abs(z_ds - mch["zeta"]) / mch["zeta"])
                    zerr["joint"].append(abs(z_j - mch["zeta"]) / mch["zeta"])
                    zerr["pinn"].append(abs(z_pinn - mch["zeta"]) / mch["zeta"])
            fn_ds, z_ds, fn_j, z_j, fn_pinn, z_pinn, ld, lp = shown      # the faulty record
            print(f"  snr {snr} machine {i}: f_n {mch['f_n']:.0f} zeta {mch['zeta']:.3f} | "
                  f"damped-sine {fn_ds:.0f}/{z_ds:.3f}  joint {fn_j:.0f}/{z_j:.3f}  "
                  f"pinn {fn_pinn:.0f}/{z_pinn:.3f} "
                  f"(losses {ld:.2e} {lp:.2e})  kurt {best[1]:.0f}   [{time.time()-t0:.0f} s]")
        cell = {}
        for k, v in stats.items():
            thr = np.quantile(v["neg"], 0.95)
            cell[k] = {"auc": auc(v["pos"], v["neg"]),
                       "detection_at_5pct_fa": float(np.mean(np.array(v["pos"]) > thr))}
        out["cells"][str(snr)] = cell
        out["resonance_error"][str(snr)] = {
            k: {"median_abs_err_hz": float(np.median(v)), "p90_abs_err_hz": float(np.quantile(v, 0.9))}
            for k, v in errs.items()}
        out["resonance_error"][str(snr)]["zeta_rel_err"] = {
            k: float(np.median(v)) for k, v in zerr.items()}
        print(f"SNR {snr:4d} dB  " + "  ".join(f"{k} AUC {c['auc']:.2f}" for k, c in cell.items()))
        print(f"           f_n |err| median Hz: " + ", ".join(
            f"{k} {v['median_abs_err_hz']:.0f}" for k, v in out["resonance_error"][str(snr)].items()
            if k != "zeta_rel_err"))
    with open("results/exp5_resonance_pinn.json", "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
