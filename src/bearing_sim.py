"""Vibration model of a localised bearing defect.

Signal model (the standard cyclostationary description of a spalled bearing):

    x(t) = sum_i  A_i * s(t - t_i)  +  shaft harmonics  +  white noise

where

  t_i     impulse arrival times, nominally i/f_defect, perturbed by slip
  A_i     impulse amplitude; constant for an outer-race defect, modulated at
          the shaft rate for an inner-race defect (the fault rotates in and
          out of the load zone)
  s(t)    impulse response of the structure between defect and sensor,
          modelled as a single decaying resonance
          s(t) = exp(-zeta * w_n * t) * sin(w_d * t),  t >= 0

Two slip models are provided, because they behave very differently and the
distinction decides how a diagnostic band should be built:

  "iid"    each impulse is displaced independently, t_i = i/f_d + e_i.
           Spectral lines stay narrow but lose amplitude; the lost energy
           reappears as a broadband pedestal.  Widening the analysis bin does
           not recover it.
  "walk"   the inter-arrival time itself fluctuates, t_i = t_{i-1} + T + e_i.
           The phase error accumulates, so lines broaden with harmonic order.
           Widening the bin does recover the energy.

Amplitudes are dimensionless; scale to engineering units at the sensor.
"""

from dataclasses import dataclass

import numpy as np


@dataclass
class Resonance:
    f_n: float = 8000.0   # undamped natural frequency, Hz
    zeta: float = 0.02    # damping ratio
    dur: float = 0.01     # truncation length of the impulse response, s

    def impulse_response(self, fs: float) -> np.ndarray:
        t = np.arange(0.0, self.dur, 1.0 / fs)
        w_n = 2 * np.pi * self.f_n
        w_d = w_n * np.sqrt(max(1.0 - self.zeta ** 2, 0.0))
        h = np.exp(-self.zeta * w_n * t) * np.sin(w_d * t)
        # unit energy so that "impulse power" is set purely by the amplitudes
        return h / np.linalg.norm(h)


def _arrival_times(f_defect, duration, slip, model, rng):
    """Impulse arrival times over [0, duration) under the chosen slip model."""
    T = 1.0 / f_defect
    n = int(np.floor(duration * f_defect)) + 2
    if model == "iid":
        t = np.arange(n) * T + rng.normal(0.0, slip * T, n)
    elif model == "walk":
        steps = T + rng.normal(0.0, slip * T, n)
        t = np.cumsum(steps) - steps[0]
    else:
        raise ValueError(f"unknown slip model {model!r}")
    t = np.sort(t)
    return t[(t >= 0) & (t < duration)]


def _place(times, amps, fs, n_samples):
    """Sparse impulse train with linear interpolation for sub-sample timing."""
    x = np.zeros(n_samples)
    pos = times * fs
    i0 = np.floor(pos).astype(int)
    frac = pos - i0
    ok = (i0 >= 0) & (i0 + 1 < n_samples)
    np.add.at(x, i0[ok], amps[ok] * (1.0 - frac[ok]))
    np.add.at(x, i0[ok] + 1, amps[ok] * frac[ok])
    return x


def simulate(
    f_defect,
    fr_hz,
    fs=48000.0,
    duration=1.0,
    snr_db=0.0,
    slip=0.01,
    slip_model="iid",
    resonance=Resonance(),
    load_modulation=0.0,
    modulation_hz=None,
    shaft_harmonics=(1.0, 0.5, 0.3),
    shaft_gain_db=12.0,
    seed=0,
):
    """Return (x, meta).

    snr_db          impulse-train power relative to the white-noise power
    load_modulation depth of amplitude modulation (0 = outer race, which is
                    fixed in the load zone; ~0.8 = inner race, which rotates
                    through it)
    modulation_hz   modulation rate; defaults to the shaft rate.  A defect on
                    a rolling element is carried around by the cage, so its
                    modulation rate is the cage rate (FTF), not the shaft rate
    shaft_gain_db   power of the shaft-harmonic interference relative to the
                    impulse train.  Positive values are the realistic case:
                    the raw spectrum is dominated by running speed, which is
                    exactly why the defect must be found by demodulation.
    """
    rng = np.random.default_rng(seed)
    n = int(round(duration * fs))
    t = np.arange(n) / fs

    times = _arrival_times(f_defect, duration, slip, slip_model, rng)
    amps = np.ones_like(times)
    f_mod = fr_hz if modulation_hz is None else modulation_hz
    if load_modulation:
        amps *= 1.0 + load_modulation * np.cos(2 * np.pi * f_mod * times)
    amps = np.maximum(amps, 0.0)

    train = _place(times, amps, fs, n)
    h = resonance.impulse_response(fs)
    impulses = np.convolve(train, h)[:n]

    p_imp = np.mean(impulses ** 2)
    if p_imp <= 0:
        raise RuntimeError("empty impulse train")
    impulses = impulses / np.sqrt(p_imp)          # unit power
    p_imp = 1.0

    shaft = np.zeros(n)
    for k, a in enumerate(shaft_harmonics, start=1):
        shaft += a * np.sin(2 * np.pi * k * fr_hz * t + rng.uniform(0, 2 * np.pi))
    if np.any(shaft):
        shaft *= np.sqrt(10 ** (shaft_gain_db / 10.0) / np.mean(shaft ** 2))

    p_noise = p_imp / (10 ** (snr_db / 10.0))
    noise = rng.normal(0.0, np.sqrt(p_noise), n)

    x = impulses + shaft + noise
    meta = dict(
        fs=fs,
        f_defect=f_defect,
        fr_hz=fr_hz,
        n_impulses=len(times),
        slip=slip,
        slip_model=slip_model,
        snr_db=snr_db,
        f_n=resonance.f_n,
        zeta=resonance.zeta,
        load_modulation=load_modulation,
    )
    return x, meta


def simulate_healthy(fr_hz, fs=48000.0, duration=1.0, snr_db=0.0,
                     shaft_harmonics=(1.0, 0.5, 0.3), shaft_gain_db=12.0,
                     seed=0):
    """Same background, no impulse train.  Power is matched to `simulate` so
    that healthy and faulty records differ only by the defect."""
    rng = np.random.default_rng(seed)
    n = int(round(duration * fs))
    t = np.arange(n) / fs
    shaft = np.zeros(n)
    for k, a in enumerate(shaft_harmonics, start=1):
        shaft += a * np.sin(2 * np.pi * k * fr_hz * t + rng.uniform(0, 2 * np.pi))
    if np.any(shaft):
        shaft *= np.sqrt(10 ** (shaft_gain_db / 10.0) / np.mean(shaft ** 2))
    noise = rng.normal(0.0, np.sqrt(10 ** (-snr_db / 10.0)), n)
    return shaft + noise, dict(fs=fs, fr_hz=fr_hz, snr_db=snr_db, f_defect=None)


def band_limited_noise(n, fs, band, power, rng):
    """Gaussian noise confined to `band`, scaled to the requested power.

    Real machine halls are not white.  Cavitation, turbulence, gear mesh and
    neighbouring plant put strong energy into a limited part of the spectrum,
    and that is the situation in which the choice of demodulation band stops
    being cosmetic.
    """
    w = rng.normal(0.0, 1.0, n)
    W = np.fft.rfft(w)
    f = np.fft.rfftfreq(n, 1.0 / fs)
    W[(f < band[0]) | (f > band[1])] = 0.0
    y = np.fft.irfft(W, n=n)
    p = np.mean(y ** 2)
    if p <= 0:
        return y
    return y * np.sqrt(power / p)
