"""Envelope demodulation and the fast kurtogram.

A localised bearing defect is not visible in the raw spectrum: the running
speed and its harmonics carry far more energy, and the defect appears only as
a low-level modulation of a high-frequency structural resonance.  The standard
remedy is to band-pass the signal around that resonance, take the amplitude
envelope, and look for the defect frequency in the spectrum of the envelope.

The open question in practice is *which* band.  The kurtogram answers it
without knowing the resonance: impulsive content raises the kurtosis of the
band-limited envelope, so the band that maximises kurtosis is the band the
impulses live in.
"""

import numpy as np


def analytic_band(x, fs, f_lo, f_hi):
    """Analytic signal restricted to [f_lo, f_hi] by frequency-domain masking."""
    n = len(x)
    X = np.fft.rfft(x)
    f = np.fft.rfftfreq(n, 1.0 / fs)
    mask = (f >= f_lo) & (f <= f_hi)
    Y = np.zeros(n, dtype=complex)
    # one-sided spectrum doubled -> analytic signal
    Y[: len(X)] = np.where(mask, X, 0.0) * 2.0
    if n % 2 == 0:
        Y[len(X) - 1] /= 2.0
    Y[0] = 0.0
    return np.fft.ifft(Y)


def spectral_kurtosis(z):
    """Kurtosis of a complex band-limited signal.

    For a circular complex Gaussian process E|z|^4 = 2 (E|z|^2)^2, so the
    estimator below is zero for noise and positive for impulsive content.
    """
    p = np.abs(z) ** 2
    m2 = p.mean()
    if m2 <= 0:
        return 0.0
    return float((p ** 2).mean() / m2 ** 2 - 2.0)


def kurtogram(x, fs, n_levels=7, f_min=500.0, f_max=None, bw_min=0.0):
    """Kurtosis over a dyadic grid of bands.

    Returns (table, best) where `table` is a list of (level, fc, bw, kurtosis)
    and `best` is the entry with the highest kurtosis.  Level L splits the
    usable band into 2^L sub-bands.

    `bw_min` rejects bands too narrow to carry the diagnostic information.
    A demodulation band of width B produces an envelope whose spectrum is
    supported up to about B; asking for K harmonics of a defect frequency f_d
    therefore requires B >= K * f_d.  Without this constraint the kurtogram of
    a healthy record wanders into narrow low-frequency bands whose envelope
    spectrum rolls off inside the search range, which biases any statistic
    measured against a global noise floor.
    """
    if f_max is None:
        f_max = fs / 2.0
    table = []
    for level in range(1, n_levels + 1):
        n_bands = 2 ** level
        edges = np.linspace(f_min, f_max, n_bands + 1)
        bw = (f_max - f_min) / n_bands
        if bw < bw_min:
            continue
        for i in range(n_bands):
            lo, hi = edges[i], edges[i + 1]
            k = spectral_kurtosis(analytic_band(x, fs, lo, hi))
            table.append((level, 0.5 * (lo + hi), hi - lo, k))
    if not table:
        raise ValueError("bw_min excludes every band in the kurtogram grid")
    best = max(table, key=lambda r: r[3])
    return table, best


def envelope_spectrum(x, fs, band=None, detrend=True):
    """Amplitude spectrum of the envelope.  Returns (freqs, amplitude)."""
    if band is None:
        z = analytic_band(x, fs, 0.0, fs / 2.0)
    else:
        z = analytic_band(x, fs, band[0], band[1])
    env = np.abs(z)
    if detrend:
        env = env - env.mean()
    w = np.hanning(len(env))
    E = np.abs(np.fft.rfft(env * w)) / (np.sum(w) / 2.0)
    f = np.fft.rfftfreq(len(env), 1.0 / fs)
    return f, E


def local_floor(freqs, amp, f_centre, halfwidth_hz, span_hz):
    """Median amplitude in a neighbourhood of `f_centre`, peak window removed.

    The reference must be local.  The envelope spectrum of a band-limited
    signal rolls off over a scale set by the demodulation bandwidth, so a
    single global median compares harmonics against a floor measured somewhere
    else on that slope and reports a defect that is not there.
    """
    near = (freqs >= f_centre - span_hz) & (freqs <= f_centre + span_hz)
    peak = (freqs >= f_centre - halfwidth_hz) & (freqs <= f_centre + halfwidth_hz)
    sel = near & ~peak & (freqs > 0)
    if sel.sum() < 8:
        return float("nan")
    return float(np.median(amp[sel]))


def harmonic_statistic(freqs, amp, f_defect, n_harm=8, halfwidth_hz=2.0):
    """Detection statistic: harmonic peaks over their local floor, in dB.

    For each harmonic k = 1..n_harm the peak amplitude within
    +/- halfwidth_hz of k*f_defect is compared with the median amplitude in
    the surrounding half-interval (which lies between harmonics, so it is not
    contaminated by them).  The statistic is the quadratic mean of the
    per-harmonic ratios, in dB.

    Returns (statistic_db, peaks, floors).
    """
    peaks = np.zeros(n_harm)
    floors = np.zeros(n_harm)
    span = 0.5 * f_defect
    for k in range(1, n_harm + 1):
        fc = k * f_defect
        sel = (freqs >= fc - halfwidth_hz) & (freqs <= fc + halfwidth_hz)
        peaks[k - 1] = amp[sel].max() if sel.any() else 0.0
        floors[k - 1] = local_floor(freqs, amp, fc, halfwidth_hz, span)
    ok = np.isfinite(floors) & (floors > 0)
    if not ok.any():
        return 0.0, peaks, floors
    ratio = peaks[ok] / floors[ok]
    return float(20.0 * np.log10(np.sqrt(np.mean(ratio ** 2)))), peaks, floors
