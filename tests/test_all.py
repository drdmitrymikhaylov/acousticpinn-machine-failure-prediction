"""Checks that the physics is right before any of the results are believed."""

import sys

sys.path.insert(0, 'src')

import numpy as np

from bearing_sim import Resonance, simulate, simulate_healthy
from envelope import envelope_spectrum, harmonic_statistic, spectral_kurtosis
from kinematics import SKF_6205
try:                                    # the discharge model is in the private repository
    from pd_features import features
    from pd_sim import simulate_looseness, simulate_pd
except ImportError:                     # public checkout: the two discharge checks are skipped
    features = None

FS = 48000.0


def test_published_multipliers():
    """The kinematic formulas must reproduce the tabulated multipliers."""
    published = {"BPFO": 3.5848, "BPFI": 5.4152, "FTF": 0.3983}
    o = SKF_6205.orders()
    for k, v in published.items():
        assert abs(o[k] - v) < 1e-4, (k, o[k], v)
    # the table's 4.7135 is the ball *impact* rate, twice the ball spin rate
    assert abs(o["BALL_IMPACT"] - 4.7135) < 1e-4, o["BALL_IMPACT"]
    assert abs(2 * o["BSF"] - o["BALL_IMPACT"]) < 1e-12


def test_frequency_sum_rule():
    """BPFO + BPFI = n * fr exactly, for any geometry."""
    o = SKF_6205.orders()
    assert abs(o["BPFO"] + o["BPFI"] - SKF_6205.n_elements) < 1e-12
    # the cage carries the rolling elements: BPFO = n * FTF
    assert abs(o["BPFO"] - SKF_6205.n_elements * o["FTF"]) < 1e-12


def test_kurtosis_zero_for_gaussian():
    rng = np.random.default_rng(0)
    z = rng.normal(size=200_000) + 1j * rng.normal(size=200_000)
    assert abs(spectral_kurtosis(z)) < 0.05


def test_jitter_law_recovers_slip():
    """The measured harmonic decay must return the slip that produced it."""
    fr = 30.0
    fd = SKF_6205.orders()["BPFO"] * fr
    band = (6000.0, 10000.0)
    k = np.arange(1, 21)

    def amps(slip):
        out = []
        for seed in range(4):
            x, _ = simulate(fd, fr, fs=FS, duration=8.0, snr_db=40.0,
                            slip=slip, seed=seed, shaft_harmonics=(),
                            shaft_gain_db=0.0)
            f, E = envelope_spectrum(x, FS, band=band)
            out.append([E[(f >= i * fd - 6) & (f <= i * fd + 6)].max()
                        for i in k])
        return np.mean(out, axis=0)

    ref = amps(0.0)
    for true_slip in (0.005, 0.02):
        a = amps(true_slip) / ref
        use = a > 0.15
        slope = np.polyfit(k[use] ** 2, np.log(a[use]), 1)[0]
        est = np.sqrt(-slope / (2 * np.pi ** 2))
        assert abs(est - true_slip) / true_slip < 0.15, (true_slip, est)


def test_statistic_separates_fault_from_healthy():
    fr = 30.0
    fd = SKF_6205.orders()["BPFO"] * fr
    res = Resonance(f_n=8000.0, zeta=0.02)
    band = (2000.0, 20000.0)

    def stat(x):
        f, E = envelope_spectrum(x, FS, band=band)
        return harmonic_statistic(f, E, fd, 6, 2.0)[0]

    fault = [stat(simulate(fd, fr, fs=FS, duration=1.0, snr_db=-6.0,
                           resonance=res, seed=s)[0]) for s in range(8)]
    ok = [stat(simulate_healthy(fr, FS, 1.0, -6.0, seed=100 + s)[0])
          for s in range(8)]
    assert min(fault) > max(ok) + 5.0, (min(fault), max(ok))


def test_healthy_statistic_is_flat_in_snr():
    """The statistic must not drift with background level on healthy records,
    or the threshold would have to be re-derived for every operating point."""
    fr = 30.0
    fd = SKF_6205.orders()["BPFO"] * fr
    vals = []
    for snr in (0.0, -10.0, -20.0):
        s = [harmonic_statistic(
            *envelope_spectrum(simulate_healthy(fr, FS, 1.0, snr, seed=s)[0],
                               FS, band=(2000.0, 20000.0)),
            fd, 6, 2.0)[0] for s in range(8)]
        vals.append(np.mean(s))
    assert max(vals) - min(vals) < 1.5, vals


def test_pd_more_impulsive_than_hum():
    if features is None:
        return
    fs = 500_000.0
    pd_ = [features(simulate_pd(0.4, fs, 0.0, seed=s)[0], fs)["crest_hf"]
           for s in range(8)]
    hum = [features(simulate_looseness(0.4, fs, 0.0, rattle_prob=0.0,
                                       seed=s)[0], fs)["crest_hf"]
           for s in range(8)]
    assert np.median(pd_) > np.median(hum)


def test_pd_counts_are_stochastic_and_hum_is_not():
    if features is None:
        return
    fs = 500_000.0
    pd_ = np.median([features(simulate_pd(1.0, fs, 6.0, seed=s)[0],
                              fs)["count_fano"] for s in range(6)])
    hum = np.median([features(simulate_looseness(1.0, fs, 6.0, rattle_prob=0.0,
                                                 seed=s)[0],
                              fs)["count_fano"] for s in range(6)])
    assert pd_ > 1.0 > hum, (pd_, hum)


def test_impact_windows_are_ringing_only():
    """The windows must contain the resonance and not the shaft harmonics:
    after the high-pass, the strongest spectral line in a window sits at the
    resonance, and there is no offset."""
    from exp5_resonance_pinn import impact_windows
    fr = 30.0
    fd = SKF_6205.frequencies(fr)["BPFO"]
    x, _ = simulate(fd, fr, fs=FS, duration=1.0, snr_db=0.0, slip=0.01,
                    resonance=Resonance(f_n=7000.0, zeta=0.03), seed=1)
    W = impact_windows(x, FS)
    assert W.shape[0] == 16
    assert abs(W.mean()) < 0.05 * np.abs(W).max()
    S = np.abs(np.fft.rfft(W * np.hanning(W.shape[1]), n=8192, axis=1)).mean(0)
    f = np.fft.rfftfreq(8192, 1 / FS)
    assert abs(f[np.argmax(S)] - 7000.0) < 300.0


def test_joint_fit_recovers_the_resonance():
    """Sixteen noisy ringing windows with one (w_n, zeta): the joint
    free-decay fit returns them to 1% and 20%."""
    from exp5_resonance_pinn import joint_fit, WIN_S
    rng = np.random.default_rng(0)
    t = np.arange(int(WIN_S * FS)) / FS
    wn, z = 2 * np.pi * 9100.0, 0.025
    wd = wn * np.sqrt(1 - z ** 2)
    W = np.stack([np.exp(-z * wn * t) * np.cos(wd * t + rng.uniform(0, 2 * np.pi))
                  + 0.1 * rng.standard_normal(t.size) for _ in range(16)])
    fn, zeta = joint_fit(W, FS)
    assert abs(fn - 9100.0) < 91.0, fn
    assert abs(zeta - z) / z < 0.2, zeta


def test_network_features_stay_below_nyquist():
    """The Fourier features of the ringing network reach every resonance the
    benchmark draws (up to 14 kHz) and stop below the 24 kHz Nyquist rate;
    above it the network could bend between samples and satisfy the equation
    of motion at any w_n."""
    from exp5_resonance_pinn import Ringing, WIN_S
    m = Ringing(1)
    f_max = float(m.k.max()) / np.pi / (2 * WIN_S)
    assert 14000.0 < f_max < FS / 2, f_max


def test_resonance_results_are_what_the_page_says():
    import json, os
    p = "results/exp5_resonance_pinn.json"
    if not os.path.exists(p):
        return
    d = json.load(open(p))
    for s in d["config"]["snrs"]:
        e = d["resonance_error"][str(s)]
        # the joint free-decay fit is never worse than the kurtogram at the median
        assert e["joint"]["median_abs_err_hz"] <= e["kurtogram"]["median_abs_err_hz"] + 1.0
        # and the network does not beat the estimator without a network
        assert e["pinn"]["median_abs_err_hz"] >= 0.5 * e["joint"]["median_abs_err_hz"]


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS  {fn.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL  {fn.__name__}: {e}")
    print(f"\n{len(fns) - failed}/{len(fns)} passed")
    sys.exit(1 if failed else 0)
