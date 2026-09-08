# Vibroacoustic diagnostics of industrial machinery

What a condition-monitoring module can and cannot decide from sound and
vibration — worked out quantitatively, for rotating machinery and for
high-voltage plant.

A predictive-diagnostics module clamped to a pump, a compressor or a
transformer tank has to answer three questions before it is worth its
installation cost:

1. **Which part of the spectrum should it listen to?** The defect is never
   where the energy is.
2. **How much of the diagnostic signature survives the machine's own
   irregularity?** Bearings slip; discharges ignite at random.
3. **When should it raise an alarm?** A system that cries wolf is switched
   off, and a system that warns an hour ahead has warned nobody.

This repository documents the measurements behind those answers.

---

## Read this first: what the numbers are measured on

Every result here comes from **physical signal models**, not from plant
recordings. Each model is written from the mechanism up — impulse trains from
a spalled raceway convolved with a structural resonance, discharge bursts
phase-locked to the mains and filtered by the propagation path — and each is
checked against an independent prediction before it is used (see
*Verification*).

That makes these results statements about **method and detection limits**,
not about any particular machine. Where the model is easier than plant, this
page says so rather than quoting the flattering number. Field validation is
the separate step, and nothing here substitutes for it.

---

## Design rules that came out of it

| # | Rule | Where it comes from |
|---|------|---------------------|
| 1 | Never look for a defect frequency in the raw spectrum. Demodulate. | §1 — detection rate 0.24 vs 1.00 at −12 dB |
| 2 | Sum no more than `k_max ≈ 0.225 / slip` harmonics. At 1% slip that is 22; at 4% it is 6. | §2, verified to 3% over a 16× range of slip |
| 3 | Widen the harmonic window only if the slip is a random walk. For independent jitter a wider window collects noise, not signal. | §2 — linewidth grows as k², or not at all |
| 4 | Do not pay for per-record adaptive band selection. Pay instead for measuring each sensor's structural resonance once, at commissioning. | §3 — the kurtogram never beats a fixed wide band; knowing the resonance beats both by 0.5 in detection rate at −18 dB |
| 5 | An ultrasonic-band energy ratio does not detect partial discharge against sensor noise. It is the phase-locked event structure that does. | §4 — ratio AUC 0.36–0.55 |
| 6 | Require M consecutive exceedances rather than a higher threshold. Three-in-a-row cut false alarms from 315/month to none observed while keeping 29 h of the 48 h warning. | §5 |

---

## 1. The defect is not where the energy is

![demodulation](figures/01_demodulation.png)

A localised raceway defect produces a train of short impacts. Each impact
rings the structure between the defect and the sensor, so the energy lands at
the structural resonance — several kilohertz — and the defect frequency
appears only as the *repetition rate* of that ringing, i.e. as a modulation.

The left panel is the plain spectrum over the frequency range where the
defect frequency lives. The running speed and its harmonics are there; the
defect comb (red lines) is not. The right panel is the spectrum of the
amplitude envelope of the band around the resonance, from the very same
record: the comb at the ball-pass frequency of the outer race is unmistakable,
out to the eighth harmonic.

Over a population of machines with unknown resonances, at 5% false alarms:

| SNR | raw spectrum | envelope, fixed band |
|-----|--------------|----------------------|
| −6 dB | 0.56 | 1.00 |
| −12 dB | 0.24 | 1.00 |
| −15 dB | 0.07 | 0.99 |

The middle panel is the kurtogram — the kurtosis of the band-limited
envelope over a grid of centre frequencies and bandwidths. Impulsive content
raises it, so the brightest cell should be the band the impacts live in. It
finds the true 8 kHz resonance here. Section 3 asks whether that is worth
anything.

**One constraint worth stating explicitly:** a demodulation band of width *B*
produces an envelope whose spectrum only reaches about *B*. Asking for six
harmonics of a 107 Hz defect frequency therefore requires at least ~640 Hz of
bandwidth. Without that constraint the kurtogram of a *healthy* record drifts
into narrow low-frequency bands whose envelope spectrum rolls off inside the
search range, and any statistic measured against a global noise floor then
reports a defect that is not there. This is the kind of bug that ships.

---

## 2. How far up the harmonic series can you go?

![slip](figures/02_slip.png)

Rolling elements do not roll perfectly. The impulse train is not exactly
periodic, and the high harmonics lose coherence first. How fast depends on
what *kind* of irregularity it is, and the two kinds behave differently in
kind, not merely in degree.

**Independent jitter.** If each impact is displaced independently by
`e ~ N(0, σ_t²)`, the expected amplitude at frequency *f* is multiplied by the
characteristic function of the jitter, `exp(−2π²f²σ_t²)`. Writing the slip as
a fraction of a period, `σ = σ_t·f_d`, the k-th harmonic is attenuated by

```
a(k) = exp(−2π² k² σ²)      →      k_max = 0.2251 / σ
```

Measured against that prediction (left panel, points vs lines), fitting the
decay back to the slip that produced it:

| true slip | recovered slip | k_max predicted | k_max measured |
|-----------|----------------|-----------------|----------------|
| 0.25% | 0.27% | 90 | 84 |
| 0.5% | 0.51% | 45 | 44 |
| 1% | 0.99% | 22.5 | 22.7 |
| 2% | 1.98% | 11.3 | 11.4 |
| 4% | 3.90% | 5.6 | 5.8 |

The energy lost from the lines does not disappear — it reappears as a
broadband pedestal, which is the grey band in the figure. Below it the
"harmonic amplitude" being measured is the pedestal, not the harmonic.

*A control that mattered:* the resonance's own impulse response has a finite
decay time (~1 ms here), and that alone rolls the harmonic series off with a
corner near 160 Hz — below the second harmonic. Fitting the jitter law
without dividing that out returns σ ≈ 0.017 for every slip from 0.25% to 4%:
it measures the impulse shape and calls it slip.

**Random-walk slip.** If instead the *interval* fluctuates, phase error
accumulates and the lines broaden rather than fade. Measured exponent of
linewidth against harmonic order: **1.98, 2.01, 2.00** for slips of 0.25%,
0.5% and 1% — the predicted k². (At 2% and 4% the fit falls to 1.64 and 1.12,
because neighbouring harmonics have merged and there is no separate line left
to measure.)

The practical consequence is rule 3 above: the two cases call for opposite
window policies, so a module that widens its harmonic windows "to be safe"
is collecting noise in half of all installations.

---

## 3. Choosing the band: the result that was not expected

![band selection](figures/03_band_selection.png)

Four strategies over a fleet whose resonances are drawn between 4 and 14 kHz
and are *not* known to the detector: no demodulation; one fixed 2–20 kHz band
for everything; the kurtogram, choosing per record; and an oracle centred on
that record's true resonance.

Adaptive band selection **does not beat the fixed wide band**. On a white
background the fixed band is at or above the kurtogram everywhere. With a
strong band-limited interferer over the lower half of the band — cavitation,
gear mesh, a neighbouring machine — the kurtogram pulls ahead only in the
last two SNR points, and not by much.

The oracle, meanwhile, is far ahead of both: at −18 dB it detects 0.97 where
the fixed band manages 0.47 and the kurtogram 0.58. The kurtogram's estimate
of the resonance has a median error of 340 Hz but a 90th percentile of
2.6 kHz, rising to 6.6 kHz once there is interference — it is right most of
the time and badly wrong often enough to matter.

So the money is not in cleverer per-record band selection. It is in **knowing
the resonance**, which is a fixed property of the sensor and its mounting and
can be measured once at commissioning with an impact test, then stored per
sensor. That is a configuration step, not an algorithm.

---

## 4. Partial discharge, and the loose bolt that looks exactly like it

![partial discharge](figures/04_partial_discharge.png)

Discharges are driven by the instantaneous voltage, so they cluster near the
peaks — twice per cycle, ~100 Hz on a 50 Hz supply. A loose clamp on the same
tank does the same thing, at the same rate, for entirely mechanical reasons.
A detector keyed on "activity at twice the mains frequency" answers
*discharge* to both, and every such answer dispatches a crew to healthy plant.

**Counting comes before classifying.** The simulator knows how many bursts it
emitted, so the event detector can be calibrated directly. Adaptive threshold
at median + 8·MAD, 40 mains cycles:

| SNR | detected / emitted |
|-----|--------------------|
| +6 dB | 1.07 |
| 0 dB | 0.98 |
| −6 dB | 0.88 |
| −12 dB | 0.58 |
| −18 dB | 0.14 |

Below about −12 dB the count is a **lower bound, not a measurement**, and
every statistic built on top of it — rate, phase concentration, Fano factor —
inherits that. A module in this regime should report "no decision", not a
confident number. Longer windows do not fix it: 10 cycles and 100 cycles give
the same recovery ratio, because the problem is per-event detectability, not
sample size.

**The single-number rule of thumb fails.** Ultrasonic-band energy over
audio-band energy — what a threshold on a broadband ultrasonic meter amounts
to — separates discharge from sensor noise at AUC 0.36–0.55, i.e. not at all,
and at the low end *backwards*: broadband sensor noise has a flatter spectrum
than a discharge that has crossed oil and steel, so the ratio can be higher
for a healthy sensor than for a live fault. What separates them is the
phase-locked event structure, which no energy ratio can see.

**Honest limit of this section.** Against the mechanical confuser, the six
features separate discharge at AUC ≈ 1.0 at every SNR tested — but so does
the crest factor on its own, and nearly so does the energy ratio. That says
the model's two classes are further apart than plant classes are, not that the
problem is solved. The transferable results here are the counting calibration
and the negative result about energy ratios; the AUC is not.

---

## 5. From a statistic to an alarm somebody will act on

![alarm policy](figures/05_alarm_policy.png)

One measurement every 10 minutes is 4 320 decisions a month. A threshold with
a 1% per-record false-alarm rate therefore fires about 43 times a month, and
the system gets switched off. The experiment runs a fleet through 24 h of
quiet history and then 48 h in which the defect grows, sets the threshold from
the first half of each machine's *own* history, and measures false alarms on
the held-out half:

| threshold | consecutive | false alarms / month | median warning before failure |
|-----------|-------------|----------------------|-------------------------------|
| 95th pct | 1 | 315 | 46.7 h |
| 99th pct | 1 | 113 | 44.3 h |
| 99.9th pct | 1 | 83 | 41.5 h |
| 95th pct | 2 | 22.5 | 34.0 h |
| 99th pct | 2 | 2.5 | 30.8 h |
| 99th pct | 3 | none observed | 27.5 h |

Two things are visible. First, **raising the threshold barely helps**: going
from the 95th to the 99.9th percentile only takes 315 false alarms a month
down to 83, because 72 baseline samples cannot resolve a 1-in-1000 tail — the
quantile is beyond the data. Second, **requiring consecutive exceedances
works**, and cheaply: three in a row removed every false alarm in 1 728
held-out opportunities while still leaving 27 of the 48 hours of warning.

"None observed" is not zero. By the rule of three, 0 events in 1 728 trials
bounds the rate at under 7.5 false alarms per month with 95% confidence. To
claim less than one a month, the baseline has to be roughly ten times longer.

The alarm fires at a median defect level of −18 dB SNR, which is where §3 put
the detection floor for a fixed band. The chain is consistent end to end.

---

## Verification

Physics and estimator checks run as a test suite (8 checks, all passing):

- kinematic frequencies reproduce the tabulated multipliers for a documented
  bearing geometry to 2 × 10⁻⁵ — including resolving an apparent factor-of-two
  discrepancy in the published ball-spin figure, which lists the *impact* rate
  (a spall strikes both races), not the ball rotation rate
- the identities `BPFO + BPFI = n·f_r` and `BPFO = n·FTF` hold exactly
- the spectral-kurtosis estimator returns zero on circular complex Gaussian
  noise
- the jitter law recovers the slip that generated the signal, to 15%
- the harmonic statistic on healthy records does not drift with background
  level (< 1.5 dB across a 20 dB range) — otherwise the threshold would have
  to be re-derived at every operating point
- discharge records are more impulsive, and their per-half-cycle counts more
  variable (Fano > 1 vs < 1), than mechanical ones

---

## What this does not show

- **No plant data.** Detection rates here are properties of the models, not
  of any installed system. Only field data can confirm the SNR at which a
  particular sensor and mounting actually operate.
- **One fault at a time.** Real machines present bearing wear, misalignment
  and looseness simultaneously; the interaction is not modelled.
- **Stationary speed.** Variable-speed drives require order tracking, which is
  a separate problem and is not addressed.
- **Single resonance.** A real transmission path has many, and the impulse
  response is not a single decaying sinusoid.
- **The discharge classes are too far apart**, as §4 says of itself.

---

## Source code

**The source code for this project is not public.** This page documents the
method, the measurements and the conclusions; the implementation is held in a
private repository and is available under NDA.

What is described here: the signal models, the demodulation and event-detection
front ends, the detection statistics, the experiment protocols, and the
figures generated from their outputs.

---

## Scope

The three regimes covered — rotating machinery, partial discharge in
high-voltage plant, and threshold policy for continuous monitoring — are the
ones where the acoustic signature has a mechanism you can write down and
therefore a detection limit you can compute. Boiler and heat-exchanger
monitoring is a fourth such area, where the diagnostic bands are set by phase
transition and combustion rather than by rotation; it is not covered here.

## Licence

Documentation and figures: CC BY 4.0.
