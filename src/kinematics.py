"""Characteristic defect frequencies of rolling-element bearings.

All frequencies are returned in Hz for a shaft rate `fr` in Hz.

Geometry convention
-------------------
n      number of rolling elements
d      rolling-element diameter
D      pitch diameter (ball-centre circle)
alpha  contact angle, radians

The four kinematic frequencies assume pure rolling (no slip):

    FTF  = (fr/2) (1 - (d/D) cos a)          cage / train
    BPFO = (n/2) fr (1 - (d/D) cos a)        ball pass, outer race
    BPFI = (n/2) fr (1 + (d/D) cos a)        ball pass, inner race
    BSF  = (D/2d) fr (1 - ((d/D) cos a)^2)   ball spin

A spall on a rolling element strikes both races, so the observed repetition
rate for a ball defect is 2*BSF.  We return BSF itself and expose the impact
rate separately.  This resolves an apparent discrepancy with published tables:
for SKF 6205-2RS JEM the formula above gives BSF = 2.3567 per revolution while
the commonly quoted table value is 4.7135 -- exactly twice, because the table
lists the impact rate, not the ball rotation rate.  See tests/test_kinematics.
"""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Bearing:
    """Rolling-element bearing geometry.

    `d` and `D` must share units; only their ratio is used.
    """

    name: str
    n_elements: int
    d: float
    D: float
    alpha: float = 0.0

    @property
    def ratio(self) -> float:
        """(d/D) cos(alpha) -- the only geometric quantity that matters."""
        return (self.d / self.D) * np.cos(self.alpha)

    # ---- per-shaft-revolution multipliers ------------------------------
    @property
    def ftf_order(self) -> float:
        return 0.5 * (1.0 - self.ratio)

    @property
    def bpfo_order(self) -> float:
        return 0.5 * self.n_elements * (1.0 - self.ratio)

    @property
    def bpfi_order(self) -> float:
        return 0.5 * self.n_elements * (1.0 + self.ratio)

    @property
    def bsf_order(self) -> float:
        return 0.5 * (self.D / self.d) * (1.0 - self.ratio ** 2)

    @property
    def ball_impact_order(self) -> float:
        """Impact rate of a rolling-element spall: it strikes both races."""
        return 2.0 * self.bsf_order

    def orders(self) -> dict:
        return {
            "FTF": self.ftf_order,
            "BPFO": self.bpfo_order,
            "BPFI": self.bpfi_order,
            "BSF": self.bsf_order,
            "BALL_IMPACT": self.ball_impact_order,
        }

    def frequencies(self, fr_hz: float) -> dict:
        """Defect frequencies in Hz at shaft rate `fr_hz`."""
        return {k: v * fr_hz for k, v in self.orders().items()}


# Two documented geometries.  SKF 6205-2RS JEM is the drive-end bearing of the
# most widely published bearing test rig; its multipliers are quoted in the
# literature as BPFO 3.5848, BPFI 5.4152, FTF 0.3983, BSF 4.7135 per shaft
# revolution, which the formulas above must reproduce (see tests/).
SKF_6205 = Bearing(name="SKF 6205-2RS JEM", n_elements=9, d=0.3126, D=1.537)
SKF_6203 = Bearing(name="SKF 6203-2RS JEM", n_elements=9, d=0.2656, D=1.122)


def sidebands(f_defect: float, fr_hz: float, n_harm: int, n_side: int = 0):
    """Harmonics of a defect frequency, optionally with shaft-rate sidebands.

    An inner-race defect rotates through the load zone, so its impulse train is
    amplitude-modulated at the shaft rate and its spectrum carries sidebands at
    k*f_defect +/- m*fr.  An outer-race defect is stationary in the load zone
    and has none.
    """
    out = []
    for k in range(1, n_harm + 1):
        for m in range(-n_side, n_side + 1):
            f = k * f_defect + m * fr_hz
            if f > 0:
                out.append(f)
    return np.array(sorted(set(out)))
