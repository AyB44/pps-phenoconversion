"""
pps/pk.py

Use of PPS in pharmacokinetic/pharmacodynamic simulations:

    Vmax_eff(t) = PPS_i(t) * Vmax_ref
    Km_eff      = Km_ref                (held constant)

Reduced PPS is interpreted as reduced effective metabolic capacity
(enzyme abundance/availability), not altered substrate affinity. Km
should only change if independent mechanistic evidence supports an
altered enzyme-substrate interaction -- this module deliberately does
not expose a way to scale Km, so that can't happen silently.
"""
from dataclasses import dataclass


@dataclass
class MichaelisMentenParams:
    Vmax_ref: float
    Km_ref: float

    def scaled(self, pps: float) -> "MichaelisMentenParams":
        """Apply PPS as a Vmax scalar; Km held constant (Section 9)."""
        if not 0.0 <= pps <= 1.0:
            raise ValueError("PPS must lie in [0, 1].")
        return MichaelisMentenParams(Vmax_ref=pps * self.Vmax_ref, Km_ref=self.Km_ref)


def mm_rate(substrate_conc: float, params: MichaelisMentenParams) -> float:
    """Standard Michaelis-Menten rate: v = Vmax * [S] / (Km + [S])."""
    return params.Vmax_ref * substrate_conc / (params.Km_ref + substrate_conc)