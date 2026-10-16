"""
pps/pk.py

Section 9 (Use of PPS in PK/PD simulations) of Supplementary Methods 2,
updated for the corrected PPS range (0 to 2, reference 1.0).

    Vmax_eff = PPS * Vmax_ref
    Km_eff   = Km_ref            (held constant)

PPS = 1.0 preserves the reference Vmax; PPS < 1 gives proportionally
lower effective Vmax; PPS > 1 gives proportionally higher effective
Vmax (induction). This is an illustrative modelling assumption: it does
not establish that every mechanism of phenoconversion acts only through
maximum metabolic capacity. Km is deliberately not scalable here.
"""
from dataclasses import dataclass

from .score import PPS_CEILING_DEFAULT


@dataclass
class MichaelisMentenParams:
    Vmax_ref: float
    Km_ref: float

    def scaled(self, pps: float, pps_max: float = PPS_CEILING_DEFAULT) -> "MichaelisMentenParams":
        """Apply PPS as a Vmax multiplier; Km held constant (Section 9).

        PPS must lie in [0, pps_max] (default 0 to 2). Pass a different
        pps_max only if a non-default ceiling was used in pps.score.
        """
        if not 0.0 <= pps <= pps_max:
            raise ValueError(f"PPS must lie in [0, {pps_max}].")
        return MichaelisMentenParams(Vmax_ref=pps * self.Vmax_ref, Km_ref=self.Km_ref)


def mm_rate(substrate_conc: float, params: MichaelisMentenParams) -> float:
    """Standard Michaelis-Menten rate: v = Vmax * [S] / (Km + [S])."""
    return params.Vmax_ref * substrate_conc / (params.Km_ref + substrate_conc)