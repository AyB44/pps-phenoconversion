"""
pps/score.py

General mathematical structure

    L_i(t) = theta_0 + sum_k theta_k * Z_k(t)      (Eq. 1)
    PPS_i(t) = 1 / (1 + exp(-L_i(t)))              (Eq. 2)

L is an unbounded, unitless latent score. PPS is the bounded (0, 1)
activity index actually used for clinical interpretation or downstream
PK/PD simulation -- negative or positive L values are intermediate
outputs, never PPS values themselves.
"""
from dataclasses import dataclass, field
from typing import Dict, Optional
import math

# Canonical module names, per Section 4.
GERMLINE_MODULES = ("G1", "G2", "G3")
DYNAMIC_MODULES = ("XEN", "MIC", "IMM", "TX", "PT", "HEP", "REN")
ALL_MODULES = GERMLINE_MODULES + DYNAMIC_MODULES


@dataclass
class PPSCoefficients:
    """theta_k in Eq. 1.

    Any module not given an explicit weight defaults to 1.0, matching
    the "unit coefficients used only for illustration" convention in
    Section 8's worked examples. Real applications should replace
    these with coefficients calibrated per Section 7 (PK/PD endpoints,
    regularized regression, Bayesian hierarchical models, etc.).
    """
    intercept: float = 0.0
    weights: Dict[str, float] = field(default_factory=dict)

    def get(self, module: str) -> float:
        return self.weights.get(module, 1.0)


def latent_score(z: Dict[str, float], theta: Optional[PPSCoefficients] = None) -> float:
    """Eq. 1: the unbounded latent composite score L_i(t).

    Parameters
    ----------
    z : mapping of module name -> standardized module score Z_k(t).
        Only modules actually entered into the model for this
        patient/timepoint should be present. Modules excluded under
        Section 6 (to avoid double counting) or simply not measured
        are omitted here, not set to 0 -- see pps.modules.DynamicState.
    theta : PPSCoefficients, optional. Defaults to intercept 0 and a
        unit weight on every module supplied in z.
    """
    if theta is None:
        theta = PPSCoefficients()
    unknown = set(z) - set(ALL_MODULES)
    if unknown:
        raise ValueError(f"Unrecognized module name(s): {sorted(unknown)}")
    return theta.intercept + sum(theta.get(k) * v for k, v in z.items())


def pps_from_latent(L: float) -> float:
    """Eq. 2: the monotonic logistic mapping L(t) -> PPS(t) in (0, 1)."""
    return 1.0 / (1.0 + math.exp(-L))


def compute_pps(z: Dict[str, float], theta: Optional[PPSCoefficients] = None) -> float:
    """Convenience wrapper: standardized module scores -> PPS (Sections 2 & 8)."""
    return pps_from_latent(latent_score(z, theta))