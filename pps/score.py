"""
pps/score.py

Section 2 (General mathematical structure) of Supplementary Methods 2,
with Eq. 2 corrected so the stated reference value is actually reached.

    L_i(t)   = theta_0 + sum_k theta_k * Z_k(t)         (Eq. 1, unchanged)
    PPS_i(t) = 2 / (1 + exp(-L_i(t))) = 2 * sigmoid(L)  (Eq. 2, corrected)

Why Eq. 2 changed: the original Eq. 2, PPS = sigmoid(L), gave 0.5 -- not
the manuscript's stated reference value of 1.0 -- at the neutral latent
score L = 0, and could never exceed 1.0, so increased capacity
(induction) could not be represented.

Properties of the corrected form:
    L = 0           -> PPS = 1.0 exactly (reference functional capacity)
    L < 0           -> PPS < 1.0 (reduced relative capacity)
    L > 0           -> PPS > 1.0 (increased relative capacity)
    L -> -infinity  -> PPS -> 0
    L -> +infinity  -> PPS -> 2.0 (an assumption of the model, not a
                       biological law: it does NOT claim capacity cannot
                       exceed twice the reference)

PPS is a relative functional-capacity index. It is not a probability.

OPTIONAL `ceiling`: the upper limit defaults to 2.0 (the agreed form).
For any ceiling c > 1 the general form is
    PPS = c * sigmoid(L - ln(c - 1))
which still gives PPS = 1.0 at L = 0 and reduces exactly to
2 * sigmoid(L) when c = 2. It exists for sensitivity analysis only;
the default should be used unless the corresponding author agrees
otherwise.

The original sigmoid(L) functions are kept, marked deprecated, only so
the Supplementary Methods 2 worked examples can still be audited.
"""
from dataclasses import dataclass, field
from typing import Dict, Optional
import math

# Canonical module names, per Section 4. Unaffected by the Eq. 2 correction.
GERMLINE_MODULES = ("G1", "G2", "G3")
DYNAMIC_MODULES = ("XEN", "MIC", "IMM", "TX", "PT", "HEP", "REN")
ALL_MODULES = GERMLINE_MODULES + DYNAMIC_MODULES

PPS_CEILING_DEFAULT = 2.0


@dataclass
class PPSCoefficients:
    """theta_k in Eq. 1. Unaffected by the Eq. 2 correction.

    Unit weights (default 1.0) are illustrative only; real weights must
    be estimated from patient data.
    """
    intercept: float = 0.0
    weights: Dict[str, float] = field(default_factory=dict)

    def get(self, module: str) -> float:
        return self.weights.get(module, 1.0)


def latent_score(z: Dict[str, float], theta: Optional[PPSCoefficients] = None) -> float:
    """Eq. 1: the unbounded latent composite score L_i(t). Unchanged."""
    if theta is None:
        theta = PPSCoefficients()
    unknown = set(z) - set(ALL_MODULES)
    if unknown:
        raise ValueError(f"Unrecognized module name(s): {sorted(unknown)}")
    return theta.intercept + sum(theta.get(k) * v for k, v in z.items())


def _sigmoid(x: float) -> float:
    """Numerically stable logistic function (no overflow for large |x|)."""
    if x >= 0:
        return 1.0 / (1.0 + math.exp(-x))
    e = math.exp(x)
    return e / (1.0 + e)


def pps_from_latent(L: float, ceiling: float = PPS_CEILING_DEFAULT) -> float:
    """Eq. 2 (corrected): PPS = 2 * sigmoid(L) by default.

    Returns exactly 1.0 at L = 0 for any ceiling > 1.
    """
    if ceiling <= 1.0:
        raise ValueError("ceiling must be > 1 so that the reference value 1.0 is attainable.")
    shift = math.log(ceiling - 1.0)  # exactly 0.0 when ceiling == 2
    return ceiling * _sigmoid(L - shift)


def compute_pps(z: Dict[str, float], theta: Optional[PPSCoefficients] = None,
                ceiling: float = PPS_CEILING_DEFAULT) -> float:
    """Standardized module scores -> PPS (Eqs. 1 and corrected 2)."""
    return pps_from_latent(latent_score(z, theta), ceiling)


# ---------------------------------------------------------------------------
# ORIGINAL transform -- DEPRECATED, audit trail only
# ---------------------------------------------------------------------------

def pps_from_latent_logistic_deprecated(L: float) -> float:
    """ORIGINAL Eq. 2: PPS = sigmoid(L). Gives 0.5 at L = 0 and cannot
    exceed 1.0. Kept only to audit the original worked examples.
    Do not use for new work.
    """
    return _sigmoid(L)


def compute_pps_logistic_deprecated(z: Dict[str, float],
                                    theta: Optional[PPSCoefficients] = None) -> float:
    """DEPRECATED wrapper around the original sigmoid(L) transform."""
    return pps_from_latent_logistic_deprecated(latent_score(z, theta))