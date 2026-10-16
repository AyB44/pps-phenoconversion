"""
pps -- reference implementation of the Phenotype-Predictive Score (PPS)
framework (Supplementary Methods 2, Porrogi 2026).

Eq. 2 is PPS = 2 * sigmoid(L): PPS = 1.0 at the neutral latent score,
values below 1 mean reduced relative capacity, values above 1 mean
increased relative capacity, and the transform is bounded in (0, 2).
PPS is a relative functional-capacity index, not a probability.

The framework is hypothesis-generating and not a validated clinical
prediction model; it must not be used to determine drug doses.
"""
from .score import (
    PPSCoefficients, latent_score, pps_from_latent, compute_pps,
    PPS_CEILING_DEFAULT,
    pps_from_latent_logistic_deprecated, compute_pps_logistic_deprecated,
)
from .modules import GermlineProfile, DynamicState, combine
from .standardize import Standardizer, code_categorical
from .residualize import fit_conditional_mean, residualize, residualize_single
from .temporal import update_pps, PPSTracker
from .pk import MichaelisMentenParams, mm_rate

__all__ = [
    "PPSCoefficients", "latent_score", "pps_from_latent", "compute_pps",
    "PPS_CEILING_DEFAULT",
    "pps_from_latent_logistic_deprecated", "compute_pps_logistic_deprecated",
    "GermlineProfile", "DynamicState", "combine",
    "Standardizer", "code_categorical",
    "fit_conditional_mean", "residualize", "residualize_single",
    "update_pps", "PPSTracker",
    "MichaelisMentenParams", "mm_rate",
]

__version__ = "0.2.0"