"""
PPS is presented in the source manuscript as a hypothesis-generating,
not-yet-prospectively-validated modeling framework. This package
implements the mathematics exactly as specified, and is
not itself a validated clinical prediction tool.
"""
from .score import PPSCoefficients, latent_score, pps_from_latent, compute_pps
from .modules import GermlineProfile, DynamicState, combine
from .standardize import Standardizer, code_categorical
from .residualize import fit_conditional_mean, residualize, residualize_single
from .temporal import update_pps, PPSTracker
from .pk import MichaelisMentenParams, mm_rate

__all__ = [
    "PPSCoefficients", "latent_score", "pps_from_latent", "compute_pps",
    "GermlineProfile", "DynamicState", "combine",
    "Standardizer", "code_categorical",
    "fit_conditional_mean", "residualize", "residualize_single",
    "update_pps", "PPSTracker",
    "MichaelisMentenParams", "mm_rate",
]

__version__ = "0.1.0"